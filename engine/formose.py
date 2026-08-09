"""
The formose reaction: the classic prebiotic route from formaldehyde to
sugars, including ribose (C5H10O5) -- a required building block of RNA,
and arguably a fairer target than hydrocarbons for "does a molecule life
actually needs have a plausible abiotic source."

This is a genuinely different chemistry regime from the rest of the
simulator. Everything in engine/reactions.py is homolytic radical
chemistry (UV/discharge breaks a bond, two unpaired electrons pair back
up). Formose proceeds through base/mineral-catalyzed **aldol addition** --
a closed-shell, ionic mechanism (an enolate carbon attacks a carbonyl
carbon) with no radicals anywhere. Rather than force that into the
radical-species machinery, it gets its own species and reaction rules
here, reusing only the generic `Reaction`/propensity model and the
Gillespie SSA algorithm itself.

**The real difficulty (the "sugar problem"):** formose chemistry is
famously *non-selective*. Aldol addition and its reverse (retro-aldol)
continuously scramble formaldehyde units in and out of sugars of every
size, so a real formose reaction produces a complex mixture ("formose
tar") of many sugars, not clean ribose. This module can't distinguish
ribose from its stereoisomers either -- sugars are tracked only by carbon
count, so "C5 sugar" here means the *entire pentose pool* (ribose,
arabinose, xylose, lyxose, and the ketopentoses), not ribose specifically.
See RIBOSE_FRACTION_ESTIMATE below for a documented, non-simulated way to
put a ballpark number on that.

**The proposed resolution in the literature** (Ricardo, Carrigan, Olcott,
Benner, "Borate Minerals Stabilize Ribose", Science 2004): borate ions
selectively bind ribose's ring geometry and protect it from further
reaction, pulling it out of the aldol/retro-aldol equilibrium before it
gets scrambled into hexoses or degraded. That's modeled here as an
optional, one-way "stabilization" reaction at a chosen carbon number
(default 5): once a sugar of that size is stabilized, it's terminal --
same "protected from further chemistry" treatment already used elsewhere
in this project for e.g. H2O or NH3.

The central experiment this module exists for: does a persistent C5 sugar
pool need that rescue mechanism, or does formose chemistry alone let
ribose-sized sugars accumulate on their own?
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .molecule import Species
from .reactions import Reaction
from .simulator import SimResult

RIBOSE_FRACTION_ESTIMATE = 1 / 4
"""A documented literature-style heuristic, NOT a simulated result: of the
four D/L aldopentose diastereomer pairs formose chemistry produces roughly
statistically without a selective catalyst (ribose, arabinose, xylose,
lyxose -- ignoring ketopentoses for this rough estimate), ribose is one.
Multiply a simulated "C5 sugar" count by this to get a ballpark
ribose-equivalent number. Used when `track_pentose_stereoisomers` is off
(the default). When it's on, the simulation reports actual ribose counts
instead of needing this -- see FormoseParams.ribose_selectivity below,
which is the real, tunable variable this estimate was standing in for."""

PENTOSE_DIASTEREOMERS = ("ribose", "arabinose", "xylose", "lyxose")
"""The four aldopentose diastereomer classes (D/L collapsed -- borate, an
achiral ion, doesn't distinguish enantiomers, only the relative
"cis/trans" diol geometry that differs between these four). Real formose
aldol addition, with no chiral catalyst present, is understood to produce
these roughly statistically -- the four get equal formation weight in
`build_formose_reactions` when stereoisomer tracking is on. What's NOT
equal in reality is how well borate binds each one afterward (ribose's
furanose form has a favorable cis-diol for borate binding relative to its
stereoisomers) -- that's what `ribose_selectivity` represents."""


class _Formaldehyde:
    """HCHO, the formose reaction's food molecule."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "HCHO"

    def formula(self) -> str:
        return "HCHO"

    @property
    def n_carbon(self) -> int:
        return 1

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, _Formaldehyde)

    def __hash__(self):
        return hash("HCHO")

    def __repr__(self):
        return "Molecule(HCHO)"


class _CannizzaroWaste:
    """Lumped methanol + formate, the product of the Cannizzaro
    disproportionation (2 HCHO -> CH3OH + HCOOH) -- a real, well-known
    side reaction that consumes formaldehyde without contributing to sugar
    growth, one of the things that limits real formose yield. n_carbon=2
    (both product carbons are tracked) so carbon accounting balances
    exactly against however much formaldehyde was seeded."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "waste"

    def formula(self) -> str:
        return "CH3OH+HCOOH"

    @property
    def n_carbon(self) -> int:
        return 2

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, _CannizzaroWaste)

    def __hash__(self):
        return hash("waste")

    def __repr__(self):
        return "Molecule(waste)"


FORMALDEHYDE = _Formaldehyde()
CANNIZZARO_WASTE = _CannizzaroWaste()


class Sugar:
    """A (CH2O)_n sugar, tracked by carbon count and, optionally, a named
    stereoisomer variant (only meaningful at the tracked target carbon
    number -- see PENTOSE_DIASTEREOMERS). Reactive: can still grow
    (aldol), shrink (retro-aldol), until stabilized."""

    __slots__ = ("n", "variant")

    def __init__(self, n: int, variant: str | None = None):
        self.n = n
        self.variant = variant

    def canonical_id(self) -> str:
        return f"C{self.n}sugar" + (f"-{self.variant}" if self.variant else "")

    def formula(self) -> str:
        base = f"(CH2O){self.n}"
        return f"{base} [{self.variant}]" if self.variant else base

    @property
    def n_carbon(self) -> int:
        return self.n

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, Sugar) and other.n == self.n and other.variant == self.variant

    def __hash__(self):
        return hash(("sugar", self.n, self.variant))

    def __repr__(self):
        return f"Molecule({self.formula()})"


class StabilizedSugar:
    """A sugar pulled out of the reactive pool by mineral (borate)
    binding -- terminal, no further reactions, same treatment already
    used elsewhere in this project for e.g. H2O or NH3."""

    __slots__ = ("n", "variant")

    def __init__(self, n: int, variant: str | None = None):
        self.n = n
        self.variant = variant

    def canonical_id(self) -> str:
        return f"C{self.n}sugar-stabilized" + (f"-{self.variant}" if self.variant else "")

    def formula(self) -> str:
        base = f"(CH2O){self.n}·borate"
        return f"{base} [{self.variant}]" if self.variant else base

    @property
    def n_carbon(self) -> int:
        return self.n

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, StabilizedSugar) and other.n == self.n and other.variant == self.variant

    def __hash__(self):
        return hash(("stabilized", self.n, self.variant))

    def __repr__(self):
        return f"Molecule({self.formula()})"


@dataclass
class FormoseParams:
    k_init: float = 0.01        # 2 HCHO -> C2 sugar (glycolaldehyde); the slow, rate-limiting "induction" step
    k_aldol: float = 1.0        # Cn + HCHO -> C(n+1); the autocatalytic chain-growth step
    k_retro_aldol: float = 0.3  # C(n+1) -> Cn + HCHO; real formose equilibria run both ways -- this is
                                 # the heart of the "sugar problem": nothing here is one-way
    k_cannizzaro: float = 0.05  # 2 HCHO -> waste; competing dead end that limits how much
                                 # formaldehyde is ever available for sugar growth at all
    k_stabilize: float = 0.0    # Cn -> stabilized Cn at n=stabilize_carbon; OFF by default,
                                 # i.e. the "no mineral rescue" case
    stabilize_carbon: int = 5   # which sugar size gets mineral protection (5 = ribose's carbon count)
    max_sugar_carbon: int = 8   # complexity ceiling, same role as max_carbon in the hydrocarbon engine
    track_pentose_stereoisomers: bool = False  # split the C5 tier into the 4 real named
                                 # diastereomers (only takes effect when stabilize_carbon == 5).
                                 # OFF by default -- preserves the simple generic-pool model.
    ribose_selectivity: float = 1.0  # multiplies k_stabilize for ribose specifically, relative to
                                 # its 3 stereoisomers (arabinose/xylose/lyxose), which all keep the
                                 # base k_stabilize. 1.0 = no selectivity (matches real borate binding
                                 # ribose somewhat preferentially -- Ricardo et al. 2004 -- so values
                                 # above 1 are the realistic direction to explore). Only matters when
                                 # track_pentose_stereoisomers is on; this is the actual tunable
                                 # variable behind RIBOSE_FRACTION_ESTIMATE.
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25
    seed: int | None = None


def _sugar_variants(n: int, track: bool) -> list[str | None]:
    """Which named variants exist at carbon number n. Only the tracked
    target tier (n == 5, when stereoisomer tracking is on) has more than
    one; every other tier is the single generic, untagged sugar."""
    if track and n == 5:
        return list(PENTOSE_DIASTEREOMERS)
    return [None]


def build_formose_reactions(p: FormoseParams) -> list[Reaction]:
    """The formose network is fully known upfront (unlike the hydrocarbon
    engine's combinatorial isomer explosion) -- no on-demand discovery
    needed, just build the whole fixed reaction list once.

    When track_pentose_stereoisomers is on, the carbon-5 tier is split
    into the 4 real named diastereomers (see PENTOSE_DIASTEREOMERS).
    Aldol addition treats them as equally likely (real formose aldol
    chemistry has no stereoselectivity without a chiral catalyst) --
    whichever reaction touches the C5 tier from a generic neighbor
    (C4 or C6) becomes 4 parallel, equally-weighted reactions, one per
    variant, so each gets ~1/4 of that flux automatically through
    ordinary Gillespie competition. Stabilization is the one place real
    selectivity enters: ribose's rate is multiplied by
    `ribose_selectivity` relative to its 3 stereoisomers."""
    track = p.track_pentose_stereoisomers and p.stabilize_carbon == 5
    out: list[Reaction] = []

    out.append(Reaction(
        kind="formose_init", reactant_ids=("HCHO", "HCHO"), products=(Sugar(2),),
        rate_constant=p.k_init, weight=1.0, key="formose_init:HCHO+HCHO->C2sugar",
    ))
    out.append(Reaction(
        kind="cannizzaro", reactant_ids=("HCHO", "HCHO"), products=(CANNIZZARO_WASTE,),
        rate_constant=p.k_cannizzaro, weight=1.0, key="cannizzaro:HCHO+HCHO->waste",
    ))

    for n in range(2, p.max_sugar_carbon):
        variants_n = _sugar_variants(n, track)
        variants_n1 = _sugar_variants(n + 1, track)

        # aldol: Cn + HCHO -> C(n+1)
        if len(variants_n1) > 1:
            for v in variants_n1:
                out.append(Reaction(
                    kind="aldol", reactant_ids=tuple(sorted((f"C{n}sugar", "HCHO"))),
                    products=(Sugar(n + 1, v),), rate_constant=p.k_aldol, weight=1.0,
                    key=f"aldol:C{n}sugar+HCHO->C{n + 1}sugar-{v}",
                ))
        elif len(variants_n) > 1:
            for v in variants_n:
                rid = f"C{n}sugar-{v}"
                out.append(Reaction(
                    kind="aldol", reactant_ids=tuple(sorted((rid, "HCHO"))),
                    products=(Sugar(n + 1),), rate_constant=p.k_aldol, weight=1.0,
                    key=f"aldol:{rid}+HCHO->C{n + 1}sugar",
                ))
        else:
            out.append(Reaction(
                kind="aldol", reactant_ids=tuple(sorted((f"C{n}sugar", "HCHO"))),
                products=(Sugar(n + 1),), rate_constant=p.k_aldol, weight=1.0,
                key=f"aldol:C{n}sugar+HCHO->C{n + 1}sugar",
            ))

        # retro-aldol: C(n+1) -> Cn + HCHO
        if len(variants_n1) > 1:
            for v in variants_n1:
                rid = f"C{n + 1}sugar-{v}"
                out.append(Reaction(
                    kind="retro_aldol", reactant_ids=(rid,),
                    products=(Sugar(n), FORMALDEHYDE), rate_constant=p.k_retro_aldol, weight=1.0,
                    key=f"retro_aldol:{rid}->C{n}sugar+HCHO",
                ))
        elif len(variants_n) > 1:
            for v in variants_n:
                out.append(Reaction(
                    kind="retro_aldol", reactant_ids=(f"C{n + 1}sugar",),
                    products=(Sugar(n, v), FORMALDEHYDE), rate_constant=p.k_retro_aldol, weight=1.0,
                    key=f"retro_aldol:C{n + 1}sugar->C{n}sugar-{v}+HCHO",
                ))
        else:
            out.append(Reaction(
                kind="retro_aldol", reactant_ids=(f"C{n + 1}sugar",),
                products=(Sugar(n), FORMALDEHYDE), rate_constant=p.k_retro_aldol, weight=1.0,
                key=f"retro_aldol:C{n + 1}sugar->C{n}sugar+HCHO",
            ))

    if p.k_stabilize > 0 and p.stabilize_carbon >= 2:
        variants = _sugar_variants(p.stabilize_carbon, track)
        if len(variants) > 1:
            for v in variants:
                rate = p.k_stabilize * (p.ribose_selectivity if v == "ribose" else 1.0)
                rid = f"C{p.stabilize_carbon}sugar-{v}"
                out.append(Reaction(
                    kind="stabilization", reactant_ids=(rid,),
                    products=(StabilizedSugar(p.stabilize_carbon, v),),
                    rate_constant=rate, weight=1.0,
                    key=f"stabilization:{rid}",
                ))
        else:
            out.append(Reaction(
                kind="stabilization", reactant_ids=(f"C{p.stabilize_carbon}sugar",),
                products=(StabilizedSugar(p.stabilize_carbon),),
                rate_constant=p.k_stabilize, weight=1.0,
                key=f"stabilization:C{p.stabilize_carbon}sugar",
            ))
    return out


def run_formose(p: FormoseParams, initial_hcho: int) -> SimResult:
    """A standalone Gillespie SSA loop over the static formose reaction
    list -- deliberately not the dynamic-discovery Simulator from
    simulator.py, since that machinery exists to handle a species set that
    isn't known in advance (the hydrocarbon isomer explosion); here it
    always is. Returns the same SimResult shape so it plugs into the same
    plotting code as the rest of the app."""
    reactions = build_formose_reactions(p)
    track = p.track_pentose_stereoisomers and p.stabilize_carbon == 5

    species: dict[str, Species] = {"HCHO": FORMALDEHYDE, "waste": CANNIZZARO_WASTE}
    for n in range(2, p.max_sugar_carbon + 1):
        variants = _sugar_variants(n, track)
        if len(variants) > 1:
            for v in variants:
                species[f"C{n}sugar-{v}"] = Sugar(n, v)
        else:
            species[f"C{n}sugar"] = Sugar(n)
    if p.k_stabilize > 0 and p.stabilize_carbon >= 2:
        variants = _sugar_variants(p.stabilize_carbon, track)
        if len(variants) > 1:
            for v in variants:
                species[f"C{p.stabilize_carbon}sugar-stabilized-{v}"] = StabilizedSugar(p.stabilize_carbon, v)
        else:
            species[f"C{p.stabilize_carbon}sugar-stabilized"] = StabilizedSugar(p.stabilize_carbon)

    counts = {sid: 0 for sid in species}
    counts["HCHO"] = initial_hcho

    rng = random.Random(p.seed)
    t = 0.0
    events = 0
    history_t = [0.0]
    history_counts = [dict(counts)]
    event_log: list[tuple[float, str]] = []
    reaction_fire_counts: dict[str, int] = {}
    stopped_reason = "t_max reached"

    while t < p.t_max and events < p.max_events:
        propensities = [r.propensity(counts) for r in reactions]
        a0 = sum(propensities)
        if a0 <= 0:
            stopped_reason = "no reaction has nonzero propensity (system exhausted)"
            break

        dt = -math.log(rng.random()) / a0
        t += dt
        if t > p.t_max:
            t = p.t_max
            break

        target = rng.random() * a0
        cum = 0.0
        chosen = reactions[-1]
        for r, a in zip(reactions, propensities):
            cum += a
            if cum >= target:
                chosen = r
                break

        for rid in chosen.reactant_ids:
            counts[rid] -= 1
        for product in chosen.products:
            counts[product.canonical_id()] += 1

        events += 1
        reaction_fire_counts[chosen.key] = reaction_fire_counts.get(chosen.key, 0) + 1
        event_log.append((t, chosen.key))

        if events % p.sample_every == 0:
            history_t.append(t)
            history_counts.append(dict(counts))
    else:
        if events >= p.max_events:
            stopped_reason = "max_events reached"

    history_t.append(t)
    history_counts.append(dict(counts))

    return SimResult(
        species=species,
        counts=counts,
        reactions=reactions,
        history_t=history_t,
        history_counts=history_counts,
        event_log=event_log,
        reaction_fire_counts=reaction_fire_counts,
        stopped_reason=stopped_reason,
    )
