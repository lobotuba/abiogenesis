"""
A more direct route to a nucleotide than "make ribose, then attach a base
to it": Powner, Gerland & Sutherland, "Synthesis of activated pyrimidine
ribonucleotides in prebiotically plausible conditions" (Nature, 2009).

The formose module (engine/formose.py) asks whether ribose can accumulate
on its own and finds that it needs a rescue mechanism (mineral
stabilization) because free sugars are stuck in a scrambling
aldol/retro-aldol equilibrium (the "sugar problem"). Sutherland's insight
was different: **don't make free ribose at all.** Attach the nucleobase
precursor to the sugar chain *before* the ribose-defining stereocenters
even form, so there's no free-sugar mixture to purify ribose out of in the
first place.

**The pathway modeled here** (real, simplified to the same
generic-species-with-a-name-tag level of detail as the formose module --
see Limitations in the README):

1. Glycolaldehyde + cyanamide -> 2-aminooxazole. Glycolaldehyde is the
   very first product of formose chemistry (`Sugar(2)` in
   engine/formose.py, imported directly here) -- this pathway branches off
   *before* the sugar problem has any chance to develop.
2. 2-aminooxazole + glyceraldehyde (`Sugar(3)`) -> a mixture of
   ribo-/arabino-/xylo-/lyxo-configured aminooxazoline, formed with equal
   likelihood (this step is still an unselective aldol-type addition, just
   like formose's -- the nucleobase precursor being attached doesn't
   change that part of the chemistry).
3. **The selection step, and the whole point of this pathway**:
   Sutherland's group found the ribo-configured aminooxazoline is
   photostable, while the other three diastereomers are destroyed by UV at
   the wavelength used. Modeled as a photolysis reaction that only touches
   the 3 non-ribo variants -- ribo is immune, same mechanism-shape as
   formose's mineral stabilization (protect/remove the wrong material) but
   inverted (here it's cheaper to destroy the losers than protect the
   winner). Off by default, for comparison.
4. Ribo-aminooxazoline + cyanoacetylene -> anhydronucleoside (productive).
   Cyanoacetylene also reacts with the 3 non-ribo diastereomers at the same
   rate, but non-productively (an off-pathway dead-end adduct) -- this is
   the real reason selection (3) matters. It's tempting to assume
   destroying the "wrong" diastereomers would matter by *recycling*
   material back to glyceraldehyde/aminooxazole, but it doesn't (photolysis
   doesn't run the aldol step backwards) -- an earlier version of this
   module modeled it that way and found selection made *no* difference to
   final yield, which is what exposed the real mechanism: cyanoacetylene is
   the scarce, prebiotically special reagent here, and every unit of it
   consumed by a non-ribo diastereomer is a unit that can't go toward an
   actual nucleotide. Selection matters by *protecting a scarce downstream
   resource from being wasted on unproductive branches*, not by rescuing
   upstream material.
5. Anhydronucleoside + phosphate -> an activated pyrimidine ribonucleotide.

The central experiment this module exists for: does the photochemical
selection step (3) actually matter for how much nucleotide accumulates, or
does the pathway work fine without it? Only shows up clearly when
cyanoacetylene is seeded as the scarce reagent relative to the
aminooxazoline pool -- see the app's default starting amounts.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .formose import PENTOSE_DIASTEREOMERS, Sugar
from .molecule import Species
from .reactions import Reaction
from .simulator import SimResult


class _NamedFood:
    """A simple, unreactive-except-as-written-here starting material or
    terminal product. carbon bookkeeping is deliberately not attempted for
    this module (see README) -- these are opaque named intermediates, not
    exact molecular formulas."""

    __slots__ = ("_id", "_formula")

    def __init__(self, id_: str, formula: str):
        self._id = id_
        self._formula = formula

    def canonical_id(self) -> str:
        return self._id

    def formula(self) -> str:
        return self._formula

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, _NamedFood) and other._id == self._id

    def __hash__(self):
        return hash(("nucfood", self._id))

    def __repr__(self):
        return f"Molecule({self._formula})"


CYANAMIDE = _NamedFood("cyanamide", "NH2CN")
CYANOACETYLENE = _NamedFood("cyanoacetylene", "HC3N")
PHOSPHATE = _NamedFood("phosphate", "Pi")
AMINOOXAZOLE = _NamedFood("aminooxazole", "2-aminooxazole")
PHOTO_WASTE = _NamedFood("nuc-photo-waste", "photo-waste")
OFFPATHWAY_ADDUCT = _NamedFood("nuc-offpathway-adduct", "off-pathway adduct (dead end)")
ANHYDRONUCLEOSIDE = _NamedFood("anhydronucleoside", "anhydronucleoside")
RIBONUCLEOTIDE = _NamedFood("ribonucleotide", "ribonucleotide (activated)")

GLYCOLALDEHYDE = Sugar(2)
GLYCERALDEHYDE = Sugar(3)


class Aminooxazoline(_NamedFood):
    """One of the 4 named stereoisomer configurations formed when
    2-aminooxazole reacts with glyceraldehyde -- same diastereomer set as
    formose's pentoses (PENTOSE_DIASTEREOMERS), since these are the same
    stereocenters, just carried on a different scaffold."""

    __slots__ = ("variant",)

    def __init__(self, variant: str):
        self.variant = variant
        super().__init__(f"aminooxazoline-{variant}", f"{variant}-configured aminooxazoline")

    def __eq__(self, other):
        return isinstance(other, Aminooxazoline) and other.variant == self.variant

    def __hash__(self):
        return hash(("aminooxazoline", self.variant))


@dataclass
class NucleotideParams:
    k_oxazole: float = 1.0          # glycolaldehyde + cyanamide -> 2-aminooxazole
    k_aminooxazoline: float = 1.0   # aminooxazole + glyceraldehyde -> each of the 4 configured variants
    k_photo_destroy: float = 0.0    # UV destroys the 3 non-ribo variants; ribo is immune (the real,
                                     # documented Sutherland finding). OFF by default -- the "no
                                     # photochemical selection" case, for comparison.
    k_anhydro: float = 1.0          # ribo-aminooxazoline + cyanoacetylene -> anhydronucleoside
    k_phosphorylate: float = 1.0    # anhydronucleoside + phosphate -> ribonucleotide
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25
    seed: int | None = None


def build_nucleotide_reactions(p: NucleotideParams) -> list[Reaction]:
    out: list[Reaction] = []

    out.append(Reaction(
        kind="oxazole_formation",
        reactant_ids=tuple(sorted((GLYCOLALDEHYDE.canonical_id(), CYANAMIDE.canonical_id()))),
        products=(AMINOOXAZOLE,), rate_constant=p.k_oxazole, weight=1.0,
        key="oxazole_formation:glycolaldehyde+cyanamide->aminooxazole",
    ))

    for variant in PENTOSE_DIASTEREOMERS:
        out.append(Reaction(
            kind="aminooxazoline_formation",
            reactant_ids=tuple(sorted((AMINOOXAZOLE.canonical_id(), GLYCERALDEHYDE.canonical_id()))),
            products=(Aminooxazoline(variant),), rate_constant=p.k_aminooxazoline, weight=1.0,
            key=f"aminooxazoline_formation:aminooxazole+glyceraldehyde->{variant}",
        ))

    for variant in PENTOSE_DIASTEREOMERS:
        if variant == "ribose":
            continue  # photostable -- the whole point of this pathway
        out.append(Reaction(
            kind="photo_selection",
            reactant_ids=(Aminooxazoline(variant).canonical_id(),),
            products=(PHOTO_WASTE,), rate_constant=p.k_photo_destroy, weight=1.0,
            key=f"photo_selection:{variant}->waste",
        ))

    out.append(Reaction(
        kind="anhydro_formation",
        reactant_ids=tuple(sorted((Aminooxazoline("ribose").canonical_id(), CYANOACETYLENE.canonical_id()))),
        products=(ANHYDRONUCLEOSIDE,), rate_constant=p.k_anhydro, weight=1.0,
        key="anhydro_formation:ribose-aminooxazoline+cyanoacetylene->anhydronucleoside",
    ))

    # The non-ribo diastereomers, if not photochemically removed first, react with
    # cyanoacetylene just as readily as ribo does -- just unproductively. This is
    # what makes selection (step 3) actually matter: it protects the scarce
    # cyanoacetylene supply from being consumed by dead-end branches.
    for variant in PENTOSE_DIASTEREOMERS:
        if variant == "ribose":
            continue
        out.append(Reaction(
            kind="offpathway_consumption",
            reactant_ids=tuple(sorted((Aminooxazoline(variant).canonical_id(), CYANOACETYLENE.canonical_id()))),
            products=(OFFPATHWAY_ADDUCT,), rate_constant=p.k_anhydro, weight=1.0,
            key=f"offpathway_consumption:{variant}+cyanoacetylene->waste",
        ))

    out.append(Reaction(
        kind="phosphorylation",
        reactant_ids=tuple(sorted((ANHYDRONUCLEOSIDE.canonical_id(), PHOSPHATE.canonical_id()))),
        products=(RIBONUCLEOTIDE,), rate_constant=p.k_phosphorylate, weight=1.0,
        key="phosphorylation:anhydronucleoside+phosphate->ribonucleotide",
    ))

    return out


def run_nucleotide(p: NucleotideParams, initial_counts: dict[str, int]) -> SimResult:
    """initial_counts keys are canonical ids -- typically
    GLYCOLALDEHYDE/CYANAMIDE/GLYCERALDEHYDE/CYANOACETYLENE/PHOSPHATE's
    .canonical_id(). Same standalone static-network Gillespie loop pattern
    as engine/formose.py's run_formose (the species set is small and fully
    known upfront, no on-demand discovery needed)."""
    reactions = build_nucleotide_reactions(p)

    species: dict[str, Species] = {
        GLYCOLALDEHYDE.canonical_id(): GLYCOLALDEHYDE,
        GLYCERALDEHYDE.canonical_id(): GLYCERALDEHYDE,
        CYANAMIDE.canonical_id(): CYANAMIDE,
        CYANOACETYLENE.canonical_id(): CYANOACETYLENE,
        PHOSPHATE.canonical_id(): PHOSPHATE,
        AMINOOXAZOLE.canonical_id(): AMINOOXAZOLE,
        PHOTO_WASTE.canonical_id(): PHOTO_WASTE,
        OFFPATHWAY_ADDUCT.canonical_id(): OFFPATHWAY_ADDUCT,
        ANHYDRONUCLEOSIDE.canonical_id(): ANHYDRONUCLEOSIDE,
        RIBONUCLEOTIDE.canonical_id(): RIBONUCLEOTIDE,
    }
    for variant in PENTOSE_DIASTEREOMERS:
        aox = Aminooxazoline(variant)
        species[aox.canonical_id()] = aox

    counts = {sid: 0 for sid in species}
    for sid, n in initial_counts.items():
        counts[sid] = n

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
