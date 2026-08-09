"""
Roadmap item 3: "Beyond the activated ribonucleotide." The nucleotide
module (engine/nucleotide.py) stops at one finished, activated
ribonucleotide monomer. This module asks the next question directly: can
those monomers polymerize into an actual strand, and can that process be
templated/copied -- i.e. does self-amplification (this project's very
first goal, back at the hydrocarbon-radical stage) show up again once
there's finally something strand-like to copy?

**The real system this is modeled on**: von Kiedrowski, "A Self-Replicating
Hexadeoxynucleotide" (Nature, 1986) -- the first experimentally
demonstrated minimal template-directed self-replicator. Its trick: make the
product strand (TEMPLATE here, a hexamer) *self-complementary*, so a single
TEMPLATE molecule can hybridize two copies of a shorter FRAGMENT (a trimer,
i.e. exactly half of TEMPLATE) side by side and template their ligation
into a brand new TEMPLATE. One species acts as substrate, catalyst, and
product all at once -- genuine autocatalysis, not a metaphor for it.

**The pathway modeled here** (same generic-species, no-explicit-sequence
abstraction level as nucleotide.py -- see Limitations in the README):

1. Oligomerization: 3 activated ribonucleotides (`nucleotide.RIBONUCLEOTIDE`,
   imported directly -- this pathway starts exactly where that module ends)
   spontaneously join into a FRAGMENT (a trimer). Deliberately slow/weak by
   default: condensation reactions release water, and in bulk aqueous
   solution that makes uncatalyzed ligation thermodynamically uphill, the
   same real difficulty noted throughout this project.
2. Background ligation: two FRAGMENTs join into a TEMPLATE (a hexamer),
   equally slow and uncatalyzed -- the "no templating exists yet" baseline.
3. **Templated ligation, the point of this module**: two FRAGMENTs plus an
   *existing* TEMPLATE ligate into two TEMPLATEs. The old TEMPLATE isn't
   consumed -- it's a genuine catalyst, so this reaction is autocatalytic
   by construction. This needed a real termolecular reaction (three
   reactants at once), which is why `Reaction.propensity` in
   engine/reactions.py was generalized from its old hand-written
   one-or-two-reactant special case to a general N-body combinatorial
   formula (`comb(n, multiplicity)` per distinct reactant id, multiplied
   together) -- a strict generalization that reduces to the exact old
   arithmetic for every reaction already in the project.
4. **Duplex formation, the other point of this module**: the same
   self-complementarity that lets one TEMPLATE catalyze step 3 also lets
   two TEMPLATEs hybridize with EACH OTHER into an inert double-stranded
   DUPLEX -- a dead end until it melts back apart. This isn't a bolted-on
   penalty; it's mechanistically the same base-pairing chemistry as step 3,
   just acting on two product strands instead of a product and two
   fragments. Von Kiedrowski's actual, famous, counterintuitive finding was
   that real minimal self-replicators of this type show *parabolic*
   (roughly sqrt(t)) growth, not exponential/Malthusian growth, precisely
   because of this self-inhibition -- the catalyst increasingly sequesters
   itself out of solution as it accumulates.
5. Duplex melting: DUPLEX -> 2 TEMPLATE, a slow reverse (thermal
   denaturation / wet-dry or day-night cycling in the real system),
   reintroducing catalytically active single strands.

The central experiment this module exists for: turning on the
autocatalytic step (3) should measurably speed up TEMPLATE accumulation
relative to background ligation alone -- but does the self-inhibition step
(4), which comes from the *same underlying chemistry*, cap that
amplification well short of a runaway? I.e., is real template-directed
self-replication actually as unconditionally powerful as "autocatalysis"
sounds, or is there a structural reason (built into the same mechanism that
makes it work at all) it can't just run away exponentially?
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .molecule import Species
from .nucleotide import RIBONUCLEOTIDE
from .reactions import Reaction
from .simulator import SimResult


class Oligomer:
    """A short RNA-like strand, tracked only by length (in ribonucleotide
    units) -- same no-explicit-sequence abstraction nucleotide.py uses.
    This module's reaction network only ever produces two lengths: FRAGMENT
    (3, a "half-strand") and TEMPLATE (6, self-complementary -- its own two
    halves are each a copy of FRAGMENT), which is the whole mechanism."""

    __slots__ = ("n",)

    def __init__(self, n: int):
        self.n = n

    def canonical_id(self) -> str:
        return f"oligo{self.n}"

    def formula(self) -> str:
        return f"(Nmp){self.n}"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, Oligomer) and other.n == self.n

    def __hash__(self):
        return hash(("oligomer", self.n))

    def __repr__(self):
        return f"Molecule({self.formula()})"


class Duplex:
    """Two TEMPLATE strands hybridized into an inert double helix --
    unreactive (no further reactions defined against it) until it melts
    back apart. The product-inhibition step: the same self-complementarity
    that lets one TEMPLATE catalyze ligation of two FRAGMENTs also lets two
    TEMPLATEs pair with EACH OTHER, taking both out of circulation."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "duplex"

    def formula(self) -> str:
        return "(Nmp)6:(Nmp)6 duplex"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list:
        return []

    def __eq__(self, other):
        return isinstance(other, Duplex)

    def __hash__(self):
        return hash("duplex")

    def __repr__(self):
        return "Molecule(duplex)"


FRAGMENT = Oligomer(3)
TEMPLATE = Oligomer(6)
DUPLEX = Duplex()


@dataclass
class PolymerParams:
    k_oligomerize: float = 0.01  # 3 ribonucleotide -> FRAGMENT; slow, uncatalyzed --
                                  # condensation releases water, disfavored in bulk solution
    k_background: float = 0.01   # FRAGMENT + FRAGMENT -> TEMPLATE, uncatalyzed; the
                                  # "no templating exists yet" baseline this module tests against
    k_template: float = 0.0      # FRAGMENT + FRAGMENT + TEMPLATE -> TEMPLATE + TEMPLATE;
                                  # OFF by default -- the autocatalytic step itself
    k_duplex: float = 0.0        # TEMPLATE + TEMPLATE -> DUPLEX; OFF by default. Not a bolted-on
                                  # penalty -- the same base-pairing chemistry that makes
                                  # k_template work at all, acting on two products instead of a
                                  # product and two fragments. Von Kiedrowski's real finding:
                                  # turning this on caps growth at parabolic, not exponential.
    k_melt: float = 0.0          # DUPLEX -> TEMPLATE + TEMPLATE; strand separation (thermal /
                                  # wet-dry cycling), reverses k_duplex
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25
    seed: int | None = None


def build_polymer_reactions(p: PolymerParams) -> list[Reaction]:
    rid = RIBONUCLEOTIDE.canonical_id()
    fid = FRAGMENT.canonical_id()
    tid = TEMPLATE.canonical_id()

    return [
        Reaction(
            kind="oligomerization", reactant_ids=(rid, rid, rid), products=(FRAGMENT,),
            rate_constant=p.k_oligomerize, weight=1.0,
            key="oligomerization:3ribonucleotide->fragment",
        ),
        Reaction(
            kind="background_ligation", reactant_ids=(fid, fid), products=(TEMPLATE,),
            rate_constant=p.k_background, weight=1.0,
            key="background_ligation:2fragment->template",
        ),
        Reaction(
            kind="templated_ligation", reactant_ids=(fid, fid, tid),
            products=(TEMPLATE, TEMPLATE), rate_constant=p.k_template, weight=1.0,
            key="templated_ligation:2fragment+template->2template",
        ),
        Reaction(
            kind="duplex_formation", reactant_ids=(tid, tid), products=(DUPLEX,),
            rate_constant=p.k_duplex, weight=1.0,
            key="duplex_formation:2template->duplex",
        ),
        Reaction(
            kind="duplex_melting", reactant_ids=(DUPLEX.canonical_id(),),
            products=(TEMPLATE, TEMPLATE), rate_constant=p.k_melt, weight=1.0,
            key="duplex_melting:duplex->2template",
        ),
    ]


def run_polymer(p: PolymerParams, initial_counts: dict[str, int]) -> SimResult:
    """initial_counts keys are canonical ids -- typically
    RIBONUCLEOTIDE.canonical_id() (the food supply) and, optionally,
    TEMPLATE.canonical_id() seeded at a small nonzero count to explore "a
    rare first copy already exists" versus "must arise from background
    ligation alone" (0, the default if omitted). Same standalone
    static-network Gillespie loop pattern as formose.py/nucleotide.py."""
    reactions = build_polymer_reactions(p)

    species: dict[str, Species] = {
        RIBONUCLEOTIDE.canonical_id(): RIBONUCLEOTIDE,
        FRAGMENT.canonical_id(): FRAGMENT,
        TEMPLATE.canonical_id(): TEMPLATE,
        DUPLEX.canonical_id(): DUPLEX,
    }

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

        for rid_ in chosen.reactant_ids:
            counts[rid_] -= 1
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
