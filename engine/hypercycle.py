"""
Item 3's last remaining piece: cooperation between distinct replicators.
selection.py showed two replicators *competing* for a shared food supply,
each relying entirely on its own copying ability. Real origin-of-life
theory also has a famous answer for the opposite arrangement: Eigen &
Schuster's **hypercycle** (*The Hypercycle: A Principle of Natural
Self-Organization*, 1979; Eigen, *Naturwissenschaften*, 1971) -- a closed
loop of replicators where each one catalyzes replication of the *next* one
in the cycle, not itself. The classic point of a hypercycle is that it lets
replicators too weak to sustain themselves alone persist and grow *only
because the loop is closed* -- break the loop into an open chain and the
same set of species, with the same rates, cannot sustain themselves the
same way.

**The reaction network**, per named variant `v` at position `i` in
`variants`, reusing `Oligomer`/`Duplex` from `engine/polymer.py` (same
technique as `engine/selection.py`):

1. Oligomerization and background ligation, shared/unbiased across
   variants, same as selection.py.
2. **Self-templated ligation** (`k_self`, per variant, 0.0 by default for
   every variant not given an explicit rate): FRAGMENT_v + FRAGMENT_v +
   TEMPLATE_v -> 2 TEMPLATE_v. Deliberately OFF by default for every
   member -- the classic hypercycle setup studies replicators that are
   individually *too weak to be self-sufficient*, so whatever growth
   happens has to come from cooperation, not solo copying.
3. **Cross-catalyzed ligation, the actual hypercycle step**: FRAGMENT_w +
   FRAGMENT_w + TEMPLATE_v -> TEMPLATE_v + TEMPLATE_w, where `w` is the
   *next* variant after `v` in the cycle. TEMPLATE_v is a genuine catalyst
   here (not consumed) for ligating two of ITS NEIGHBOR's fragments into a
   new copy of that neighbor -- the defining move of a hypercycle, and
   structurally the same three-reactant reaction shape as
   `templated_ligation` (self) and `mutation` (error) elsewhere in this
   project's replicator modules, just pointed at a different target.
4. Duplex formation/melting, same self-inhibition mechanism as
   polymer.py/selection.py, shared across variants.

`closed` (bool) controls whether the LAST variant's cross-catalyzed
reaction wraps back around to the FIRST variant (a true closed loop,
A -> B -> C -> A) or is simply omitted (an open chain, A -> B -> C, where
nothing ever catalyzes A). This is the whole experiment: build the exact
same set of species and rates, and ask whether closing that one edge
changes what the system can do.

**The central experiment this module exists for**: with every member's
self-templating rate at 0 (nobody can replicate alone) and background
ligation also at 0 (so a TEMPLATE's count can only ever change by being
produced as a reaction product), does closing the loop let a member that
would otherwise be permanently frozen at its starting count -- because
literally nothing in an open chain ever produces more of it -- actually
grow? This is about as sharp a test as a stochastic model can give: in the
open-chain case the very first variant's count is not just "slower to
grow," it is mathematically incapable of changing at all, since no
reaction in the whole network has it as a product.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from .molecule import Species
from .nucleotide import RIBONUCLEOTIDE
from .polymer import Duplex, Oligomer
from .reactions import Reaction
from .simulator import SimResult


@dataclass
class HypercycleParams:
    variants: tuple[str, ...] = ("A", "B", "C")
    closed: bool = True    # True: A->B->C->A (a real cycle). False: A->B->C (open chain, nothing helps A)
    k_oligomerize: float = 0.0005   # shared across variants -- unbiased access to food
    k_background: float = 0.0       # shared across variants; keep 0 for the cleanest "only the
                                      # cycle can grow anything" comparison (see module docstring)
    k_self: dict[str, float] = field(default_factory=dict)   # per-variant self-templating rate;
                                      # a variant not present as a key gets 0.0 -- individually
                                      # non-self-sufficient by default, the classic hypercycle setup
    k_cross: float = 1.0            # shared cross-catalysis rate: v catalyzes ligation of the
                                      # NEXT variant's own fragments into a new copy of that neighbor
    k_duplex: float = 0.0
    k_melt: float = 0.0
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25
    seed: int | None = None


def build_hypercycle_reactions(p: HypercycleParams) -> list[Reaction]:
    rid = RIBONUCLEOTIDE.canonical_id()
    n = len(p.variants)
    out: list[Reaction] = []

    for i, v in enumerate(p.variants):
        fragment_v = Oligomer(3, v)
        template_v = Oligomer(6, v)
        duplex_v = Duplex(v)
        fid, tid = fragment_v.canonical_id(), template_v.canonical_id()

        out.append(Reaction(
            kind="oligomerization", reactant_ids=(rid, rid, rid), products=(fragment_v,),
            rate_constant=p.k_oligomerize, weight=1.0,
            key=f"oligomerization:3ribonucleotide->fragment-{v}",
        ))
        out.append(Reaction(
            kind="background_ligation", reactant_ids=(fid, fid), products=(template_v,),
            rate_constant=p.k_background, weight=1.0,
            key=f"background_ligation:2fragment-{v}->template-{v}",
        ))
        out.append(Reaction(
            kind="self_templated_ligation", reactant_ids=(fid, fid, tid),
            products=(template_v, template_v), rate_constant=p.k_self.get(v, 0.0), weight=1.0,
            key=f"self_templated_ligation:2fragment-{v}+template-{v}->2template-{v}",
        ))
        out.append(Reaction(
            kind="duplex_formation", reactant_ids=(tid, tid), products=(duplex_v,),
            rate_constant=p.k_duplex, weight=1.0,
            key=f"duplex_formation:2template-{v}->duplex-{v}",
        ))
        out.append(Reaction(
            kind="duplex_melting", reactant_ids=(duplex_v.canonical_id(),),
            products=(template_v, template_v), rate_constant=p.k_melt, weight=1.0,
            key=f"duplex_melting:duplex-{v}->2template-{v}",
        ))

        is_wraparound_edge = (i == n - 1)
        if n >= 2 and (not is_wraparound_edge or p.closed):
            w = p.variants[(i + 1) % n]
            fragment_w = Oligomer(3, w)
            template_w = Oligomer(6, w)
            out.append(Reaction(
                kind="cross_catalyzed_ligation",
                reactant_ids=(fragment_w.canonical_id(), fragment_w.canonical_id(), tid),
                products=(template_v, template_w), rate_constant=p.k_cross, weight=1.0,
                key=f"cross_catalyzed_ligation:2fragment-{w}+template-{v}->template-{v}+template-{w}",
            ))

    return out


def run_hypercycle(p: HypercycleParams, initial_counts: dict[str, int]) -> SimResult:
    """initial_counts keys are canonical ids -- typically
    RIBONUCLEOTIDE.canonical_id() (the shared food supply) and
    Oligomer(6, variant).canonical_id() for each member's starting
    TEMPLATE count. Same standalone static-network Gillespie loop pattern
    as the rest of this project's chemistry-regime modules."""
    reactions = build_hypercycle_reactions(p)

    species: dict[str, Species] = {RIBONUCLEOTIDE.canonical_id(): RIBONUCLEOTIDE}
    for v in p.variants:
        for sp in (Oligomer(3, v), Oligomer(6, v), Duplex(v)):
            species[sp.canonical_id()] = sp

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
