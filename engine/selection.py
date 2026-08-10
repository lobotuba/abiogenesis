"""
The question this whole project was built to ask, made literal at last:
given two heritable replicator sequences competing for the same finite
food supply, does the faster-copying one actually win -- and can that
faster variant even *arise* from copying error in the first place, rather
than having to be seeded in from outside?

polymer.py has one self-complementary replicator species. Real Darwinian
selection needs three ingredients: (1) heritable **variation** -- more than
one distinguishable replicating type, (2) a copying process that can
produce that variation from imperfect replication of an existing type
(**mutation**), and (3) **differential survival/reproduction** under a
shared constraint. This module adds all three on top of polymer.py's
already-tested chemistry, at the same "generic species, no explicit base
sequence" abstraction level formose.py and nucleotide.py use for their
diastereomers -- named variant tags ("A", "B", ...) stand in for distinct
heritable sequences, the same way "ribose"/"arabinose"/etc. stand in for
distinct 3D stereoisomers elsewhere in this project.

**The reaction network**, per named variant `v`, reusing `Oligomer`/
`Duplex` from `engine/polymer.py` with `variant=v`:

1. Oligomerization: 3 RIBONUCLEOTIDE -> FRAGMENT_v. Rate is the SAME
   `k_oligomerize` for every variant -- deliberately unbiased, so every
   lineage has equal, fair access to the shared food supply. This is what
   makes it a genuine competition rather than a rigged one.
2. Background ligation: FRAGMENT_v + FRAGMENT_v -> TEMPLATE_v. Also
   variant-unbiased (`k_background`) -- ordinary condensation chemistry
   doesn't care about self-complementary structure, only templating does.
3. **Templated ligation, where fitness actually lives**: FRAGMENT_v +
   FRAGMENT_v + TEMPLATE_v -> 2 TEMPLATE_v, at `k_template[v]` -- each
   variant's OWN rate. This is the one place a "better" sequence (a
   structure that hybridizes and ligates its own fragments more
   efficiently) is allowed to matter, matching the real biophysical claim:
   fitness differences between self-replicators come from how well each
   one catalyzes its own copying, not from how easily its raw material
   forms.
4. Duplex formation / melting: TEMPLATE_v + TEMPLATE_v -> DUPLEX_v -> back,
   same self-inhibition mechanism as polymer.py, `k_duplex`/`k_melt` shared
   across variants. Deliberately no cross-variant duplexes (TEMPLATE_A +
   TEMPLATE_B -> hybrid) -- modeling that as absent is the simplifying
   assumption that distinct heritable sequences are different enough not
   to cross-hybridize, the same role sequence *specificity* plays in real
   template-directed replication (a real replicator generally can't use a
   mismatched template).
5. **Mutation** (off by default, `k_mutation`): during templated ligation,
   an existing TEMPLATE_v occasionally produces a TEMPLATE_w of a
   *different* tracked variant instead of a faithful copy of itself
   (FRAGMENT_v + FRAGMENT_v + TEMPLATE_v -> TEMPLATE_v + TEMPLATE_w, one
   such reaction per ordered pair v != w). This is the actual source of
   heritable variation this project's Limitations section flagged as
   missing: with mutation on, a variant that was never seeded can still
   *appear* through imperfect copying of an existing one, and then rise or
   fall on its own merits exactly like any other replicator. (For the
   cleanest version of that experiment, also set `k_background = 0` --
   with it nonzero, undirected background chemistry could independently
   stumble into forming any tracked variant's TEMPLATE from scratch too,
   which is a real if slow phenomenon in its own right, but it muddies a
   clean "mutation is what introduced this variant" comparison.)

The central experiment this module exists for: seed only one variant, give
a second (never-seeded) variant a faster `k_template`, and turn mutation
on (with `k_background = 0`, per the note above). Does the faster variant
emerge from copying error and then take over the population -- the literal
minimal signature of Darwinian selection
this whole project set out to test for?
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from itertools import permutations

from .molecule import Species
from .nucleotide import RIBONUCLEOTIDE
from .polymer import Duplex, Oligomer
from .reactions import Reaction
from .simulator import SimResult


@dataclass
class SelectionParams:
    variants: tuple[str, ...] = ("A", "B")
    k_oligomerize: float = 0.01   # shared across variants -- unbiased access to food
    k_background: float = 0.01    # shared across variants -- ordinary, non-selective ligation
    k_template: dict[str, float] = field(default_factory=lambda: {"A": 1.0, "B": 1.0})
    # per-variant templated-ligation rate -- this is where fitness differences live.
    # Equal by default (1.0/1.0 = no selection pressure), matching the project's
    # convention of shipping selective mechanisms at their neutral/off setting.
    k_duplex: float = 0.0         # shared across variants; off by default, same as polymer.py
    k_melt: float = 0.0           # shared across variants
    k_mutation: float = 0.0       # off by default -- with it off, only seeded variants can
                                    # ever exist; turning it on lets copying itself generate
                                    # variants that were never seeded at all
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25
    seed: int | None = None


def build_selection_reactions(p: SelectionParams) -> list[Reaction]:
    rid = RIBONUCLEOTIDE.canonical_id()
    out: list[Reaction] = []

    for v in p.variants:
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
            kind="templated_ligation", reactant_ids=(fid, fid, tid),
            products=(template_v, template_v), rate_constant=p.k_template.get(v, 0.0), weight=1.0,
            key=f"templated_ligation:2fragment-{v}+template-{v}->2template-{v}",
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

    if p.k_mutation > 0:
        for v, w in permutations(p.variants, 2):
            fragment_v = Oligomer(3, v)
            template_v = Oligomer(6, v)
            template_w = Oligomer(6, w)
            fid, tid = fragment_v.canonical_id(), template_v.canonical_id()
            out.append(Reaction(
                kind="mutation", reactant_ids=(fid, fid, tid),
                products=(template_v, template_w), rate_constant=p.k_mutation, weight=1.0,
                key=f"mutation:2fragment-{v}+template-{v}->template-{v}+template-{w}",
            ))

    return out


def run_selection(p: SelectionParams, initial_counts: dict[str, int]) -> SimResult:
    """initial_counts keys are canonical ids -- typically
    RIBONUCLEOTIDE.canonical_id() (the shared food supply) and, for
    whichever variants should start already existing,
    Oligomer(6, variant).canonical_id() (a seed TEMPLATE count). A variant
    with no seed and k_mutation == 0 can never appear at all -- that's the
    control case this module's central experiment compares against. Same
    standalone static-network Gillespie loop pattern as the rest of this
    project's chemistry-regime modules."""
    reactions = build_selection_reactions(p)

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
