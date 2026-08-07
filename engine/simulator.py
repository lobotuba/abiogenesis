"""
Gillespie stochastic simulation (SSA) over a reaction network that grows
on the fly: reactions are only generated for a species once, the first time
it appears with nonzero count, and are generated against a shared registry
of Species objects (see engine/reactions.py). This avoids recomputing the
full O(n^2) candidate reaction list every step while still letting the
network's chemical complexity expand as new molecules are discovered.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from .molecule import Species
from .reactions import Reaction, generate_reactions_for_species


@dataclass
class SimParams:
    k_photo: float = 0.02      # UV photolysis rate constant (per equivalent C-H bond)
    k_comb: float = 5.0        # radical-radical combination rate constant
    k_abstr: float = 0.5       # H-abstraction base rate constant
    max_carbon: int = 6        # complexity ceiling: cap on carbons per molecule
    t_max: float = 500.0
    max_events: int = 200_000
    sample_every: int = 25     # record a history snapshot every N fired events
    seed: int | None = None


@dataclass
class SimResult:
    species: dict[str, Species]
    counts: dict[str, int]
    reactions: list[Reaction]
    history_t: list[float]
    history_counts: list[dict[str, int]]
    event_log: list[tuple[float, str]]        # (time, reaction.key) per fired event
    reaction_fire_counts: dict[str, int]       # reaction.key -> times fired
    stopped_reason: str


class Simulator:
    def __init__(self, params: SimParams):
        self.p = params
        self.rng = random.Random(params.seed)
        self.species: dict[str, Species] = {}
        self.counts: dict[str, int] = {}
        self.reactions: dict[str, Reaction] = {}
        self.reaction_fire_counts: dict[str, int] = {}

    def seed_species(self, sp: Species, count: int) -> None:
        self._add_species(sp, count)

    def _add_species(self, sp: Species, count: int) -> None:
        sid = sp.canonical_id()
        if sid not in self.species:
            new_reactions = generate_reactions_for_species(
                sid, sp, self.species,
                max_carbon=self.p.max_carbon,
                k_photo=self.p.k_photo,
                k_comb=self.p.k_comb,
                k_abstr=self.p.k_abstr,
            )
            for r in new_reactions:
                self.reactions.setdefault(r.key, r)
            self.species[sid] = sp
            self.counts[sid] = 0
        self.counts[sid] += count

    def run(self) -> SimResult:
        t = 0.0
        events = 0
        history_t = [0.0]
        history_counts = [dict(self.counts)]
        event_log: list[tuple[float, str]] = []
        stopped_reason = "t_max reached"

        while t < self.p.t_max and events < self.p.max_events:
            rxns = list(self.reactions.values())
            propensities = [r.propensity(self.counts) for r in rxns]
            a0 = sum(propensities)
            if a0 <= 0:
                stopped_reason = "no reaction has nonzero propensity (system exhausted)"
                break

            dt = -math.log(self.rng.random()) / a0
            t += dt
            if t > self.p.t_max:
                t = self.p.t_max
                break

            target = self.rng.random() * a0
            cum = 0.0
            chosen = rxns[-1]
            for r, a in zip(rxns, propensities):
                cum += a
                if cum >= target:
                    chosen = r
                    break

            self._fire(chosen)
            events += 1
            self.reaction_fire_counts[chosen.key] = self.reaction_fire_counts.get(chosen.key, 0) + 1
            event_log.append((t, chosen.key))

            if events % self.p.sample_every == 0:
                history_t.append(t)
                history_counts.append(dict(self.counts))
        else:
            if events >= self.p.max_events:
                stopped_reason = "max_events reached"

        history_t.append(t)
        history_counts.append(dict(self.counts))

        return SimResult(
            species=dict(self.species),
            counts=dict(self.counts),
            reactions=list(self.reactions.values()),
            history_t=history_t,
            history_counts=history_counts,
            event_log=event_log,
            reaction_fire_counts=dict(self.reaction_fire_counts),
            stopped_reason=stopped_reason,
        )

    def _fire(self, r: Reaction) -> None:
        for rid in r.reactant_ids:
            self.counts[rid] -= 1
        for product in r.products:
            self._add_species(product, 1)
