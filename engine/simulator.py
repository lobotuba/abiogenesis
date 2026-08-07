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
    k_photo_o2: float = 0.03   # O2 -> 2 O* photolysis rate constant
    k_photo_co2: float = 0.005  # CO2 -> CO + O* photolysis rate constant (needs harder UV, so lower by default)
    k_o2_scavenge: float = 5.0  # R* + O2 -> ROO* rate constant; same order as k_comb by
                                 # default so the O2-vs-self-combination outcome is driven
                                 # by relative *concentration*, not a thumb on the scale
    k_o3_formation: float = 5.0  # O* + O2 -> O3 rate constant (near diffusion-limited, like k_comb)
    k_photo_o3: float = 0.08     # O3 -> O2 + O* photolysis rate constant (weaker bond than O2, so higher than k_photo_o2)
    k_o3_scavenge: float = 8.0   # R* + O3 -> RO* + O2; ozone reacts with radicals faster than O2 does in reality
    k_photo_h2o: float = 0.004   # H2O -> H* + OH* photolysis rate constant (weak absorber, needs harder UV than O2/CO2)
    k_oh_disprop: float = 2.0    # OH* + OH* -> H2O + O*; the secondary step that lets water photolysis
                                  # bootstrap free O2/O3 without any O2 present to begin with
    k_escape_h2: float = 0.0     # H2 -> escapes to space (irreversible sink); off by default, opt-in
    k_escape_h: float = 0.0      # H* -> escapes to space (irreversible sink); off by default, opt-in
    k_discharge_n2: float = 0.0  # N2 -> 2 N* electric discharge (NOT UV); off by default, opt-in.
                                  # This is "electricity" as a variable distinct from k_photo: UV in
                                  # this model can never touch N2 (its bond is too strong), but a
                                  # spark can -- see engine/reactions.py.
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
                k_photo_o2=self.p.k_photo_o2,
                k_photo_co2=self.p.k_photo_co2,
                k_o2_scavenge=self.p.k_o2_scavenge,
                k_o3_formation=self.p.k_o3_formation,
                k_photo_o3=self.p.k_photo_o3,
                k_o3_scavenge=self.p.k_o3_scavenge,
                k_photo_h2o=self.p.k_photo_h2o,
                k_oh_disprop=self.p.k_oh_disprop,
                k_escape_h2=self.p.k_escape_h2,
                k_escape_h=self.p.k_escape_h,
                k_discharge_n2=self.p.k_discharge_n2,
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
