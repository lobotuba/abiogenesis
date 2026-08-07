"""
Heuristic detection of candidate autocatalytic / chain-propagation loops in
a *realized* reaction network (i.e. built only from reactions that actually
fired during a simulation run, not the full set of candidate reactions).

This is deliberately approximate. True autocatalysis requires a species to
be net-produced faster because of its own presence; what we detect here is
weaker but chemically meaningful: directed cycles in the reactant->product
flow graph, weighted by how often each contributing reaction actually fired.
A short, high-flux cycle is a candidate chain-propagation loop -- e.g. a
radical that gets regenerated a few steps after it's consumed, the closest
thing this hydrocarbon-only model has to "a catalyst that reproduces itself."

Longer-term (see README roadmap), genuine template-directed self-replication
would need a richer chemistry (e.g. heteroatom/backbone polymers) -- this
module is a diagnostic for whether the current run's network has the right
shape (cycles) to be a substrate for that, not a claim that it *is* alive.
"""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from .reactions import Reaction


@dataclass
class Cycle:
    node_ids: tuple[str, ...]      # species ids, in cycle order
    formulas: tuple[str, ...]      # matching display formulas
    bottleneck_flux: int           # min realized fire-count among edges in the cycle
    length: int


def build_flow_graph(
    reactions: list[Reaction],
    reaction_fire_counts: dict[str, int],
) -> nx.DiGraph:
    """Directed species graph: edge u->v weighted by summed fire-counts of
    every realized reaction where u was a reactant and v a product."""
    g = nx.DiGraph()
    for r in reactions:
        fired = reaction_fire_counts.get(r.key, 0)
        if fired <= 0:
            continue
        for rid in set(r.reactant_ids):
            for p in r.products:
                pid = p.canonical_id()
                if rid == pid:
                    continue
                if g.has_edge(rid, pid):
                    g[rid][pid]["weight"] += fired
                else:
                    g.add_edge(rid, pid, weight=fired)
    return g


def find_candidate_cycles(
    reactions: list[Reaction],
    reaction_fire_counts: dict[str, int],
    species: dict,
    *,
    max_cycle_length: int = 6,
    max_results: int = 25,
) -> list[Cycle]:
    g = build_flow_graph(reactions, reaction_fire_counts)
    if g.number_of_edges() == 0:
        return []

    cycles: list[Cycle] = []
    try:
        raw = nx.simple_cycles(g, length_bound=max_cycle_length)
    except TypeError:
        # older networkx without length_bound support
        raw = (c for c in nx.simple_cycles(g) if len(c) <= max_cycle_length)

    for i, cycle_nodes in enumerate(raw):
        if i >= 5000:  # hard safety cap on how many candidates we even score
            break
        if len(cycle_nodes) < 2:
            continue
        weights = []
        ok = True
        for a, b in zip(cycle_nodes, cycle_nodes[1:] + cycle_nodes[:1]):
            if not g.has_edge(a, b):
                ok = False
                break
            weights.append(g[a][b]["weight"])
        if not ok:
            continue
        formulas = tuple(species[n].formula() if n in species else n for n in cycle_nodes)
        cycles.append(Cycle(
            node_ids=tuple(cycle_nodes),
            formulas=formulas,
            bottleneck_flux=min(weights),
            length=len(cycle_nodes),
        ))

    cycles.sort(key=lambda c: (-c.bottleneck_flux, c.length))
    return cycles[:max_results]
