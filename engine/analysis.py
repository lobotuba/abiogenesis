"""
Chain-length distribution helpers, shared by the single-run UI view and the
concentration/UV parameter-sweep tool.
"""
from __future__ import annotations

from dataclasses import dataclass

from .molecule import Molecule
from .simulator import SimResult


@dataclass
class ChainLengthStats:
    counts_by_carbon: dict[int, int]  # n_carbon -> total molecule count (hydrocarbons only)
    mean_carbon: float                # molecule-count-weighted mean carbon number
    max_carbon_present: int           # largest n_carbon with count > 0
    c3plus_carbon_fraction: float     # fraction of all hydrocarbon carbon atoms sitting in C3+ species


def chain_length_stats(result: SimResult) -> ChainLengthStats:
    """Carbon-number distribution over hydrocarbon species only (Molecule
    instances) -- excludes H*/H2/O2/N2/etc, which have no chain length."""
    counts_by_carbon: dict[int, int] = {}
    total_carbon = 0
    total_c3plus_carbon = 0
    total_molecules = 0
    weighted_carbon = 0

    for sid, sp in result.species.items():
        if not isinstance(sp, Molecule):
            continue
        n = result.counts.get(sid, 0)
        if n <= 0:
            continue
        counts_by_carbon[sp.n_carbon] = counts_by_carbon.get(sp.n_carbon, 0) + n
        total_molecules += n
        weighted_carbon += sp.n_carbon * n
        total_carbon += sp.n_carbon * n
        if sp.n_carbon >= 3:
            total_c3plus_carbon += sp.n_carbon * n

    mean_carbon = weighted_carbon / total_molecules if total_molecules else 0.0
    max_carbon_present = max(counts_by_carbon) if counts_by_carbon else 0
    c3plus_fraction = total_c3plus_carbon / total_carbon if total_carbon else 0.0

    return ChainLengthStats(
        counts_by_carbon=counts_by_carbon,
        mean_carbon=mean_carbon,
        max_carbon_present=max_carbon_present,
        c3plus_carbon_fraction=c3plus_fraction,
    )
