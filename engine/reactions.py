"""
Rule-based reaction generator.

Three generic rules, applied to whatever species currently exist in the
pool, produce all candidate elementary reactions:

  1. Photolysis      M + hv -> M* (radical) + H*      (UV C-H homolysis)
  2. Combination      R* + R'* -> R-R'                 (radical-radical)
  3. Abstraction       R* + M-H -> R-H + M*             (H-atom transfer)

Rule (3) is what lets the network self-amplify: a radical can be regenerated
by a different chain a few steps later, which is the chemically real analog
of autocatalysis this project is built around.

Rates are deliberately not calibrated to real photon flux / Arrhenius
kinetics -- they're relative, adjustable knobs (see engine/simulator.py)
meant for exploring qualitative network behavior, not for quantitative
photochemistry.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .molecule import ATOMIC_H, Species, combine

ReactionKind = Literal["photolysis", "combination", "abstraction"]

# Crude C-H bond-strength proxy: more substituted carbons (higher skeleton
# degree) hold weaker C-H bonds and are abstracted/photolyzed preferentially.
# This is what gives the network a realistic-flavored bias toward branched,
# more-substituted radicals rather than a uniform random walk over isomers.
_DEGREE_WEIGHT = {0: 1.0, 1: 1.0, 2: 2.5, 3: 4.5, 4: 4.5}


def _degree_weight(degree: int) -> float:
    return _DEGREE_WEIGHT.get(degree, 4.5)


@dataclass(frozen=True)
class Reaction:
    kind: ReactionKind
    reactant_ids: tuple[str, ...]
    products: tuple[Species, ...]   # actual Species objects (ids alone can't rebuild a molecule graph)
    rate_constant: float            # base rate constant (a UI-adjustable knob, by kind)
    weight: float                   # statistical multiplicity * bond-strength factor
    key: str

    @property
    def product_ids(self) -> tuple[str, ...]:
        return tuple(p.canonical_id() for p in self.products)

    def propensity(self, counts: dict[str, int]) -> float:
        a, b = (self.reactant_ids + (None,))[:2]
        na = counts.get(a, 0)
        if len(self.reactant_ids) == 1:
            return self.rate_constant * self.weight * na
        nb = counts.get(b, 0)
        if a == b:
            pairs = na * (na - 1) / 2.0
        else:
            pairs = na * nb
        return self.rate_constant * self.weight * pairs


def _mk_key(kind: str, reactants: tuple[str, ...], products: tuple[str, ...]) -> str:
    return f"{kind}:{'+'.join(sorted(reactants))}->{'+'.join(sorted(products))}"


def generate_reactions_for_species(
    new_id: str,
    new_species: Species,
    pool: dict[str, Species],
    *,
    max_carbon: int,
    k_photo: float,
    k_comb: float,
    k_abstr: float,
) -> list[Reaction]:
    """Generate all new reactions enabled by `new_species` joining the pool.

    Called once per newly-discovered species against the full existing pool
    (including itself), so the candidate reaction list only grows when the
    species set grows -- not every SSA step.
    """
    out: list[Reaction] = []

    # 1. Photolysis: only closed-shell, carbon-containing species absorb UV
    #    and homolyze a C-H bond in this model.
    if new_species.n_carbon >= 1 and not new_species.is_radical:
        for site in new_species.radicalizable_sites():
            if site.product.n_carbon > max_carbon:
                continue
            products = (site.product, ATOMIC_H)
            product_ids = tuple(p.canonical_id() for p in products)
            out.append(Reaction(
                kind="photolysis",
                reactant_ids=(new_id,),
                products=products,
                rate_constant=k_photo,
                weight=float(site.multiplicity),
                key=_mk_key("photolysis", (new_id,), product_ids),
            ))

    # 2. Combination: new_species (if a radical) with every radical in the
    #    pool (including itself), and every existing radical with new_species
    #    if new_species is closed-shell -- covered by iterating all radicals.
    all_ids = list(pool.keys())
    if new_id not in all_ids:
        all_ids.append(new_id)
    radicals = [sid for sid in all_ids if _species_for(sid, new_id, new_species, pool).is_radical]
    if new_species.is_radical:
        for other_id in radicals:
            other = _species_for(other_id, new_id, new_species, pool)
            product = combine(new_species, other)
            if product.n_carbon > max_carbon:
                continue
            reactants = tuple(sorted((new_id, other_id)))
            out.append(Reaction(
                kind="combination",
                reactant_ids=reactants,
                products=(product,),
                rate_constant=k_comb,
                weight=1.0,
                key=_mk_key("combination", reactants, (product.canonical_id(),)),
            ))

    # 3. Abstraction: radicals pull an H off closed-shell H-bearing species.
    #    Cover both directions: new_species as the abstracting radical, and
    #    new_species as the H-donor.
    def add_abstraction(radical_id: str, radical: Species, donor_id: str, donor: Species):
        if radical is donor or radical_id == donor_id:
            return
        if donor.is_radical:
            return  # v1: no biradical formation, donors must be closed-shell
        for site in donor.radicalizable_sites():
            if site.product.n_carbon > max_carbon:
                continue
            saturated = combine(radical, ATOMIC_H)
            products = (saturated, site.product)
            product_ids = tuple(sorted((saturated.canonical_id(), site.product.canonical_id())))
            reactants = (radical_id, donor_id)
            out.append(Reaction(
                kind="abstraction",
                reactant_ids=reactants,
                products=products,
                rate_constant=k_abstr,
                weight=float(site.multiplicity) * _degree_weight(site.degree),
                key=_mk_key("abstraction", reactants, product_ids),
            ))

    if new_species.is_radical:
        for other_id, other in pool.items():
            add_abstraction(new_id, new_species, other_id, other)
    else:
        for other_id, other in pool.items():
            if other.is_radical:
                add_abstraction(other_id, other, new_id, new_species)

    # dedupe within this batch (symmetric cases can double-generate)
    seen = {}
    for r in out:
        seen[r.key] = r
    return list(seen.values())


def _species_for(sid: str, new_id: str, new_species: Species, pool: dict[str, Species]) -> Species:
    return new_species if sid == new_id else pool[sid]
