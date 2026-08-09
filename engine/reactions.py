"""
Rule-based reaction generator.

Generic rules, applied to whatever species currently exist in the pool,
produce all candidate elementary reactions:

  1. Photolysis        M + hv -> M* (radical) + H*      (UV C-H homolysis)
  2. Combination         R* + R'* -> R-R'                 (radical-radical)
  3. Abstraction          R* + M-H -> R-H + M*             (H-atom transfer)
  4. O2 scavenging         R* + O2 -> ROO*                   (radical trapped by O2)
  5. Ozone formation        O* + O2 -> O3                      (fixed pair)
  6. Ozone scavenging        R* + O3 -> RO* + O2                (radical trapped by O3)
  7. Escape                   H2 -> (nothing) / H* -> (nothing)   (irreversible loss to space)
  8. OH disproportionation     OH* + OH* -> H2O + O*                (fixed self-pair)
  9. Electric discharge         N2 + spark -> 2 N*                    (fixed pair; NOT UV)

Rule 9 is the point of the whole "electricity" knob: N2's triple bond
(~9.8 eV) is well beyond what UV drives in this model (C-H ~4.5 eV, O2's
O=O ~5.2 eV) -- no photolysis rate constant, however high, ever touches
N2. A separate `k_discharge_n2` rate constant represents a qualitatively
different energy source (electric spark, as in Miller-Urey) that *can*
reach that bond. Once N* exists, it plugs into the *existing* generic
combination/abstraction rules for free (nothing N-specific needed there):
N* + N* -> N2 (recombination), N* + H* -> NH* (imidogen), N* + hydrocarbon
radical -> a closed-shell amine (R-NH2, skipping the aminyl-radical
intermediate a fully valence-rigorous treatment would need -- N's valence
of 3 doesn't split as cleanly as O's valence of 2 into "one bond + one
radical site", so amine formation is a deliberate one-step simplification;
see engine/molecule.py). This is a minimal, focused addition -- no N-O
cross chemistry (real NOx chemistry) is modeled.

H2O also gets a fixed unimolecular photolysis reaction (H2O -> H* + OH*,
rule 1c below), the real mechanism by which a wet planet's own atmospheric
water can be a *source* of O atoms/O2/O3 with no biology or hydrocarbons
involved at all. On its own, though, H2O photolysis's *primary* products
(H* and OH*) have no path to free O2/O3 -- something has to convert OH*
into an O atom first. Rule 8 is that missing link: two OH radicals meeting
and disproportionating into H2O + O* is a real, well-known secondary step,
and it's what lets water photolysis bootstrap O2/O3 production without any
O2 seeded beforehand. (This model doesn't go further into full HOx
chemistry -- no H2O2/HO2 -- rule 8 is the one secondary step needed to
answer "can a wet planet make its own O2 from UV alone.")

Rule 7 (escape) is what makes *sustained* O2/O3 buildup possible: without
some one-way loss of hydrogen, every O atom liberated from H2O eventually
finds its way back into H2O via the same radical chemistry (abstraction,
combination), and the system just cycles at a steady state instead of
trending toward net oxidation. Real planetary atmospheres do lose hydrogen
to space (it's the lightest gas), which is the accepted mechanism for
abiotic O2/O3 buildup on a wet planet with no life -- studied as a "false
positive" biosignature scenario for exoplanets. Escape is modeled as a
simple first-order sink (no altitude/exobase physics) so the sweep tool can
ask directly: at what escape rate does O2/O3 buildup start to dominate
hydrocarbon chemistry?

Rule (3) is what lets the network self-amplify: a radical can be regenerated
by a different chain a few steps later, which is the chemically real analog
of autocatalysis this project is built around. Rules (4) and (6) are what
let O2 and O3 *compete* with rule 2 for the same radical pool -- whether a
hydrocarbon radical dimerizes (rule 2, e.g. CH3* + CH3* -> C2H6) or gets
scavenged (rule 4, CH3* + O2 -> CH3OO*; rule 6, CH3* + O3 -> CH3O* + O2) is
decided purely by relative propensities, i.e. by how much O2/O3 is actually
present -- not by a thumb on the scale.

O2 and CO2 each get one fixed unimolecular photolysis reaction (O2 -> 2 O*,
CO2 -> CO + O*) generated directly (rule 1b below) rather than through rule
1, since they have no C-H bond to homolyze. O3 does too (O3 -> O2 + O*).
Ozone itself forms from O* + O2 (rule 5) -- a simplified bimolecular stand-in
for the real termolecular O + O2 + M -> O3 + M (no third-body/pressure
dependence modeled) -- and, being "highly reactive" as ozone actually is,
attacks hydrocarbon radicals directly (rule 6) rather than sitting inert.

Rates are deliberately not calibrated to real photon flux / Arrhenius
kinetics -- they're relative, adjustable knobs (see engine/simulator.py)
meant for exploring qualitative network behavior, not for quantitative
photochemistry.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import comb
from typing import Literal

from .molecule import (
    ATOMIC_H,
    ATOMIC_N,
    ATOMIC_O,
    AlkoxyRadical,
    HYDROXYL_OH,
    MOLECULAR_CO,
    MOLECULAR_H2O,
    MOLECULAR_O2,
    MOLECULAR_O3,
    Molecule,
    PeroxyRadical,
    Species,
    combine,
)

ReactionKind = Literal[
    "photolysis", "combination", "abstraction",
    "o2_scavenging", "ozone_formation", "ozone_scavenging", "escape",
    "oh_disproportionation", "discharge",
]


def _peroxy(radical: Molecule) -> PeroxyRadical:
    return PeroxyRadical(radical)


def _alkoxy(radical: Molecule) -> AlkoxyRadical:
    return AlkoxyRadical(radical)

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
        # General Gillespie combinatorial term: for each distinct reactant
        # id appearing with multiplicity m in this reaction, the number of
        # ways to draw m molecules of it from the current pool is
        # comb(n, m) (n choose m; 0 if n < m). Multiplying across distinct
        # ids handles unimolecular, bimolecular (both same and different
        # species -- the two cases this used to hand-special-case), and
        # genuinely termolecular reactions (e.g. polymer.py's templated
        # ligation, A + A + T -> T + T) with one formula. For every
        # reactant_ids shape used before this generalization it reduces to
        # exactly the old na / na*nb / na*(na-1)/2 arithmetic.
        factor = 1.0
        for rid, multiplicity in Counter(self.reactant_ids).items():
            factor *= comb(counts.get(rid, 0), multiplicity)
        return self.rate_constant * self.weight * factor


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
    k_photo_o2: float = 0.0,
    k_photo_co2: float = 0.0,
    k_o2_scavenge: float = 0.0,
    k_o3_formation: float = 0.0,
    k_photo_o3: float = 0.0,
    k_o3_scavenge: float = 0.0,
    k_photo_h2o: float = 0.0,
    k_escape_h2: float = 0.0,
    k_escape_h: float = 0.0,
    k_oh_disprop: float = 0.0,
    k_discharge_n2: float = 0.0,
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

    # 1b. Fixed unimolecular photolysis for O2 and CO2 -- neither has a C-H
    #     bond to homolyze, so rule 1 above never fires for them.
    if new_id == "O2":
        products = (ATOMIC_O, ATOMIC_O)
        out.append(Reaction(
            kind="photolysis", reactant_ids=(new_id,), products=products,
            rate_constant=k_photo_o2, weight=1.0,
            key=_mk_key("photolysis", (new_id,), tuple(p.canonical_id() for p in products)),
        ))
    if new_id == "CO2":
        products = (MOLECULAR_CO, ATOMIC_O)
        out.append(Reaction(
            kind="photolysis", reactant_ids=(new_id,), products=products,
            rate_constant=k_photo_co2, weight=1.0,
            key=_mk_key("photolysis", (new_id,), tuple(p.canonical_id() for p in products)),
        ))
    if new_id == "O3":
        products = (MOLECULAR_O2, ATOMIC_O)
        out.append(Reaction(
            kind="photolysis", reactant_ids=(new_id,), products=products,
            rate_constant=k_photo_o3, weight=1.0,
            key=_mk_key("photolysis", (new_id,), tuple(p.canonical_id() for p in products)),
        ))
    if new_id == "H2O":
        products = (ATOMIC_H, HYDROXYL_OH)
        out.append(Reaction(
            kind="photolysis", reactant_ids=(new_id,), products=products,
            rate_constant=k_photo_h2o, weight=1.0,
            key=_mk_key("photolysis", (new_id,), tuple(p.canonical_id() for p in products)),
        ))

    # 8. OH disproportionation: OH* + OH* -> H2O + O*. A fixed self-pair
    #    (both reactants are the same species), generated once when OH is
    #    first discovered -- no bidirectional pool-checking needed since
    #    there's no second species to wait on.
    if new_id == "OH":
        reactants = ("OH", "OH")
        products = (MOLECULAR_H2O, ATOMIC_O)
        out.append(Reaction(
            kind="oh_disproportionation", reactant_ids=reactants, products=products,
            rate_constant=k_oh_disprop, weight=1.0,
            key=_mk_key("oh_disproportionation", reactants, tuple(p.canonical_id() for p in products)),
        ))

    # 9. Electric discharge: N2 + spark -> 2 N*. Not a photolysis reaction
    #    (deliberately not folded into rule 1/1b) -- it represents a
    #    qualitatively different energy source (electricity, not UV) that
    #    can reach N2's triple bond when k_discharge_n2 is turned on.
    if new_id == "N2":
        products = (ATOMIC_N, ATOMIC_N)
        out.append(Reaction(
            kind="discharge", reactant_ids=(new_id,), products=products,
            rate_constant=k_discharge_n2, weight=1.0,
            key=_mk_key("discharge", (new_id,), tuple(p.canonical_id() for p in products)),
        ))

    # 7. Escape: H2 and H* are the lightest species in the system and the
    #    ones real planetary atmospheres actually lose to space. Modeled as
    #    unimolecular decay to nothing (empty products) -- an irreversible
    #    sink, unlike every other reaction here which just reshuffles atoms
    #    among tracked species.
    if new_id == "H2":
        out.append(Reaction(
            kind="escape", reactant_ids=("H2",), products=(),
            rate_constant=k_escape_h2, weight=1.0, key=_mk_key("escape", ("H2",), ()),
        ))
    if new_id == "H":
        out.append(Reaction(
            kind="escape", reactant_ids=("H",), products=(),
            rate_constant=k_escape_h, weight=1.0, key=_mk_key("escape", ("H",), ()),
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
            if product is None or product.n_carbon > max_carbon:
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

    # 2b. O2 scavenging: any hydrocarbon radical + O2 -> peroxy radical. This
    #     is the pathway that competes with rule 2's radical self-combination
    #     for the same CH3*-type radical pool. Covers both directions: a new
    #     hydrocarbon radical appearing while O2 is already present, and O2
    #     appearing while hydrocarbon radicals already exist.
    if new_species.is_radical and isinstance(new_species, Molecule) and "O2" in pool:
        product = _peroxy(new_species)
        if product.n_carbon <= max_carbon:
            reactants = tuple(sorted((new_id, "O2")))
            out.append(Reaction(
                kind="o2_scavenging", reactant_ids=reactants, products=(product,),
                rate_constant=k_o2_scavenge, weight=1.0,
                key=_mk_key("o2_scavenging", reactants, (product.canonical_id(),)),
            ))
    if new_id == "O2":
        for other_id, other in pool.items():
            if other.is_radical and isinstance(other, Molecule):
                product = _peroxy(other)
                if product.n_carbon > max_carbon:
                    continue
                reactants = tuple(sorted((other_id, "O2")))
                out.append(Reaction(
                    kind="o2_scavenging", reactant_ids=reactants, products=(product,),
                    rate_constant=k_o2_scavenge, weight=1.0,
                    key=_mk_key("o2_scavenging", reactants, (product.canonical_id(),)),
                ))

    # 5. Ozone formation: O* + O2 -> O3. A fixed pair of specific species
    #    (not a generic radical-pool loop like rules 2/2b), so just check
    #    both discovery orders directly.
    if new_id == "O" and "O2" in pool:
        reactants = ("O", "O2")
        out.append(Reaction(
            kind="ozone_formation", reactant_ids=reactants, products=(MOLECULAR_O3,),
            rate_constant=k_o3_formation, weight=1.0,
            key=_mk_key("ozone_formation", reactants, ("O3",)),
        ))
    if new_id == "O2" and "O" in pool:
        reactants = ("O", "O2")
        out.append(Reaction(
            kind="ozone_formation", reactant_ids=reactants, products=(MOLECULAR_O3,),
            rate_constant=k_o3_formation, weight=1.0,
            key=_mk_key("ozone_formation", reactants, ("O3",)),
        ))

    # 6. Ozone scavenging: any hydrocarbon radical + O3 -> alkoxy radical +
    #    O2 (O3 is "highly reactive" -- this is the fast, real R* + O3 ->
    #    RO* + O2 step, a second radical sink alongside rule 2b, and it
    #    regenerates O2 in the process). Same bidirectional-discovery pattern
    #    as O2 scavenging.
    if new_species.is_radical and isinstance(new_species, Molecule) and "O3" in pool:
        product = _alkoxy(new_species)
        if product.n_carbon <= max_carbon:
            reactants = tuple(sorted((new_id, "O3")))
            out.append(Reaction(
                kind="ozone_scavenging", reactant_ids=reactants,
                products=(product, MOLECULAR_O2),
                rate_constant=k_o3_scavenge, weight=1.0,
                key=_mk_key("ozone_scavenging", reactants, (product.canonical_id(), "O2")),
            ))
    if new_id == "O3":
        for other_id, other in pool.items():
            if other.is_radical and isinstance(other, Molecule):
                product = _alkoxy(other)
                if product.n_carbon > max_carbon:
                    continue
                reactants = tuple(sorted((other_id, "O3")))
                out.append(Reaction(
                    kind="ozone_scavenging", reactant_ids=reactants,
                    products=(product, MOLECULAR_O2),
                    rate_constant=k_o3_scavenge, weight=1.0,
                    key=_mk_key("ozone_scavenging", reactants, (product.canonical_id(), "O2")),
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
