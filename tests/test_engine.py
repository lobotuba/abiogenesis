"""Quick sanity checks, runnable without pytest: `python tests/test_engine.py`."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from engine.molecule import ATOMIC_H, MOLECULAR_H2, combine, ethane, methane, propane
from engine.reactions import generate_reactions_for_species
from engine.simulator import SimParams, Simulator
from engine.autocatalysis import find_candidate_cycles


def check(label, cond):
    status = "OK  " if cond else "FAIL"
    print(f"[{status}] {label}")
    assert cond, label


def test_molecule_formulas():
    m = methane()
    check("methane formula CH4", m.formula() == "C1H4")
    e = ethane()
    check("ethane formula C2H6", e.formula() == "C2H6")
    p = propane()
    check("propane formula C3H8", p.formula() == "C3H8")


def test_radicalization_and_combination():
    m = methane()
    sites = m.radicalizable_sites()
    check("methane has 1 site (all 4 H equivalent)", len(sites) == 1)
    check("methane site multiplicity is 4", sites[0].multiplicity == 4)
    methyl = sites[0].product
    check("methyl radical formula CH3*", methyl.formula() == "C1H3•")

    combined = combine(methyl, methyl)
    check("2x methyl -> ethane", combined.formula() == "C2H6")

    saturated = combine(methyl, ATOMIC_H)
    check("methyl + H -> methane", saturated.formula() == "C1H4")

    hh = combine(ATOMIC_H, ATOMIC_H)
    check("H + H -> H2", hh.formula() == "H2")


def test_reaction_generation_from_methane():
    pool = {}
    m = methane()
    reactions = generate_reactions_for_species(
        m.canonical_id(), m, pool,
        max_carbon=6, k_photo=1.0, k_comb=1.0, k_abstr=1.0,
    )
    kinds = {r.kind for r in reactions}
    check("methane alone only yields photolysis reactions", kinds == {"photolysis"})
    check("exactly one photolysis reaction (symmetric H's)", len(reactions) == 1)
    r = reactions[0]
    check("photolysis weight == 4 (equivalent C-H bonds)", r.weight == 4.0)
    prod_formulas = sorted(p.formula() for p in r.products)
    check("products are CH3* and H*", prod_formulas == ["C1H3•", "H•"])


def test_simulation_runs_and_produces_ethane():
    params = SimParams(
        k_photo=0.05, k_comb=5.0, k_abstr=0.5,
        max_carbon=4, t_max=200.0, max_events=5000, sample_every=10, seed=42,
    )
    sim = Simulator(params)
    sim.seed_species(methane(), 200)
    result = sim.run()

    ethane_id = ethane().canonical_id()
    check("simulation produced ethane", result.counts.get(ethane_id, 0) > 0)
    check("simulation produced at least one event", len(result.event_log) > 0)
    check("species pool grew beyond the seed", len(result.species) > 2)

    total_c_seed = 200
    total_c_now = sum(sp.n_carbon * result.counts.get(sid, 0) for sid, sp in result.species.items())
    check(f"carbon count conserved ({total_c_now} == {total_c_seed})", total_c_now == total_c_seed)

    cycles = find_candidate_cycles(result.reactions, result.reaction_fire_counts, result.species)
    print(f"  -> {len(cycles)} candidate autocatalytic/chain cycles found")
    check("cycle search runs without error", isinstance(cycles, list))


if __name__ == "__main__":
    test_molecule_formulas()
    test_radicalization_and_combination()
    test_reaction_generation_from_methane()
    test_simulation_runs_and_produces_ethane()
    print("\nAll sanity checks passed.")
