"""Quick sanity checks, runnable without pytest: `python tests/test_engine.py`."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from engine.molecule import (
    ATMOSPHERIC_GASES,
    ATOMIC_H,
    ATOMIC_N,
    ATOMIC_O,
    HYDROXYL_OH,
    IMIDOGEN_NH,
    MOLECULAR_H2,
    MOLECULAR_H2O,
    MOLECULAR_N2,
    MOLECULAR_NH3,
    MOLECULAR_O2,
    MOLECULAR_O3,
    combine,
    ethane,
    methane,
    propane,
)
from engine.reactions import generate_reactions_for_species
from engine.simulator import SimParams, Simulator
from engine.autocatalysis import find_candidate_cycles
from engine.analysis import chain_length_stats


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


def test_oxygen_chemistry_combinations():
    check("O + O -> O2", combine(ATOMIC_O, ATOMIC_O).formula() == "O2")
    check("O + H -> OH", combine(ATOMIC_O, ATOMIC_H).formula() == "OH•")
    check("H + O -> OH (order independent)", combine(ATOMIC_H, ATOMIC_O).formula() == "OH•")
    check("OH + H -> H2O", combine(HYDROXYL_OH, ATOMIC_H).formula() == "H2O")

    methyl = methane().radicalizable_sites()[0].product
    methoxy = combine(methyl, ATOMIC_O)
    check("CH3* + O* -> CH3O* (alkoxy radical)", methoxy.formula() == "C1H3O•")
    check("alkoxy radical is a radical", methoxy.is_radical)

    methanol = combine(methoxy, ATOMIC_H)
    check("CH3O* + H* -> CH3OH (alcohol)", methanol.formula() == "C1H3OH")
    check("alcohol is closed-shell", not methanol.is_radical)

    from engine.molecule import PeroxyRadical
    peroxy = PeroxyRadical(methyl)
    check("CH3* + O2 -> CH3OO* formula", peroxy.formula() == "C1H3OO•")
    hydroperoxide = combine(peroxy, ATOMIC_H)
    check("CH3OO* + H* -> CH3OOH", hydroperoxide.formula() == "C1H3OOH")
    check("hydroperoxide is closed-shell", not hydroperoxide.is_radical)

    check("two peroxy radicals -> not modeled (None)", combine(peroxy, PeroxyRadical(methyl)) is None)


def test_o2_photolysis_and_scavenging_reaction_generation():
    pool = {"O2": MOLECULAR_O2}
    reactions = generate_reactions_for_species(
        "O2", MOLECULAR_O2, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_photo_o2=2.0, k_photo_co2=0.1, k_o2_scavenge=3.0,
    )
    check("O2 alone yields exactly its photolysis reaction", len(reactions) == 1)
    r = reactions[0]
    check("O2 photolysis kind", r.kind == "photolysis")
    check("O2 photolysis uses k_photo_o2", r.rate_constant == 2.0)
    check("O2 photolysis products are 2x O*", sorted(p.formula() for p in r.products) == ["O•", "O•"])

    # Now a methyl radical discovered while O2 is already in the pool should
    # get an o2_scavenging reaction generated against it.
    methyl = methane().radicalizable_sites()[0].product
    reactions2 = generate_reactions_for_species(
        methyl.canonical_id(), methyl, pool, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_photo_o2=2.0, k_photo_co2=0.1, k_o2_scavenge=3.0,
    )
    scavenging = [r for r in reactions2 if r.kind == "o2_scavenging"]
    check("methyl radical gets an o2_scavenging reaction against existing O2", len(scavenging) == 1)
    check("o2_scavenging product is a peroxy radical", scavenging[0].products[0].formula() == "C1H3OO•")


def test_o2_competes_with_self_combination_for_c2h6():
    """The central question: with O2 present, does CH3* preferentially form
    C2H6 (self-combination) or CH3OO* (O2 scavenging)? Run twice with the
    same seed/rates but very different O2 abundance and check the balance
    of *fired* reactions swings the way real chemistry says it should."""

    def run(o2_count):
        params = SimParams(
            k_photo=0.05, k_comb=5.0, k_abstr=0.5,
            k_photo_o2=0.0,  # isolate the scavenging competition; no O2 photolysis noise
            k_o2_scavenge=5.0,  # same order as k_comb on purpose -- see SimParams docstring
            max_carbon=4, t_max=300.0, max_events=8000, sample_every=50, seed=7,
        )
        sim = Simulator(params)
        sim.seed_species(methane(), 300)
        if o2_count > 0:
            sim.seed_species(ATMOSPHERIC_GASES["O2"](), o2_count)
        result = sim.run()
        comb_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("combination:"))
        scav_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("o2_scavenging:"))
        return comb_fires, scav_fires, result

    comb_no_o2, scav_no_o2, _ = run(0)
    check("with no O2, zero scavenging events fire", scav_no_o2 == 0)
    check("with no O2, some self-combination fires", comb_no_o2 > 0)

    comb_hi_o2, scav_hi_o2, result_hi = run(4000)  # O2-rich: abundant relative to trace CH3*
    check("with abundant O2, scavenging fires far more than self-combination",
          scav_hi_o2 > comb_hi_o2 * 5)
    c2h6_id = ethane().canonical_id()
    print(f"  -> low-O2: {comb_no_o2} combination fires, {scav_no_o2} scavenging fires")
    print(f"  -> high-O2: {comb_hi_o2} combination fires, {scav_hi_o2} scavenging fires, "
          f"final C2H6 count = {result_hi.counts.get(c2h6_id, 0)}")

    # Carbon conservation must still hold across the extended species zoo
    # (peroxy/alkoxy/alcohol/hydroperoxide all delegate n_carbon to their parent).
    total_c = sum(sp.n_carbon * result_hi.counts.get(sid, 0) for sid, sp in result_hi.species.items())
    check(f"carbon conserved with O2 chemistry active ({total_c} == 300)", total_c == 300)


def test_ozone_formation_photolysis_and_scavenging():
    # O* discovered while O2 already present -> ozone_formation reaction.
    pool = {"O2": MOLECULAR_O2}
    reactions = generate_reactions_for_species(
        "O", ATOMIC_O, pool, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0,
        k_o3_formation=4.0, k_photo_o3=0.5, k_o3_scavenge=6.0,
    )
    formation = [r for r in reactions if r.kind == "ozone_formation"]
    check("O* + O2 in pool yields exactly one ozone_formation reaction", len(formation) == 1)
    check("ozone_formation product is O3", formation[0].products[0].formula() == "O3")
    check("ozone_formation uses k_o3_formation", formation[0].rate_constant == 4.0)

    # O3 alone photolyzes back to O2 + O*.
    reactions_o3 = generate_reactions_for_species(
        "O3", MOLECULAR_O3, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_photo_o3=0.5,
    )
    photo = [r for r in reactions_o3 if r.kind == "photolysis"]
    check("O3 alone yields exactly one photolysis reaction", len(photo) == 1)
    check("O3 photolysis uses k_photo_o3", photo[0].rate_constant == 0.5)
    check("O3 photolysis products are O2 + O*",
          sorted(p.formula() for p in photo[0].products) == ["O2", "O•"])

    # A methyl radical discovered while O3 is already present -> ozone_scavenging.
    methyl = methane().radicalizable_sites()[0].product
    reactions3 = generate_reactions_for_species(
        methyl.canonical_id(), methyl, {"O3": MOLECULAR_O3}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_o3_scavenge=6.0,
    )
    scav = [r for r in reactions3 if r.kind == "ozone_scavenging"]
    check("methyl radical gets an ozone_scavenging reaction against existing O3", len(scav) == 1)
    check("ozone_scavenging uses k_o3_scavenge", scav[0].rate_constant == 6.0)
    prod_formulas = sorted(p.formula() for p in scav[0].products)
    check("CH3* + O3 -> CH3O* + O2", prod_formulas == ["C1H3O•", "O2"])


def test_full_atmosphere_ethane_accumulation():
    """End-to-end check with default (nonzero) O2/O3 chemistry active:
    seeding abundant O2 alongside methane should still suppress C2H6
    accumulation once ozone's additional scavenging pathway is included."""
    params = SimParams(
        k_photo=0.05, k_comb=5.0, k_abstr=0.5,
        max_carbon=4, t_max=300.0, max_events=8000, sample_every=50, seed=11,
    )  # k_photo_o2/k_o2_scavenge/k_o3_formation/k_photo_o3/k_o3_scavenge all at their SimParams defaults
    sim = Simulator(params)
    sim.seed_species(methane(), 300)
    sim.seed_species(ATMOSPHERIC_GASES["O2"](), 3000)
    result = sim.run()

    c2h6_id = ethane().canonical_id()
    o2_sink_fires = sum(n for k, n in result.reaction_fire_counts.items()
                         if k.startswith("o2_scavenging:") or k.startswith("ozone_scavenging:"))
    comb_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("combination:"))
    print(f"  -> full atmosphere run: {comb_fires} self-combinations, {o2_sink_fires} O2/O3 scavenging events, "
          f"final C2H6 = {result.counts.get(c2h6_id, 0)}, O3 formed = {result.counts.get('O3', 0)}")
    check("O2/O3 scavenging dominates over self-combination with abundant O2",
          o2_sink_fires > comb_fires)
    check("some ozone actually formed", any(k.startswith("ozone_formation:") for k in result.reaction_fire_counts))


def test_water_photolysis_and_escape_reaction_generation():
    reactions = generate_reactions_for_species(
        "H2O", MOLECULAR_H2O, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_photo_h2o=0.7,
    )
    photo = [r for r in reactions if r.kind == "photolysis"]
    check("H2O alone yields exactly one photolysis reaction", len(photo) == 1)
    check("H2O photolysis uses k_photo_h2o", photo[0].rate_constant == 0.7)
    check("H2O photolysis products are H* + OH*",
          sorted(p.formula() for p in photo[0].products) == ["H•", "OH•"])

    h2_reactions = generate_reactions_for_species(
        "H2", MOLECULAR_H2, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_escape_h2=2.5,
    )
    escape = [r for r in h2_reactions if r.kind == "escape"]
    check("H2 alone yields exactly one escape reaction", len(escape) == 1)
    check("H2 escape uses k_escape_h2", escape[0].rate_constant == 2.5)
    check("H2 escape has no products (irreversible sink)", escape[0].products == ())

    h_reactions = generate_reactions_for_species(
        "H", ATOMIC_H, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_escape_h=1.5,
    )
    escape_h = [r for r in h_reactions if r.kind == "escape"]
    check("H* alone yields exactly one escape reaction", len(escape_h) == 1)
    check("H* escape uses k_escape_h", escape_h[0].rate_constant == 1.5)

    oh_reactions = generate_reactions_for_species(
        "OH", HYDROXYL_OH, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_oh_disprop=3.0,
    )
    disprop = [r for r in oh_reactions if r.kind == "oh_disproportionation"]
    check("OH* alone yields exactly one disproportionation reaction", len(disprop) == 1)
    check("OH disproportionation is a true self-pair (OH + OH)", disprop[0].reactant_ids == ("OH", "OH"))
    check("OH disproportionation uses k_oh_disprop", disprop[0].rate_constant == 3.0)
    check("2 OH* -> H2O + O*", sorted(p.formula() for p in disprop[0].products) == ["H2O", "O•"])


def test_escape_is_an_irreversible_sink():
    """With escape active and nothing else happening, H2 count should just
    drain away with no compensating species appearing anywhere."""
    params = SimParams(k_photo=0.0, k_comb=0.0, k_abstr=0.0, k_photo_h2o=0.0,
                        k_escape_h2=5.0, max_carbon=4, t_max=50.0, max_events=2000,
                        sample_every=10, seed=1)
    sim = Simulator(params)
    sim.seed_species(MOLECULAR_H2, 200)
    result = sim.run()
    check("H2 count dropped via escape", result.counts.get("H2", 0) < 200)
    check("only H2 and its escape reaction exist (no compensating species)",
          set(result.species.keys()) == {"H2"})
    check("escape events actually fired", any(k.startswith("escape:H2") for k in result.reaction_fire_counts))


def test_wet_planet_o2_buildup_needs_hydrogen_escape():
    """The central question: on a wet planet (CH4 + H2O), does UV alone build
    up enough O2/O3 to strangle hydrocarbon chemistry? Real planetary science
    says this requires hydrogen ESCAPE (an irreversible sink) -- without it,
    every O atom pulled off H2O eventually finds its way back via abstraction
    /combination, and the system just cycles. Test both conditions."""

    def run(k_escape):
        params = SimParams(
            k_photo=0.05, k_comb=5.0, k_abstr=0.5, k_photo_h2o=0.05,
            k_escape_h2=k_escape, k_escape_h=k_escape,
            max_carbon=4, t_max=60.0, max_events=800_000, sample_every=1000, seed=13,
        )  # the O2/O3 rate constants stay at their SimParams defaults; k_photo_h2o
        # bumped up from its (deliberately weak, realistic) default so this test
        # converges within a reasonable event budget rather than testing patience
        sim = Simulator(params)
        sim.seed_species(methane(), 300)
        sim.seed_species(MOLECULAR_H2O, 300)
        result = sim.run()
        check(f"run reached t_max (k_escape={k_escape})", result.stopped_reason == "t_max reached")
        o2_now = result.counts.get("O2", 0)
        o3_now = result.counts.get("O3", 0)
        comb_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("combination:"))
        scav_fires = sum(n for k, n in result.reaction_fire_counts.items()
                          if k.startswith("o2_scavenging:") or k.startswith("ozone_scavenging:"))
        return o2_now, o3_now, comb_fires, scav_fires

    o2_closed, o3_closed, comb_closed, scav_closed = run(k_escape=0.0)
    o2_open, o3_open, comb_open, scav_open = run(k_escape=2.0)
    print(f"  -> no escape:     O2={o2_closed} O3={o3_closed} self-comb={comb_closed} scavenging={scav_closed}")
    print(f"  -> with escape:   O2={o2_open} O3={o3_open} self-comb={comb_open} scavenging={scav_open}")

    check("hydrogen escape leads to more free O2 than a closed (no-escape) system",
          o2_open + o3_open > o2_closed + o3_closed)
    check("without escape, O2/O3 scavenging does not dominate hydrocarbon self-combination",
          scav_closed <= comb_closed * 3)


def test_nitrogen_chemistry_combinations():
    check("N* + N* -> N2", combine(ATOMIC_N, ATOMIC_N).formula() == "N2")
    check("N* + H* -> NH* (imidogen)", combine(ATOMIC_N, ATOMIC_H).formula() == "NH•")
    check("H* + N* -> NH* (order independent)", combine(ATOMIC_H, ATOMIC_N).formula() == "NH•")
    nh3 = combine(IMIDOGEN_NH, ATOMIC_H)
    check("NH* + H* -> NH3 (skips aminyl-radical intermediate)", nh3.formula() == "NH3")
    check("NH* + H* product is the shared MOLECULAR_NH3 singleton", nh3.canonical_id() == MOLECULAR_NH3.canonical_id())
    check("NH3 is closed-shell", not MOLECULAR_NH3.is_radical)

    methyl = methane().radicalizable_sites()[0].product
    amine_via_n = combine(methyl, ATOMIC_N)
    check("CH3* + N* -> CH3NH2 (amine)", amine_via_n.formula() == "C1H3NH2")
    check("amine is closed-shell", not amine_via_n.is_radical)
    amine_via_nh = combine(methyl, IMIDOGEN_NH)
    check("CH3* + NH* -> same amine (convergent path)",
          amine_via_nh.canonical_id() == amine_via_n.canonical_id())


def test_discharge_reaction_generation_and_uv_cannot_touch_n2():
    """The core claim behind adding electricity as a variable: N2 should be
    completely unreactive under UV alone (rule 1 never applies to it, and
    no other UV-driven rule touches it either), and only get a live
    discharge reaction when k_discharge_n2 is nonzero."""
    reactions_uv_only = generate_reactions_for_species(
        "N2", MOLECULAR_N2, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0,  # no k_discharge_n2 passed -> defaults to 0.0
    )
    check("N2 with no discharge rate still gets a (rate=0) discharge reaction slot",
          len(reactions_uv_only) == 1 and reactions_uv_only[0].kind == "discharge")
    check("that reaction has zero rate constant (i.e. inert) when discharge is off",
          reactions_uv_only[0].rate_constant == 0.0)

    reactions_discharge = generate_reactions_for_species(
        "N2", MOLECULAR_N2, {}, max_carbon=6,
        k_photo=1.0, k_comb=1.0, k_abstr=1.0, k_discharge_n2=4.0,
    )
    disch = [r for r in reactions_discharge if r.kind == "discharge"]
    check("N2 alone yields exactly one discharge reaction", len(disch) == 1)
    check("discharge reaction uses k_discharge_n2", disch[0].rate_constant == 4.0)
    check("N2 discharge products are 2x N*",
          sorted(p.formula() for p in disch[0].products) == ["N•", "N•"])


def test_electricity_cracks_n2_but_uv_alone_never_does():
    """End-to-end: seed CH4 + N2 together. With discharge off, N2 must stay
    at exactly its seeded count forever, no matter how much UV/time -- UV
    genuinely cannot reach it in this model. With discharge on, N2 gets
    consumed and nitrogen shows up incorporated into hydrocarbon radicals
    as an amine."""

    def run(k_discharge):
        params = SimParams(
            k_photo=0.05, k_comb=5.0, k_abstr=0.5, k_discharge_n2=k_discharge,
            max_carbon=4, t_max=100.0, max_events=200_000, sample_every=1000, seed=17,
        )
        sim = Simulator(params)
        sim.seed_species(methane(), 300)
        sim.seed_species(MOLECULAR_N2, 300)
        return sim.run()

    result_uv_only = run(k_discharge=0.0)
    check("with no discharge, N2 count is completely unchanged by UV",
          result_uv_only.counts.get("N2", 0) == 300)
    check("with no discharge, zero discharge events fired",
          not any(k.startswith("discharge:") for k in result_uv_only.reaction_fire_counts))
    check("with no discharge, no amine species ever appears",
          not any(sid.startswith("RNH2:") for sid in result_uv_only.species))

    result_discharge = run(k_discharge=3.0)
    check("with discharge on, N2 count drops", result_discharge.counts.get("N2", 0) < 300)
    check("with discharge on, discharge events fired",
          any(k.startswith("discharge:") for k in result_discharge.reaction_fire_counts))
    amine_formed = any(sid.startswith("RNH2:") for sid in result_discharge.species)
    print(f"  -> UV-only: N2={result_uv_only.counts.get('N2', 0)}/300 remaining, amine formed=False (by construction)")
    print(f"  -> with discharge: N2={result_discharge.counts.get('N2', 0)}/300 remaining, amine formed={amine_formed}")
    check("with discharge on, at least one amine species (nitrogen incorporated into a hydrocarbon) formed",
          amine_formed)


def test_chain_length_stats_basic():
    from engine.molecule import Molecule as _Molecule

    params = SimParams(k_photo=0.05, k_comb=5.0, k_abstr=0.5,
                        max_carbon=4, t_max=200.0, max_events=5000, sample_every=50, seed=3)
    sim = Simulator(params)
    sim.seed_species(methane(), 200)
    result = sim.run()
    stats = chain_length_stats(result)
    check("chain stats mean_carbon is positive", stats.mean_carbon > 0)
    check("chain stats max_carbon_present <= max_carbon cap", stats.max_carbon_present <= 4)
    expected_total = sum(n for sid, n in result.counts.items()
                          if isinstance(result.species.get(sid), _Molecule) and n > 0)
    check("chain stats counts_by_carbon sums to total hydrocarbon molecule count",
          sum(stats.counts_by_carbon.values()) == expected_total)


def test_higher_uv_yields_longer_chains_at_fixed_time():
    """Regression test for the UV-vs-chain-length finding: at fixed real
    simulated time and fixed starting concentration, more UV should produce
    longer hydrocarbon chains (more radicals in circulation -> radical
    combination, which scales with [R*]^2, becomes relatively more likely
    than non-growing abstraction)."""
    def mean_carbon_after(k_photo):
        params = SimParams(k_photo=k_photo, k_comb=5.0, k_abstr=0.5,
                            max_carbon=6, t_max=5.0, max_events=2_000_000, sample_every=1000, seed=5)
        sim = Simulator(params)
        sim.seed_species(methane(), 500)
        result = sim.run()
        check("run reached t_max (not truncated by max_events)", result.stopped_reason == "t_max reached")
        return chain_length_stats(result).mean_carbon

    low = mean_carbon_after(0.01)
    mid = mean_carbon_after(0.05)
    high = mean_carbon_after(0.2)
    print(f"  -> mean carbon number at k_photo=0.01/0.05/0.2: {low:.2f} / {mid:.2f} / {high:.2f}")
    check("mean chain length increases monotonically with UV level", low < mid < high)


if __name__ == "__main__":
    test_molecule_formulas()
    test_radicalization_and_combination()
    test_reaction_generation_from_methane()
    test_simulation_runs_and_produces_ethane()
    test_oxygen_chemistry_combinations()
    test_o2_photolysis_and_scavenging_reaction_generation()
    test_o2_competes_with_self_combination_for_c2h6()
    test_ozone_formation_photolysis_and_scavenging()
    test_full_atmosphere_ethane_accumulation()
    test_water_photolysis_and_escape_reaction_generation()
    test_escape_is_an_irreversible_sink()
    test_wet_planet_o2_buildup_needs_hydrogen_escape()
    test_nitrogen_chemistry_combinations()
    test_discharge_reaction_generation_and_uv_cannot_touch_n2()
    test_electricity_cracks_n2_but_uv_alone_never_does()
    test_chain_length_stats_basic()
    test_higher_uv_yields_longer_chains_at_fixed_time()
    print("\nAll sanity checks passed.")
