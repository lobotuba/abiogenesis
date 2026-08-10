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
from engine.reactions import Reaction, generate_reactions_for_species
from engine.simulator import SimParams, Simulator
from engine.autocatalysis import find_candidate_cycles
from engine.analysis import chain_length_stats
from engine.formose import PENTOSE_DIASTEREOMERS, FormoseParams, build_formose_reactions, run_formose
from engine.nucleotide import (
    ANHYDRONUCLEOSIDE,
    CYANAMIDE,
    CYANOACETYLENE,
    GLYCERALDEHYDE,
    GLYCOLALDEHYDE,
    PHOSPHATE,
    RIBONUCLEOTIDE,
    NucleotideParams,
    build_nucleotide_reactions,
    run_nucleotide,
)
from engine.polymer import (
    DUPLEX,
    FRAGMENT,
    TEMPLATE,
    Duplex,
    Oligomer,
    PolymerParams,
    build_polymer_reactions,
    run_polymer,
)
from engine.selection import SelectionParams, build_selection_reactions, run_selection
from engine.hypercycle import HypercycleParams, build_hypercycle_reactions, run_hypercycle


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


def test_electricity_indirectly_enables_o2_via_hydrogen_sink():
    """Follow-up finding: discharge doesn't touch water/O2 chemistry at all
    (it only cracks N2), but nitrogen fixation consumes H* (into NH3/amines),
    starving the H* + OH* -> H2O recombination that normally keeps a closed
    system oxidant-free. So a CLOSED system (no hydrogen escape) should still
    show more O2+O3 with discharge on than off, purely as a side effect --
    and that carbon should end up in three clean, mutually exclusive buckets:
    plain hydrocarbons, amines, and O-containing oxidation products."""
    from engine.molecule import Molecule as _Molecule

    def run(k_discharge):
        params = SimParams(
            k_photo=0.05, k_comb=5.0, k_abstr=0.5, k_discharge_n2=k_discharge,
            # no k_escape_h2/k_escape_h at all -- closed system on purpose
            max_carbon=4, t_max=100.0, max_events=1_000_000, sample_every=2000, seed=23,
        )
        sim = Simulator(params)
        sim.seed_species(methane(), 300)
        sim.seed_species(MOLECULAR_H2O, 300)
        sim.seed_species(MOLECULAR_N2, 300)
        return sim.run()

    result_off = run(0.0)
    result_on = run(3.0)
    o2o3_off = result_off.counts.get("O2", 0) + result_off.counts.get("O3", 0)
    o2o3_on = result_on.counts.get("O2", 0) + result_on.counts.get("O3", 0)
    print(f"  -> closed system, discharge off: O2+O3={o2o3_off}")
    print(f"  -> closed system, discharge on:  O2+O3={o2o3_on}")
    check("discharge off, closed system: zero free O2/O3 (matches the no-N2 wet-planet result)",
          o2o3_off == 0)
    check("discharge on, closed system: some free O2/O3 appears purely as a hydrogen-sink side effect",
          o2o3_on > 0)

    # Carbon accounting: every carbon atom is in exactly one of three buckets.
    total_carbon = sum(sp.n_carbon * result_on.counts.get(sid, 0) for sid, sp in result_on.species.items())
    hydrocarbon_carbon = sum(sp.n_carbon * result_on.counts.get(sid, 0) for sid, sp in result_on.species.items()
                              if isinstance(sp, _Molecule))
    amine_carbon = sum(sp.n_carbon * result_on.counts.get(sid, 0) for sid, sp in result_on.species.items()
                        if sid.startswith("RNH2:"))
    oxidized_carbon = sum(sp.n_carbon * result_on.counts.get(sid, 0) for sid, sp in result_on.species.items()
                           if sid.startswith(("ROO:", "RO:", "ROOH:", "ROH:")))
    check(f"hydrocarbon + amine + oxidized carbon accounts for all {total_carbon} carbon atoms",
          hydrocarbon_carbon + amine_carbon + oxidized_carbon == total_carbon)
    check("with discharge on, nitrogen fixation captured the large majority of carbon",
          amine_carbon > total_carbon * 0.5)


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


def test_formose_reaction_structure():
    params = FormoseParams(max_sugar_carbon=5, k_stabilize=0.0)
    reactions = build_formose_reactions(params)
    kinds = {r.kind for r in reactions}
    check("formose network has init/cannizzaro/aldol/retro_aldol kinds",
          kinds == {"formose_init", "cannizzaro", "aldol", "retro_aldol"})
    check("no stabilization reaction when k_stabilize is 0",
          not any(r.kind == "stabilization" for r in reactions))
    # max_sugar_carbon=5 -> aldol/retro_aldol pairs for n=2,3,4 (growing to 3,4,5)
    aldol = [r for r in reactions if r.kind == "aldol"]
    check("aldol reactions generated for each growth step up to the carbon ceiling", len(aldol) == 3)

    params_stab = FormoseParams(max_sugar_carbon=5, k_stabilize=1.0, stabilize_carbon=5)
    reactions_stab = build_formose_reactions(params_stab)
    stab = [r for r in reactions_stab if r.kind == "stabilization"]
    check("stabilization reaction appears when k_stabilize > 0", len(stab) == 1)
    check("stabilization reaction targets C5", stab[0].reactant_ids == ("C5sugar",))
    check("stabilization product is the stabilized C5 species", stab[0].products[0].formula() == "(CH2O)5·borate")


def test_formose_carbon_conservation_and_chain_growth():
    params = FormoseParams(k_init=0.05, k_aldol=1.0, k_retro_aldol=0.3, k_cannizzaro=0.05,
                            k_stabilize=0.0, max_sugar_carbon=7, t_max=100.0, max_events=100_000,
                            sample_every=500, seed=11)
    result = run_formose(params, initial_hcho=500)

    total_now = (
        result.counts.get("HCHO", 0)
        + 2 * result.counts.get("waste", 0)
        + sum(n * result.counts.get(f"C{n}sugar", 0) for n in range(2, params.max_sugar_carbon + 1))
    )
    check(f"carbon conserved through formose chemistry ({total_now} == 500)", total_now == 500)

    reached = {n for n in range(2, params.max_sugar_carbon + 1) if result.counts.get(f"C{n}sugar", 0) > 0}
    print(f"  -> sugar sizes with nonzero count at end of run: {sorted(reached)}")
    check("aldol addition actually grows sugars past glycolaldehyde (reaches C4+ at some point)",
          any(n >= 4 for n in reached) or any(
              any(k.startswith(f"aldol:C{n}sugar") for k in result.reaction_fire_counts) for n in range(3, 6)
          ))


def test_formose_ribose_needs_stabilization():
    """The central question: does a persistent C5 (ribose-sized) sugar
    pool require mineral stabilization, or does formose chemistry alone
    let it accumulate? Real formose is dominated by retro-aldol scrambling
    (the "sugar problem") -- expect the unprotected C5 pool to stay small/
    transient while the stabilized pool (once protected) only grows."""

    def run(k_stabilize):
        params = FormoseParams(
            k_init=0.05, k_aldol=1.0, k_retro_aldol=0.3, k_cannizzaro=0.05,
            k_stabilize=k_stabilize, stabilize_carbon=5, max_sugar_carbon=8,
            t_max=300.0, max_events=300_000, sample_every=1000, seed=7,
        )
        return run_formose(params, initial_hcho=1000)

    result_off = run(0.0)
    result_on = run(0.5)

    c5_off = result_off.counts.get("C5sugar", 0)
    stab_off = result_off.counts.get("C5sugar-stabilized", 0)
    c5_on = result_on.counts.get("C5sugar", 0)
    stab_on = result_on.counts.get("C5sugar-stabilized", 0)

    print(f"  -> no stabilization: free C5={c5_off}, stabilized C5={stab_off}")
    print(f"  -> with stabilization: free C5={c5_on}, stabilized C5={stab_on}")

    check("without stabilization, no stabilized C5 sugar exists (mechanism is off)", stab_off == 0)
    check("with stabilization on, a protected C5 pool accumulates", stab_on > 0)
    check("stabilization captures more C5-equivalent sugar than the free pool ever holds without it",
          stab_on > c5_off)


def test_formose_stereoisomer_tracking_structure():
    params = FormoseParams(max_sugar_carbon=7, k_stabilize=1.0, stabilize_carbon=5,
                            track_pentose_stereoisomers=True, ribose_selectivity=3.0)
    reactions = build_formose_reactions(params)

    aldol_into_c5 = [r for r in reactions if r.kind == "aldol" and r.reactant_ids == ("C4sugar", "HCHO")]
    check("C4 -> C5 aldol growth splits into all 4 pentose diastereomers", len(aldol_into_c5) == 4)
    check("each variant produced by aldol growth into C5 is a real named diastereomer",
          {p.canonical_id() for r in aldol_into_c5 for p in r.products}
          == {f"C5sugar-{v}" for v in PENTOSE_DIASTEREOMERS})

    stab = [r for r in reactions if r.kind == "stabilization"]
    check("stabilization reaction exists for each of the 4 diastereomers", len(stab) == 4)
    ribose_stab = next(r for r in stab if "ribose" in r.reactant_ids[0])
    other_stab = [r for r in stab if "ribose" not in r.reactant_ids[0]]
    check("ribose's stabilization rate is boosted by ribose_selectivity", ribose_stab.rate_constant == 3.0)
    check("the other 3 diastereomers keep the base stabilization rate", all(r.rate_constant == 1.0 for r in other_stab))

    # without stereoisomer tracking, the network should be identical to the plain generic model
    params_generic = FormoseParams(max_sugar_carbon=7, k_stabilize=1.0, stabilize_carbon=5,
                                    track_pentose_stereoisomers=False)
    reactions_generic = build_formose_reactions(params_generic)
    check("stereoisomer tracking off falls back to a single generic C5 stabilization reaction",
          len([r for r in reactions_generic if r.kind == "stabilization"]) == 1)


def test_formose_ribose_selectivity_is_a_real_correction_factor():
    """Answers the actual question this feature exists for: is there a
    variable that lets ribose be preferentially captured, rather than
    treating the whole C5 pool as generic? Real chemistry: aldol addition
    itself is unselective (no chiral catalyst), but borate's binding
    affinity is NOT uniform across pentose stereoisomers -- that's what
    ribose_selectivity represents, applied only to the stabilization step."""

    def stabilized_fractions(selectivity):
        params = FormoseParams(
            k_init=0.05, k_aldol=1.0, k_retro_aldol=0.3, k_cannizzaro=0.05,
            k_stabilize=0.5, stabilize_carbon=5, max_sugar_carbon=8,
            track_pentose_stereoisomers=True, ribose_selectivity=selectivity,
            t_max=300.0, max_events=300_000, sample_every=1000, seed=7,
        )
        result = run_formose(params, initial_hcho=1000)
        counts = {v: result.counts.get(f"C5sugar-stabilized-{v}", 0) for v in PENTOSE_DIASTEREOMERS}
        total = sum(counts.values())
        return counts, (counts["ribose"] / total if total else 0.0)

    counts_unselective, frac_unselective = stabilized_fractions(1.0)
    counts_selective, frac_selective = stabilized_fractions(5.0)
    print(f"  -> selectivity=1.0: {counts_unselective} (ribose fraction {frac_unselective:.1%})")
    print(f"  -> selectivity=5.0: {counts_selective} (ribose fraction {frac_selective:.1%})")

    check("with no selectivity, ribose's share of the stabilized pool is roughly the naive 1/4 baseline",
          0.10 < frac_unselective < 0.45)
    check("boosting ribose_selectivity meaningfully raises ribose's share of the stabilized pool",
          frac_selective > frac_unselective)
    check("boosted selectivity produces more absolute ribose than the unselective case",
          counts_selective["ribose"] > counts_unselective["ribose"])


def test_nucleotide_reaction_structure():
    reactions = build_nucleotide_reactions(NucleotideParams(k_photo_destroy=1.0))
    kinds = {r.kind for r in reactions}
    check("nucleotide network has all expected reaction kinds", kinds == {
        "oxazole_formation", "aminooxazoline_formation", "photo_selection",
        "anhydro_formation", "offpathway_consumption", "phosphorylation",
    })
    check("aminooxazoline formation splits into all 4 diastereomers",
          len([r for r in reactions if r.kind == "aminooxazoline_formation"]) == 4)
    check("photo_selection only touches the 3 non-ribo diastereomers (ribo is photostable)",
          {r.reactant_ids[0] for r in reactions if r.kind == "photo_selection"}
          == {"aminooxazoline-arabinose", "aminooxazoline-xylose", "aminooxazoline-lyxose"})
    check("only ribo-aminooxazoline reacts productively with cyanoacetylene",
          [r for r in reactions if r.kind == "anhydro_formation"][0].reactant_ids
          == tuple(sorted(("aminooxazoline-ribose", "cyanoacetylene"))))
    check("the 3 non-ribo diastereomers also consume cyanoacetylene, but unproductively",
          len([r for r in reactions if r.kind == "offpathway_consumption"]) == 3)

    reactions_off = build_nucleotide_reactions(NucleotideParams(k_photo_destroy=0.0))
    check("photo_selection reactions still exist at k_photo_destroy=0 (present but inert)",
          len([r for r in reactions_off if r.kind == "photo_selection"]) == 3)
    check("k_photo_destroy=0 means zero propensity for photo_selection",
          all(r.rate_constant == 0.0 for r in reactions_off if r.kind == "photo_selection"))


def test_nucleotide_forms_and_photoselection_protects_scarce_cyanoacetylene():
    """The central question: does the photochemical selection step actually
    matter for nucleotide yield? It's tempting to assume selection helps by
    recycling material back to reactants, but it doesn't (photolysis isn't
    a retro-reaction) -- what it actually does is stop the 3 non-ribo
    diastereomers from wastefully consuming cyanoacetylene, which only
    shows up when cyanoacetylene is the scarce reagent (seeded lower than
    the aminooxazoline pool it competes against) and destruction is fast
    relative to that competing bimolecular consumption."""

    def run(k_photo_destroy):
        params = NucleotideParams(
            k_oxazole=1.0, k_aminooxazoline=1.0, k_photo_destroy=k_photo_destroy,
            k_anhydro=1.0, k_phosphorylate=1.0, t_max=50.0, max_events=200_000,
            sample_every=1000, seed=11,
        )
        init = {
            GLYCOLALDEHYDE.canonical_id(): 2000,
            GLYCERALDEHYDE.canonical_id(): 2000,
            CYANAMIDE.canonical_id(): 2000,
            CYANOACETYLENE.canonical_id(): 500,  # deliberately scarce relative to the aminooxazoline pool
            PHOSPHATE.canonical_id(): 2000,
        }
        return run_nucleotide(params, init)

    result_off = run(0.0)
    result_on = run(2000.0)  # must be fast relative to k_anhydro * initial [cyanoacetylene] (~500) to win
    nuc_off = result_off.counts.get(RIBONUCLEOTIDE.canonical_id(), 0)
    nuc_on = result_on.counts.get(RIBONUCLEOTIDE.canonical_id(), 0)
    print(f"  -> no selection: ribonucleotide={nuc_off}")
    print(f"  -> selection on: ribonucleotide={nuc_on}")

    check("some nucleotide forms even without selection (the pathway works at all)", nuc_off > 0)
    check("fast photochemical selection substantially increases nucleotide yield "
          "by protecting scarce cyanoacetylene from the unproductive diastereomers",
          nuc_on > nuc_off * 1.5)


def test_reaction_propensity_generalizes_to_n_body():
    """Reaction.propensity used to hand-special-case 1 vs 2 reactants. It
    was generalized to a Counter-based comb(n, multiplicity) product to
    support polymer.py's genuine 3-reactant templated-ligation step. This
    checks the generalization reduces to the exact old arithmetic for
    unimolecular and both bimolecular shapes, and gets termolecular cases
    right too (both all-distinct and repeated-species)."""
    unimolecular = Reaction(kind="x", reactant_ids=("A",), products=(), rate_constant=2.0, weight=1.0, key="u")
    check("unimolecular: k * weight * n", unimolecular.propensity({"A": 5}) == 2.0 * 5)

    bimolecular_distinct = Reaction(kind="x", reactant_ids=("A", "B"), products=(), rate_constant=2.0, weight=1.0, key="bd")
    check("bimolecular distinct species: k * weight * nA * nB",
          bimolecular_distinct.propensity({"A": 5, "B": 3}) == 2.0 * 5 * 3)

    bimolecular_same = Reaction(kind="x", reactant_ids=("A", "A"), products=(), rate_constant=2.0, weight=1.0, key="bs")
    check("bimolecular same species: k * weight * nA*(nA-1)/2",
          bimolecular_same.propensity({"A": 5}) == 2.0 * (5 * 4 / 2.0))
    check("bimolecular same species with only 1 present: zero propensity (can't pick a pair)",
          bimolecular_same.propensity({"A": 1}) == 0.0)

    termolecular_distinct = Reaction(kind="x", reactant_ids=("A", "B", "C"), products=(), rate_constant=1.0, weight=1.0, key="td")
    check("termolecular, 3 distinct species: k * weight * nA * nB * nC",
          termolecular_distinct.propensity({"A": 5, "B": 3, "C": 2}) == 1.0 * 5 * 3 * 2)

    termolecular_repeated = Reaction(kind="x", reactant_ids=("A", "A", "B"), products=(), rate_constant=1.0, weight=1.0, key="tr")
    check("termolecular, 2 copies of A + 1 B: k * weight * comb(nA,2) * nB",
          termolecular_repeated.propensity({"A": 5, "B": 3}) == 1.0 * (5 * 4 / 2.0) * 3)
    check("termolecular, 2 copies of A but only 1 present: zero propensity",
          termolecular_repeated.propensity({"A": 1, "B": 3}) == 0.0)


def test_polymer_reaction_structure():
    reactions = build_polymer_reactions(PolymerParams(k_template=1.0, k_duplex=1.0, k_melt=1.0))
    kinds = {r.kind for r in reactions}
    check("polymer network has all expected reaction kinds", kinds == {
        "oligomerization", "background_ligation", "templated_ligation",
        "duplex_formation", "duplex_melting",
    })

    oligo = next(r for r in reactions if r.kind == "oligomerization")
    check("oligomerization consumes 3 ribonucleotides", oligo.reactant_ids == (RIBONUCLEOTIDE.canonical_id(),) * 3)
    check("oligomerization produces one FRAGMENT", oligo.products == (FRAGMENT,))

    templated = next(r for r in reactions if r.kind == "templated_ligation")
    check("templated_ligation is genuinely termolecular: 2 FRAGMENT + 1 TEMPLATE",
          sorted(templated.reactant_ids) == sorted((FRAGMENT.canonical_id(), FRAGMENT.canonical_id(), TEMPLATE.canonical_id())))
    check("templated_ligation produces 2 TEMPLATE (the existing one regenerated, plus a new copy)",
          templated.products == (TEMPLATE, TEMPLATE))

    duplex_form = next(r for r in reactions if r.kind == "duplex_formation")
    check("duplex_formation consumes 2 TEMPLATE", duplex_form.reactant_ids == (TEMPLATE.canonical_id(), TEMPLATE.canonical_id()))
    check("duplex_formation produces DUPLEX", duplex_form.products == (DUPLEX,))

    melt = next(r for r in reactions if r.kind == "duplex_melting")
    check("duplex_melting reverses duplex_formation", melt.reactant_ids == (DUPLEX.canonical_id(),) and melt.products == (TEMPLATE, TEMPLATE))


def test_polymer_mass_conservation():
    """Every reaction in this module is a pure reorganization of
    ribonucleotide units (1 per RIBONUCLEOTIDE, 3 per FRAGMENT, 6 per
    TEMPLATE, 12 per DUPLEX) -- nothing is created or destroyed. That
    conserved quantity should equal the initial seeded ribonucleotide count
    at every sampled point in the run, the same correctness technique used
    for carbon accounting elsewhere in this project."""
    params = PolymerParams(k_oligomerize=0.02, k_background=0.02, k_template=1.0,
                            k_duplex=0.5, k_melt=0.05, t_max=200.0, max_events=200_000,
                            sample_every=500, seed=5)
    initial_ribo = 3000
    result = run_polymer(params, {RIBONUCLEOTIDE.canonical_id(): initial_ribo, TEMPLATE.canonical_id(): 2})
    initial_units = initial_ribo + 2 * 6  # the 2 seed TEMPLATE count as 6 units each

    for snap in result.history_counts:
        total = (snap.get(RIBONUCLEOTIDE.canonical_id(), 0) * 1
                 + snap.get(FRAGMENT.canonical_id(), 0) * 3
                 + snap.get(TEMPLATE.canonical_id(), 0) * 6
                 + snap.get(DUPLEX.canonical_id(), 0) * 12)
        check(f"ribonucleotide-unit total conserved at every sampled point (got {total}, want {initial_units})",
              total == initial_units)
    print(f"  -> conserved at {len(result.history_counts)} sampled points, final state: {result.counts}")


def test_polymer_templating_accelerates_strand_formation():
    """The core autocatalysis question: does adding the templated pathway
    actually make TEMPLATE accumulate faster than background ligation
    alone, within a fixed time window? (Given infinite time a closed
    fragment pool eventually fully converts either way -- the real
    signature of autocatalysis is rate, not final yield, so this compares
    at a fixed t_max rather than waiting for exhaustion.)"""

    def run(k_template):
        params = PolymerParams(k_oligomerize=0.0, k_background=0.0002, k_template=k_template,
                                k_duplex=0.0, k_melt=0.0, t_max=10.0, max_events=200_000,
                                sample_every=1000, seed=9)
        init = {FRAGMENT.canonical_id(): 200, TEMPLATE.canonical_id(): 2}
        return run_polymer(params, init)

    result_off = run(0.0)
    result_on = run(2.0)
    tmpl_off = result_off.counts.get(TEMPLATE.canonical_id(), 0)
    tmpl_on = result_on.counts.get(TEMPLATE.canonical_id(), 0)
    print(f"  -> background ligation only: template={tmpl_off}")
    print(f"  -> templated (autocatalytic) pathway on: template={tmpl_on}")

    check("background ligation alone still produces some template (the pathway works without templating)",
          tmpl_off > 0)
    check("the templated pathway accumulates template far faster within the same time window",
          tmpl_on > tmpl_off * 2)


def test_polymer_duplex_formation_causes_self_inhibition():
    """Von Kiedrowski's real, counterintuitive 1986 finding: the same
    self-complementarity that makes TEMPLATE able to catalyze its own
    formation also lets two TEMPLATEs hybridize with each other, taking
    both out of circulation as inert DUPLEX. This checks that turning on
    duplex formation collapses the catalytically-active (single-stranded)
    TEMPLATE pool, while the total amount of material that got converted
    out of FRAGMENT (active TEMPLATE + 2*DUPLEX, i.e. "how much ligation
    chemistry actually happened") stays the same -- i.e. duplex formation
    doesn't reduce how much gets made, it just locks up what already
    formed."""

    def run(k_duplex):
        params = PolymerParams(k_oligomerize=0.0, k_background=0.01, k_template=2.0,
                                k_duplex=k_duplex, k_melt=0.0, t_max=100.0, max_events=300_000,
                                sample_every=1000, seed=13)
        init = {FRAGMENT.canonical_id(): 1000, TEMPLATE.canonical_id(): 2}
        return run_polymer(params, init)

    result_off = run(0.0)
    result_on = run(1.0)

    frag_left_off = result_off.counts.get(FRAGMENT.canonical_id(), 0)
    frag_left_on = result_on.counts.get(FRAGMENT.canonical_id(), 0)
    check("without duplex formation, fragment pool is fully converted by t_max", frag_left_off == 0)
    check("with duplex formation on, fragment pool is also fully converted (same conversion capacity)", frag_left_on == 0)

    active_off = result_off.counts.get(TEMPLATE.canonical_id(), 0)
    active_on = result_on.counts.get(TEMPLATE.canonical_id(), 0)
    duplex_off = result_off.counts.get(DUPLEX.canonical_id(), 0)
    duplex_on = result_on.counts.get(DUPLEX.canonical_id(), 0)
    equiv_off = active_off + 2 * duplex_off
    equiv_on = active_on + 2 * duplex_on
    print(f"  -> duplex off: active_template={active_off} duplex={duplex_off} (template-equivalent={equiv_off})")
    print(f"  -> duplex on:  active_template={active_on} duplex={duplex_on} (template-equivalent={equiv_on})")

    check("without duplex formation, essentially all converted material is catalytically active free template",
          active_off > 400)
    check("with duplex formation on, the catalytically active free template pool is largely sequestered away",
          active_on < active_off * 0.2)
    check("total template-equivalent material converted is the same either way -- "
          "duplex formation locks up product, it doesn't reduce total conversion",
          abs(equiv_on - equiv_off) <= 5)


def test_selection_reaction_structure():
    reactions_no_mut = build_selection_reactions(SelectionParams(
        variants=("A", "B"), k_template={"A": 1.0, "B": 1.0}, k_mutation=0.0,
    ))
    kinds_no_mut = {r.kind for r in reactions_no_mut}
    check("with mutation off, no mutation reactions are generated at all",
          "mutation" not in kinds_no_mut)
    check("5 reaction kinds per variant x 2 variants = 10 reactions with mutation off",
          len(reactions_no_mut) == 10)

    reactions_mut = build_selection_reactions(SelectionParams(
        variants=("A", "B"), k_template={"A": 1.0, "B": 1.0}, k_mutation=0.5,
    ))
    mutation_reactions = [r for r in reactions_mut if r.kind == "mutation"]
    check("with mutation on, one mutation reaction exists per ordered (v, w) pair, v != w",
          len(mutation_reactions) == 2)
    check("mutation reaction for A->B consumes 2 FRAGMENT-A + 1 TEMPLATE-A",
          sorted(next(r.reactant_ids for r in mutation_reactions if r.key.endswith("A->template-A+template-B")))
          == sorted((Oligomer(3, "A").canonical_id(), Oligomer(3, "A").canonical_id(), Oligomer(6, "A").canonical_id())))
    check("mutation reaction for A->B produces TEMPLATE-A (regenerated) + TEMPLATE-B (the mutant)",
          next(r.products for r in mutation_reactions if r.key.endswith("A->template-A+template-B"))
          == (Oligomer(6, "A"), Oligomer(6, "B")))

    reactions_3var = build_selection_reactions(SelectionParams(
        variants=("A", "B", "C"), k_template={"A": 1.0, "B": 1.0, "C": 1.0}, k_mutation=0.5,
    ))
    check("with 3 variants and mutation on, 6 ordered mutation pairs exist (3 * 2)",
          len([r for r in reactions_3var if r.kind == "mutation"]) == 6)


def test_selection_mass_conservation():
    """Same technique as polymer.py's conservation test, extended across a
    shared RIBONUCLEOTIDE pool and 2 competing variants: every reaction,
    including mutation (which reallocates units between variants but
    creates or destroys none), is a pure reorganization of
    ribonucleotide-equivalent units."""
    params = SelectionParams(
        variants=("A", "B"), k_oligomerize=0.01, k_background=0.01,
        k_template={"A": 1.0, "B": 2.0}, k_duplex=0.3, k_melt=0.05, k_mutation=0.1,
        t_max=200.0, max_events=300_000, sample_every=500, seed=5,
    )
    initial_ribo = 3000
    result = run_selection(params, {RIBONUCLEOTIDE.canonical_id(): initial_ribo,
                                     Oligomer(6, "A").canonical_id(): 2})
    initial_units = initial_ribo + 2 * 6

    for snap in result.history_counts:
        total = snap.get(RIBONUCLEOTIDE.canonical_id(), 0)
        for v in params.variants:
            total += snap.get(Oligomer(3, v).canonical_id(), 0) * 3
            total += snap.get(Oligomer(6, v).canonical_id(), 0) * 6
            total += snap.get(Duplex(v).canonical_id(), 0) * 12
        check(f"ribonucleotide-unit total conserved across variants at every sampled point "
              f"(got {total}, want {initial_units})", total == initial_units)
    print(f"  -> conserved at {len(result.history_counts)} sampled points, final state: {result.counts}")


def test_selection_mutation_is_required_for_unseeded_variant_to_appear():
    """With k_background=0, the ONLY way a variant that was never seeded can
    come into existence at all is a mutation event off an existing
    different variant's copying. This is the control this module's central
    experiment depends on."""

    def run(k_mutation, seed):
        params = SelectionParams(
            variants=("A", "B"), k_oligomerize=0.0005, k_background=0.0,
            k_template={"A": 1.0, "B": 10.0}, k_duplex=0.0, k_melt=0.0, k_mutation=k_mutation,
            t_max=400.0, max_events=600_000, sample_every=1000, seed=seed,
        )
        init = {RIBONUCLEOTIDE.canonical_id(): 20000, Oligomer(6, "A").canonical_id(): 2}
        return run_selection(params, init)

    for seed in (1, 2, 3):
        result = run(0.0, seed)
        check(f"with mutation off (seed={seed}), never-seeded variant B stays exactly 0",
              result.counts.get(Oligomer(6, "B").canonical_id(), 0) == 0)

    result_mut = run(0.05, 1)
    b_count = result_mut.counts.get(Oligomer(6, "B").canonical_id(), 0)
    print(f"  -> mutation on: variant B (never seeded) reached {b_count}")
    check("with mutation on, variant B appears despite never being seeded", b_count > 0)


def test_selection_fitter_variant_wins_direct_competition():
    """The cleanest isolation of fitness: both variants pre-seeded with
    identical fragment/template stock (no oligomerization/background, so
    no shared-resource or emergence-timing confound), differing only in
    k_template. Compared mid-transient (both are still actively growing,
    not yet saturated against their shared fixed fragment ceiling -- given
    infinite time both would fully consume their own stock regardless of
    rate, the same lesson polymer.py's own acceleration test needed)."""
    params = SelectionParams(
        variants=("A", "B"), k_oligomerize=0.0, k_background=0.0,
        k_template={"A": 0.0003, "B": 0.0009}, k_duplex=0.0, k_melt=0.0, k_mutation=0.0,
        t_max=0.25, max_events=300_000, sample_every=200, seed=1,
    )
    init = {
        Oligomer(3, "A").canonical_id(): 200, Oligomer(3, "B").canonical_id(): 200,
        Oligomer(6, "A").canonical_id(): 2, Oligomer(6, "B").canonical_id(): 2,
    }
    result = run_selection(params, init)
    a = result.counts.get(Oligomer(6, "A").canonical_id(), 0)
    b = result.counts.get(Oligomer(6, "B").canonical_id(), 0)
    print(f"  -> equal starting stock, 3x k_template: A={a} B={b}")
    check("the 3x-faster-templating variant pulls decisively ahead from identical starting conditions",
          b > a * 2)


def test_selection_fitter_mutant_overtakes_established_lineage():
    """The capstone experiment: variant B is never seeded (only exists if
    mutation creates it from A's copying), starts at a strict timing
    disadvantage, but is 10x fitter. Does it not only appear but actually
    overtake the already-established, pre-seeded variant A by the time the
    shared food supply runs out?"""
    params = SelectionParams(
        variants=("A", "B"), k_oligomerize=0.0005, k_background=0.0,
        k_template={"A": 1.0, "B": 10.0}, k_duplex=0.0, k_melt=0.0, k_mutation=0.05,
        t_max=400.0, max_events=600_000, sample_every=1000, seed=1,
    )
    init = {RIBONUCLEOTIDE.canonical_id(): 20000, Oligomer(6, "A").canonical_id(): 2}
    result = run_selection(params, init)
    a = result.counts.get(Oligomer(6, "A").canonical_id(), 0)
    b = result.counts.get(Oligomer(6, "B").canonical_id(), 0)
    print(f"  -> A (seeded from the start): {a}   B (arose only via mutation): {b}")
    check("B actually exists (mutation successfully introduced it)", b > 0)
    check("the fitter variant, despite starting from zero and arising later, overtakes "
          "the established lineage by the time the shared food supply is exhausted",
          b > a)


def test_hypercycle_reaction_structure():
    variants = ("A", "B", "C")
    reactions_closed = build_hypercycle_reactions(HypercycleParams(variants=variants, closed=True, k_cross=1.0))
    reactions_open = build_hypercycle_reactions(HypercycleParams(variants=variants, closed=False, k_cross=1.0))

    cross_closed = [r for r in reactions_closed if r.kind == "cross_catalyzed_ligation"]
    cross_open = [r for r in reactions_open if r.kind == "cross_catalyzed_ligation"]
    check("a closed 3-member cycle has 3 cross-catalysis edges (A->B, B->C, C->A)", len(cross_closed) == 3)
    check("an open 3-member chain has only 2 cross-catalysis edges (A->B, B->C, no C->A)", len(cross_open) == 2)
    check("closing the loop adds exactly one reaction (the wraparound edge) and nothing else",
          len(reactions_closed) == len(reactions_open) + 1)

    wraparound_key = "cross_catalyzed_ligation:2fragment-A+template-C->template-C+template-A"
    check("the wraparound edge (C catalyzes A) exists when closed", any(r.key == wraparound_key for r in reactions_closed))
    check("the wraparound edge does not exist when open", not any(r.key == wraparound_key for r in reactions_open))

    a_to_b = next(r for r in reactions_closed if r.key == "cross_catalyzed_ligation:2fragment-B+template-A->template-A+template-B")
    check("A catalyzing B's replication consumes 2 FRAGMENT-B + 1 TEMPLATE-A",
          sorted(a_to_b.reactant_ids) == sorted((Oligomer(3, "B").canonical_id(), Oligomer(3, "B").canonical_id(), Oligomer(6, "A").canonical_id())))
    check("A catalyzing B's replication regenerates TEMPLATE-A (catalyst) and produces TEMPLATE-B",
          a_to_b.products == (Oligomer(6, "A"), Oligomer(6, "B")))


def test_hypercycle_mass_conservation():
    variants = ("A", "B", "C")
    params = HypercycleParams(
        variants=variants, closed=True, k_oligomerize=0.01, k_background=0.01,
        k_self={"A": 0.5}, k_cross=0.3, k_duplex=0.2, k_melt=0.05,
        t_max=200.0, max_events=300_000, sample_every=500, seed=5,
    )
    initial_ribo = 3000
    init = {RIBONUCLEOTIDE.canonical_id(): initial_ribo}
    for v in variants:
        init[Oligomer(6, v).canonical_id()] = 2
    result = run_hypercycle(params, init)
    initial_units = initial_ribo + len(variants) * 2 * 6

    for snap in result.history_counts:
        total = snap.get(RIBONUCLEOTIDE.canonical_id(), 0)
        for v in variants:
            total += snap.get(Oligomer(3, v).canonical_id(), 0) * 3
            total += snap.get(Oligomer(6, v).canonical_id(), 0) * 6
            total += snap.get(Duplex(v).canonical_id(), 0) * 12
        check(f"ribonucleotide-unit total conserved across the cycle at every sampled point "
              f"(got {total}, want {initial_units})", total == initial_units)
    print(f"  -> conserved at {len(result.history_counts)} sampled points, final state: {result.counts}")


def test_hypercycle_closing_the_loop_lets_an_unhelped_member_grow():
    """The central experiment: with every member unable to replicate alone
    (k_self all 0) and background ligation off (so a TEMPLATE's count can
    ONLY change by being produced as a reaction product), variant A has no
    possible way to grow in an open chain (A->B->C, nothing points to A) --
    its count is mathematically pinned at its seed value. Closing the loop
    (C->A) is the only structural change; does it actually let A grow?"""
    variants = ("A", "B", "C")

    def run(closed, seed):
        params = HypercycleParams(
            variants=variants, closed=closed, k_oligomerize=0.001, k_background=0.0,
            k_self={}, k_cross=0.001, k_duplex=0.0, k_melt=0.0,
            t_max=500.0, max_events=400_000, sample_every=1000, seed=seed,
        )
        init = {RIBONUCLEOTIDE.canonical_id(): 20000}
        for v in variants:
            init[Oligomer(6, v).canonical_id()] = 5
        return run_hypercycle(params, init)

    for seed in (1, 2, 3):
        result_open = run(False, seed)
        a_open = result_open.counts.get(Oligomer(6, "A").canonical_id(), 0)
        check(f"open chain (seed={seed}): variant A, which nothing catalyzes, stays exactly at its seed value",
              a_open == 5)

    result_closed = run(True, 1)
    a_closed = result_closed.counts.get(Oligomer(6, "A").canonical_id(), 0)
    b_closed = result_closed.counts.get(Oligomer(6, "B").canonical_id(), 0)
    c_closed = result_closed.counts.get(Oligomer(6, "C").canonical_id(), 0)
    print(f"  -> open chain: A stays at seed value (5)")
    print(f"  -> closed cycle: A={a_closed} B={b_closed} C={c_closed}")
    check("closing the loop lets A -- which the open chain could never help -- grow far past its seed value",
          a_closed > 100)


def test_hypercycle_parasite_reaction_structure():
    variants = ("A", "B", "C")
    reactions_no_parasite = build_hypercycle_reactions(HypercycleParams(variants=variants, closed=True, k_cross=1.0))
    reactions_with_parasite = build_hypercycle_reactions(HypercycleParams(
        variants=variants, closed=True, k_cross=1.0, parasite="P", parasite_catalyst="A", k_parasite=0.5,
    ))
    check("with no parasite configured, no parasite-related reactions or species show up at all",
          not any("-P" in r.key for r in reactions_no_parasite))
    check("setting parasite adds exactly 5 new reactions (oligomerization, background_ligation, "
          "duplex_formation, duplex_melting, parasite_catalyzed_ligation)",
          len(reactions_with_parasite) == len(reactions_no_parasite) + 5)

    parasite_lig = next(r for r in reactions_with_parasite if r.kind == "parasite_catalyzed_ligation")
    check("parasite_catalyzed_ligation consumes 2 FRAGMENT-P + 1 TEMPLATE-A (the catalyst)",
          sorted(parasite_lig.reactant_ids) == sorted((Oligomer(3, "P").canonical_id(), Oligomer(3, "P").canonical_id(), Oligomer(6, "A").canonical_id())))
    check("parasite_catalyzed_ligation regenerates TEMPLATE-A and produces TEMPLATE-P",
          parasite_lig.products == (Oligomer(6, "A"), Oligomer(6, "P")))
    check("the parasite has NO outgoing cross_catalyzed_ligation reaction of its own (it can only be a product)",
          not any(r.kind == "cross_catalyzed_ligation" and r.reactant_ids[-1] == Oligomer(6, "P").canonical_id()
                  for r in reactions_with_parasite))
    check("the parasite gets NO self_templated_ligation reaction either (it cannot replicate itself)",
          not any(r.kind == "self_templated_ligation" and "P" in r.key for r in reactions_with_parasite))


def test_hypercycle_parasite_mass_conservation():
    variants = ("A", "B", "C")
    params = HypercycleParams(
        variants=variants, closed=True, k_oligomerize=0.01, k_background=0.01,
        k_self={"A": 0.3}, k_cross=0.3, k_duplex=0.2, k_melt=0.05,
        parasite="P", parasite_catalyst="A", k_parasite=0.3,
        t_max=200.0, max_events=300_000, sample_every=500, seed=5,
    )
    initial_ribo = 3000
    init = {RIBONUCLEOTIDE.canonical_id(): initial_ribo}
    for v in variants:
        init[Oligomer(6, v).canonical_id()] = 2
    result = run_hypercycle(params, init)
    initial_units = initial_ribo + len(variants) * 2 * 6

    all_names = list(variants) + ["P"]
    for snap in result.history_counts:
        total = snap.get(RIBONUCLEOTIDE.canonical_id(), 0)
        for v in all_names:
            total += snap.get(Oligomer(3, v).canonical_id(), 0) * 3
            total += snap.get(Oligomer(6, v).canonical_id(), 0) * 6
            total += snap.get(Duplex(v).canonical_id(), 0) * 12
        check(f"ribonucleotide-unit total conserved with parasite active at every sampled point "
              f"(got {total}, want {initial_units})", total == initial_units)
    print(f"  -> conserved at {len(result.history_counts)} sampled points, final state: {result.counts}")


def test_hypercycle_parasite_costs_the_cycle_regardless_of_virulence():
    """The actual (not the naively assumed) finding: the parasite's mere
    existence costs the legitimate cycle a share of the shared food supply
    -- but how EFFICIENT it is at converting that share (k_parasite) barely
    matters, because at full resource exhaustion final extent stops
    depending on rate (the same principle polymer.py's and selection.py's
    own tests ran into, in a new context)."""
    variants = ("A", "B", "C")

    def legit_total(parasite, k_parasite, seed):
        params = HypercycleParams(
            variants=variants, closed=True, k_oligomerize=0.001, k_background=0.0,
            k_self={}, k_cross=0.001, k_duplex=0.0, k_melt=0.0,
            parasite=parasite, parasite_catalyst="A" if parasite else None, k_parasite=k_parasite,
            t_max=500.0, max_events=400_000, sample_every=1000, seed=seed,
        )
        init = {RIBONUCLEOTIDE.canonical_id(): 20000}
        for v in variants:
            init[Oligomer(6, v).canonical_id()] = 5
        result = run_hypercycle(params, init)
        return sum(result.counts.get(Oligomer(6, v).canonical_id(), 0) for v in variants)

    no_parasite = legit_total(None, 0.0, seed=1)
    with_parasite = [legit_total("P", k, seed=1) for k in (0.001, 0.02, 0.5, 2.0)]
    print(f"  -> no parasite: legit total={no_parasite}")
    print(f"  -> with parasite, k_parasite in (0.001, 0.02, 0.5, 2.0): legit totals={with_parasite}")

    check("a parasite's mere presence costs the legitimate cycle a substantial share of output",
          all(v < no_parasite * 0.85 for v in with_parasite))
    check("but that cost barely depends on how virulent (efficient) the parasite is -- "
          "the spread across a 2000x k_parasite range is small",
          max(with_parasite) - min(with_parasite) < no_parasite * 0.05)


def test_hypercycle_parasite_wastes_resources_a_legitimate_member_would_have_used():
    """The sharper, fairer comparison: is a parasite actually WORSE than an
    equally hungry legitimate competitor? Total output (legitimate members
    + parasite) should be about the same whether the 4th claimant on the
    shared food supply is a true parasite or a fully legitimate 4th member
    closing a real 4-membered loop -- a parasite doesn't claim MORE than a
    cooperator would. What it does is waste what it claims: legitimate-only
    output should be substantially lower with a parasite present than with
    a legitimate 4th member in its place."""

    def run_with_parasite(seed):
        params = HypercycleParams(
            variants=("A", "B", "C"), closed=True, k_oligomerize=0.001, k_background=0.0,
            k_self={}, k_cross=0.001, k_duplex=0.0, k_melt=0.0,
            parasite="P", parasite_catalyst="A", k_parasite=0.5,
            t_max=500.0, max_events=400_000, sample_every=1000, seed=seed,
        )
        init = {RIBONUCLEOTIDE.canonical_id(): 20000}
        for v in ("A", "B", "C"):
            init[Oligomer(6, v).canonical_id()] = 5
        return run_hypercycle(params, init)

    def run_all_legit(seed):
        params = HypercycleParams(
            variants=("A", "B", "C", "D"), closed=True, k_oligomerize=0.001, k_background=0.0,
            k_self={}, k_cross=0.001, k_duplex=0.0, k_melt=0.0,
            t_max=500.0, max_events=400_000, sample_every=1000, seed=seed,
        )
        init = {RIBONUCLEOTIDE.canonical_id(): 20000}
        for v in ("A", "B", "C", "D"):
            init[Oligomer(6, v).canonical_id()] = 5
        return run_hypercycle(params, init)

    for seed in (1, 2, 3):
        r_parasite = run_with_parasite(seed)
        r_legit = run_all_legit(seed)
        legit3 = sum(r_parasite.counts.get(Oligomer(6, v).canonical_id(), 0) for v in ("A", "B", "C"))
        parasite_count = r_parasite.counts.get(Oligomer(6, "P").canonical_id(), 0)
        legit4 = sum(r_legit.counts.get(Oligomer(6, v).canonical_id(), 0) for v in ("A", "B", "C", "D"))
        print(f"  -> seed={seed}: 3-legit+parasite total={legit3 + parasite_count} (legit-only={legit3}); "
              f"4-legit total={legit4}")
        check(f"seed={seed}: total claimed resource is about the same either way "
              f"(parasite doesn't out-compete a legitimate member for raw share)",
              abs((legit3 + parasite_count) - legit4) < legit4 * 0.05)
        check(f"seed={seed}: but legitimate-only output is substantially lower with a parasite "
              f"present than with a legitimate 4th member in its place",
              legit3 < legit4 * 0.85)


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
    test_electricity_indirectly_enables_o2_via_hydrogen_sink()
    test_chain_length_stats_basic()
    test_higher_uv_yields_longer_chains_at_fixed_time()
    test_formose_reaction_structure()
    test_formose_carbon_conservation_and_chain_growth()
    test_formose_ribose_needs_stabilization()
    test_formose_stereoisomer_tracking_structure()
    test_formose_ribose_selectivity_is_a_real_correction_factor()
    test_nucleotide_reaction_structure()
    test_nucleotide_forms_and_photoselection_protects_scarce_cyanoacetylene()
    test_reaction_propensity_generalizes_to_n_body()
    test_polymer_reaction_structure()
    test_polymer_mass_conservation()
    test_polymer_templating_accelerates_strand_formation()
    test_polymer_duplex_formation_causes_self_inhibition()
    test_selection_reaction_structure()
    test_selection_mass_conservation()
    test_selection_mutation_is_required_for_unseeded_variant_to_appear()
    test_selection_fitter_variant_wins_direct_competition()
    test_selection_fitter_mutant_overtakes_established_lineage()
    test_hypercycle_reaction_structure()
    test_hypercycle_mass_conservation()
    test_hypercycle_closing_the_loop_lets_an_unhelped_member_grow()
    test_hypercycle_parasite_reaction_structure()
    test_hypercycle_parasite_mass_conservation()
    test_hypercycle_parasite_costs_the_cycle_regardless_of_virulence()
    test_hypercycle_parasite_wastes_resources_a_legitimate_member_would_have_used()
    print("\nAll sanity checks passed.")
