"""
Streamlit app: stochastic UV-photolysis hydrocarbon chemistry, starting from
methane, to explore whether growing reaction-network complexity produces
candidate autocatalytic (self-amplifying) loops -- a prerequisite substrate
for anything Darwinian.

Run with: streamlit run app.py
"""
import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from engine.analysis import chain_length_stats
from engine.autocatalysis import build_flow_graph, find_candidate_cycles
from engine.formose import PENTOSE_DIASTEREOMERS, RIBOSE_FRACTION_ESTIMATE, FormoseParams, run_formose
from engine.molecule import ATMOSPHERIC_GASES, Molecule, SEED_MOLECULES, ethane, methane
from engine.simulator import SimParams, Simulator

st.set_page_config(page_title="Abiogenesis: hydrocarbon photochemistry", layout="wide")

st.title("Methane photochemistry -> reaction-network complexity")
st.caption(
    "CH4 + UV -> CH3• + H•, then radical combination and H-abstraction let the "
    "network grow on its own. Add N2/O2/CO2/Ar/H2O to test the central question: does "
    "the C2H6 that forms actually *accumulate*, or does O2 (direct, or self-generated from "
    "water photolysis) outcompete radical self-combination and divert it into oxidized "
    "products instead? Add electricity (discharge) to see UV's blind spot: N2's triple bond "
    "is a bond UV in this model can never break, but a spark can -- the Miller-Urey question. "
    "This is a toy, non-calibrated model -- see README for what it does and doesn't capture."
)

ATMOSPHERE_PRESETS = {
    "Pure methane (baseline, no atmosphere)": dict(n2=0, o2=0, co2=0, ar=0, h2o=0),
    "Reducing, early-Earth-like (trace O2)": dict(n2=200, o2=5, co2=20, ar=5, h2o=0),
    "Modern oxidizing atmosphere (O2-rich)": dict(n2=780, o2=210, co2=4, ar=10, h2o=0),
    "Wet planet, no free O2 (water photolysis only)": dict(n2=0, o2=0, co2=0, ar=0, h2o=300),
    "Custom": dict(n2=0, o2=0, co2=0, ar=0, h2o=0),
}

with st.sidebar:
    st.header("Setup")
    seed_name = st.selectbox("Starting molecule", list(SEED_MOLECULES.keys()), index=0)
    seed_count = st.slider("Starting molecule count", 10, 2000, 300, step=10)

    st.header("Atmosphere")
    preset_name = st.selectbox("Preset", list(ATMOSPHERE_PRESETS.keys()), index=0)
    defaults = ATMOSPHERE_PRESETS[preset_name]
    n2_count = st.slider("N2 count", 0, 2000, defaults["n2"], step=10, key=f"n2_{preset_name}")
    o2_count = st.slider("O2 count", 0, 2000, defaults["o2"], step=10, key=f"o2_{preset_name}")
    co2_count = st.slider("CO2 count", 0, 500, defaults["co2"], step=5, key=f"co2_{preset_name}")
    ar_count = st.slider("Ar count", 0, 500, defaults["ar"], step=5, key=f"ar_{preset_name}")
    h2o_count = st.slider("H2O count (a 'wet planet')", 0, 2000, defaults["h2o"], step=10, key=f"h2o_{preset_name}")
    st.caption(
        "Ar never reacts in this model (it's a noble gas). O2 and CO2 photolyze; O2 "
        "additionally *scavenges* hydrocarbon radicals directly, and separately builds up "
        "O3 (ozone) with UV-generated O atoms -- O3 is itself highly reactive and scavenges "
        "radicals too. H2O also photolyzes (H2O -> H* + OH*); two OH* can disproportionate "
        "into H2O + O*, which is how a wet planet can bootstrap its *own* O2/O3 from water "
        "alone, no free O2 required. N2 is UV-inert (its triple bond needs far harder UV "
        "than this sim implies) *unless* electricity is turned on below. See the callouts "
        "below the run for how this plays out."
    )

    st.header("Hydrogen escape")
    k_escape = st.slider(
        "H2 / H* escape rate to space (k_escape)", 0.0, 20.0, 0.0, step=0.1,
        help="Off (0) by default: a closed system just cycles oxygen back into water and "
             "never builds up free O2/O3, no matter how much UV or how much time. Turning "
             "escape on gives the system a one-way loss of hydrogen -- the real mechanism "
             "proposed for abiotic O2/O3 buildup on a wet planet with no biology at all.",
    )
    st.caption(
        "Applies equally to H2 and H* (both escape in this simplified model). In testing, "
        "even a nonzero escape rate only produces a modest, self-limited O2/O3 steady state "
        "(ozone's own photolysis balances its formation) unless escape and/or water-photolysis "
        "rates are pushed well above realistic values -- see README."
    )

    st.header("Electricity (discharge)")
    k_discharge_n2 = st.slider(
        "N2 + spark -> 2 N* (k_discharge_n2)", 0.0, 20.0, 0.0, step=0.1,
        help="Off (0) by default. UV in this model, at any intensity, never touches N2 -- "
             "its triple bond (~9.8 eV) is beyond what's driving C-H/O2 photolysis here. This "
             "is a qualitatively different energy source (electric spark, as in Miller-Urey) "
             "that can reach it. Turn it on to let nitrogen enter the hydrocarbon chemistry.",
    )
    st.caption(
        "Once N* exists it plugs into the same combination/abstraction machinery as "
        "everything else for free: N* + N* -> N2 (recombination), N* + H* -> NH*, and "
        "N* (or NH*) + a hydrocarbon radical -> a closed-shell amine (R-NH2) -- nitrogen "
        "incorporated into organic chemistry, the way Miller-Urey's spark got nitrogen into "
        "amino acid precursors. No N-O cross chemistry (real NOx) is modeled."
    )

    st.header("Rate constants (relative units)")
    k_photo = st.slider("UV photolysis, hydrocarbons (k_photo)", 0.001, 1.0, 0.05, step=0.001, format="%.3f")
    k_comb = st.slider("Radical combination (k_comb)", 0.1, 50.0, 5.0, step=0.1)
    k_abstr = st.slider("H-abstraction (k_abstr)", 0.01, 10.0, 0.5, step=0.01)
    k_photo_o2 = st.slider("UV photolysis, O2 (k_photo_o2)", 0.0, 1.0, 0.03, step=0.001, format="%.3f")
    k_photo_co2 = st.slider("UV photolysis, CO2 (k_photo_co2)", 0.0, 1.0, 0.005, step=0.001, format="%.3f")
    k_o2_scavenge = st.slider(
        "Radical + O2 scavenging (k_o2_scavenge)", 0.0, 50.0, 5.0, step=0.1,
        help="Same order of magnitude as k_comb by default, so which pathway wins is "
             "driven by relative O2 vs. radical *concentration*, not a thumb on the scale.",
    )
    k_o3_formation = st.slider("O* + O2 -> O3 formation (k_o3_formation)", 0.0, 50.0, 5.0, step=0.1)
    k_photo_o3 = st.slider("UV photolysis, O3 (k_photo_o3)", 0.0, 1.0, 0.08, step=0.001, format="%.3f")
    k_o3_scavenge = st.slider(
        "Radical + O3 scavenging (k_o3_scavenge)", 0.0, 50.0, 8.0, step=0.1,
        help="Ozone is set to react with radicals a bit faster than O2 by default, matching "
             "how reactive O3 actually is toward radical species in real atmospheric chemistry.",
    )
    k_photo_h2o = st.slider(
        "UV photolysis, H2O (k_photo_h2o)", 0.0, 1.0, 0.004, step=0.001, format="%.3f",
        help="Water is a weak UV absorber -- needs harder UV than O2/CO2 in reality, so this "
             "defaults low. Raise it to see the wet-planet mechanism proceed faster.",
    )
    k_oh_disprop = st.slider(
        "OH* + OH* -> H2O + O* (k_oh_disprop)", 0.0, 20.0, 2.0, step=0.1,
        help="The secondary step that lets water photolysis bootstrap free O atoms (and "
             "hence O2/O3) without any O2 present to begin with.",
    )

    st.header("Complexity & run length")
    max_carbon = st.slider("Max carbons per molecule (complexity ceiling)", 2, 10, 6)
    t_max = st.slider("Simulated time (t_max)", 10.0, 2000.0, 400.0, step=10.0)
    max_events = st.select_slider(
        "Max reaction events", options=[1000, 2000, 5000, 10000, 20000, 50000], value=10000
    )
    seed = st.number_input("Random seed", value=42, step=1)

    run_clicked = st.button("Run simulation", type="primary", use_container_width=True)

if run_clicked:
    params = SimParams(
        k_photo=k_photo, k_comb=k_comb, k_abstr=k_abstr,
        k_photo_o2=k_photo_o2, k_photo_co2=k_photo_co2, k_o2_scavenge=k_o2_scavenge,
        k_o3_formation=k_o3_formation, k_photo_o3=k_photo_o3, k_o3_scavenge=k_o3_scavenge,
        k_photo_h2o=k_photo_h2o, k_oh_disprop=k_oh_disprop,
        k_escape_h2=k_escape, k_escape_h=k_escape,
        k_discharge_n2=k_discharge_n2,
        max_carbon=max_carbon, t_max=t_max, max_events=max_events,
        sample_every=max(1, max_events // 500), seed=int(seed),
    )
    sim = Simulator(params)
    sim.seed_species(SEED_MOLECULES[seed_name](), seed_count)
    for gas_name, count in [("N2", n2_count), ("O2", o2_count), ("CO2", co2_count),
                             ("Ar", ar_count), ("H2O", h2o_count)]:
        if count > 0:
            sim.seed_species(ATMOSPHERIC_GASES[gas_name](), count)
    with st.spinner("Running Gillespie simulation..."):
        result = sim.run()
    st.session_state["result"] = result
    st.session_state["params"] = params

if "result" not in st.session_state:
    st.info("Set parameters in the sidebar and click **Run simulation** to begin.")
    st.stop()

result = st.session_state["result"]
params = st.session_state["params"]

c1, c2, c3, c4 = st.columns(4)
c1.metric("Species discovered", len(result.species))
c2.metric("Reactions in network", len(result.reactions))
c3.metric("Events fired", len(result.event_log))
c4.metric("Stopped because", result.stopped_reason)

# -- Does C2H6 accumulate? ---------------------------------------------------
comb_fires = sum(
    n for k, n in result.reaction_fire_counts.items()
    if k.startswith("combination:") and all(
        isinstance(result.species.get(rid), Molecule)
        for rid in k.split(":", 1)[1].split("->")[0].split("+")
    )
)
o2_scav_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("o2_scavenging:"))
o3_scav_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("ozone_scavenging:"))
scav_fires = o2_scav_fires + o3_scav_fires
o3_now = result.counts.get("O3", 0)
c2h6_id = ethane().canonical_id()
c2h6_now = result.counts.get(c2h6_id, 0)
c2h6_peak = max((snap.get(c2h6_id, 0) for snap in result.history_counts), default=0)

ec1, ec2, ec3, ec4 = st.columns(4)
ec1.metric("Hydrocarbon self-combinations fired", comb_fires)
ec2.metric("O2 + O3 scavenging events fired", scav_fires, help=f"O2: {o2_scav_fires}, O3: {o3_scav_fires}")
ec3.metric("O3 (ozone) now", o3_now)
ec4.metric("C2H6 now / peak", f"{c2h6_now} / {c2h6_peak}")

o2_now = result.counts.get("O2", 0)
oxidizer_present = o2_count > 0 or o2_now > 0 or o3_now > 0

if h2o_count > 0 and o2_count == 0 and k_escape == 0.0:
    st.success(
        f"**Wet planet, closed system**: {o2_now} O2 / {o3_now} O3 despite water photolysis "
        f"being active -- with hydrogen escape off, every O atom pulled off H2O finds its way "
        f"back into H2O via this same radical chemistry, so no free oxidant can persist. "
        f"Self-combination proceeds essentially unopposed ({comb_fires} events, {scav_fires} scavenging). "
        f"Turn on hydrogen escape (sidebar) to let this system actually build up an oxidizing atmosphere."
    )
elif not oxidizer_present:
    st.info(
        "No O2 (direct or water-generated) in this run -- self-combination has no "
        "competition, so any C2H6 formed just accumulates."
    )
elif scav_fires > comb_fires * 3:
    st.warning(
        f"**O2/O3 win this competition** ({scav_fires} scavenging events -- {o2_scav_fires} direct O2, "
        f"{o3_scav_fires} via ozone -- vs {comb_fires} self-combinations). Hydrocarbon radicals are being "
        "captured into peroxy/alkoxy radicals faster than they can dimerize into C2H6 -- consistent with why "
        "an oxidizing atmosphere is hostile to hydrocarbon accumulation, and why prebiotic chemistry is "
        "generally thought to need a *reducing* (low-O2) atmosphere. Note O2 exposed to UV doesn't just sit "
        "there either: it builds up O3, which is itself a fast, independent radical sink (see O3 count above)."
    )
elif comb_fires > scav_fires * 3:
    st.success(
        f"**Self-combination wins** ({comb_fires} vs {scav_fires} O2/O3 scavenging events) -- O2 is too scarce "
        "here to intercept much of the radical pool, so C2H6 forms largely unopposed."
    )
else:
    st.info(f"The two pathways are roughly balanced ({comb_fires} self-combination vs {scav_fires} O2/O3 scavenging).")

# -- Does electricity reach nitrogen where UV can't? -------------------------
if n2_count > 0:
    n2_now = result.counts.get("N2", 0)
    n2_consumed = n2_count - n2_now
    discharge_fires = sum(n for k, n in result.reaction_fire_counts.items() if k.startswith("discharge:"))
    amine_count = sum(
        n for sid, n in result.counts.items() if sid.startswith("RNH2:")
    )
    dc1, dc2, dc3 = st.columns(3)
    dc1.metric("N2 remaining", f"{n2_now} / {n2_count}")
    dc2.metric("Discharge events fired", discharge_fires)
    dc3.metric("Amine (R-NH2) molecules formed", amine_count)
    if k_discharge_n2 == 0.0:
        st.info(
            f"**Electricity is off**: N2 stayed at exactly {n2_now}/{n2_count}, completely untouched, no "
            "matter how much UV or simulated time ran. This is the point -- UV in this model genuinely "
            "cannot break N2's triple bond. Turn on the discharge rate (sidebar) to change that."
        )
    elif n2_consumed > 0:
        st.success(
            f"**Electricity cracked N2**: {n2_consumed} of {n2_count} N2 molecules dissociated, feeding "
            f"{amine_count} amine molecules (nitrogen incorporated into hydrocarbon radicals) -- the Miller-"
            "Urey mechanism in miniature. UV alone, at any intensity, cannot do this in this model."
        )
    else:
        st.warning("Discharge is on but no N2 has dissociated yet -- try a higher rate or longer run.")

    # -- Three-way competition: self-combination vs O2/O3 scavenging vs amination --
    amine_fires = sum(
        n for k, n in result.reaction_fire_counts.items()
        if k.startswith("combination:") and "RNH2:" in k.split("->", 1)[1]
    )
    total_carbon = sum(sp.n_carbon * result.counts.get(sid, 0) for sid, sp in result.species.items())
    hydrocarbon_carbon = sum(
        sp.n_carbon * result.counts.get(sid, 0) for sid, sp in result.species.items() if isinstance(sp, Molecule)
    )
    amine_carbon = sum(
        sp.n_carbon * result.counts.get(sid, 0) for sid, sp in result.species.items() if sid.startswith("RNH2:")
    )
    oxidized_carbon = sum(
        sp.n_carbon * result.counts.get(sid, 0) for sid, sp in result.species.items()
        if sid.startswith(("ROO:", "RO:", "ROOH:", "ROH:"))
    )

    st.subheader("Three-way competition: where do hydrocarbon radicals actually go?")
    st.caption(
        "Every hydrocarbon radical has up to three fates: combine with another hydrocarbon "
        "radical (grows a chain), get scavenged by O2/O3 (oxidized product), or combine with "
        "N*/NH* (amine). Electricity doesn't touch water chemistry directly, but nitrogen "
        "fixation consumes H* (locking it into NH3/amines), which starves the H* + OH* -> H2O "
        "recombination that normally keeps a closed system oxidant-free -- so electricity can "
        "indirectly produce a little O2/O3 even with hydrogen escape off, as a side effect."
    )
    fc1, fc2, fc3 = st.columns(3)
    fc1.metric("Self-combination events", comb_fires)
    fc2.metric("O2/O3 scavenging events", scav_fires)
    fc3.metric("Amine-forming events", amine_fires)

    if total_carbon > 0:
        fig = go.Figure(go.Bar(
            x=["Plain hydrocarbons", "Amines (N)", "Oxidized (O)"],
            y=[hydrocarbon_carbon, amine_carbon, oxidized_carbon],
            marker_color=["#2a6fdb", "#2ca02c", "#d62728"],
        ))
        fig.update_layout(yaxis_title="carbon atoms currently held", height=350)
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            f"{hydrocarbon_carbon}/{total_carbon} ({hydrocarbon_carbon/total_carbon:.0%}) of all carbon is in "
            f"plain hydrocarbons, {amine_carbon}/{total_carbon} ({amine_carbon/total_carbon:.0%}) in amines, "
            f"{oxidized_carbon}/{total_carbon} ({oxidized_carbon/total_carbon:.0%}) in O-containing products."
        )

    if k_discharge_n2 > 0 and h2o_count > 0:
        disprop_fires = sum(n for k, n in result.reaction_fire_counts.items()
                             if k.startswith("oh_disproportionation:"))
        h2_now = result.counts.get("H2", 0)
        h_now = result.counts.get("H", 0)
        if amine_fires > max(comb_fires, scav_fires) * 1.5:
            st.warning(
                f"**Nitrogen fixation dominates** ({amine_fires} amine-forming events vs {comb_fires} "
                f"self-combination and {scav_fires} O2/O3 scavenging) -- this, not the {o3_now} O3 that "
                "formed as a side effect, is what's actually suppressing other molecules. H2/H* remaining: "
                f"{h2_now}/{h_now} -- nitrogen fixation has consumed the hydrogen that would otherwise "
                f"reform H2O, which is why {disprop_fires} OH-disproportionation events fired and produced "
                "free O atoms even with hydrogen escape off."
            )
        elif disprop_fires > 0:
            st.info(
                f"{disprop_fires} OH-disproportionation events fired ({o3_now} O3 now) -- nitrogen fixation "
                "is providing a small hydrogen-sink effect here, but hasn't come to dominate the radical "
                "pool the way it can at higher discharge rates or longer run times."
            )

tab_pop, tab_chain, tab_net, tab_cycles, tab_species = st.tabs(
    ["Population dynamics", "Chain-length distribution", "Reaction network",
     "Autocatalytic cycles", "Species inventory"]
)

# -- Population dynamics ---------------------------------------------------
with tab_pop:
    peak = {}
    for snap in result.history_counts:
        for sid, n in snap.items():
            peak[sid] = max(peak.get(sid, 0), n)
    top_ids = sorted(peak, key=peak.get, reverse=True)[:12]

    rows = []
    for t, snap in zip(result.history_t, result.history_counts):
        for sid in top_ids:
            rows.append({"t": t, "species": result.species[sid].formula(), "count": snap.get(sid, 0)})
    df = pd.DataFrame(rows)

    fig = go.Figure()
    for formula, sub in df.groupby("species"):
        fig.add_trace(go.Scatter(x=sub["t"], y=sub["count"], mode="lines", name=formula))
    fig.update_layout(
        yaxis_title="molecule count", xaxis_title="time (arbitrary units)",
        yaxis_type="log", height=500, legend_title="species (top 12 by peak count)",
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Log-scale y-axis. Only the 12 species that reached the highest peak count are shown.")

# -- Chain-length distribution ------------------------------------------------
with tab_chain:
    stats = chain_length_stats(result)
    cc1, cc2, cc3 = st.columns(3)
    cc1.metric("Mean carbon number (hydrocarbons)", f"{stats.mean_carbon:.2f}")
    cc2.metric("Longest chain reached", f"C{stats.max_carbon_present}")
    cc3.metric("Carbon in C3+ chains", f"{stats.c3plus_carbon_fraction:.1%}")

    if stats.counts_by_carbon:
        carbons = sorted(stats.counts_by_carbon)
        fig = go.Figure(go.Bar(
            x=[f"C{c}" for c in carbons], y=[stats.counts_by_carbon[c] for c in carbons],
        ))
        fig.update_layout(
            yaxis_title="molecule count (log scale)", xaxis_title="carbon chain length",
            yaxis_type="log", height=420,
        )
        st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "Total molecule count grouped by carbon number, across all hydrocarbon species "
        "(isomers pooled together) at the end of the run. Longer chains only appear here if "
        "combination events had enough radical concentration to outrun abstraction's "
        "radical-hopping -- see the Parameter sweep section below to see how UV level and "
        "starting concentration each push on this."
    )

# -- Reaction network -------------------------------------------------------
with tab_net:
    g = build_flow_graph(result.reactions, result.reaction_fire_counts)
    st.caption(
        f"{g.number_of_nodes()} species and {g.number_of_edges()} realized reaction-flow edges "
        "(only reactions that actually fired at least once)."
    )
    if g.number_of_nodes() == 0:
        st.warning("No reactions fired -- try raising rate constants or t_max/max_events.")
    else:
        pos = nx.spring_layout(g, seed=1, k=1.5 / max(1, g.number_of_nodes() ** 0.5))
        edge_x, edge_y = [], []
        for u, v in g.edges():
            edge_x += [pos[u][0], pos[v][0], None]
            edge_y += [pos[u][1], pos[v][1], None]
        node_x = [pos[n][0] for n in g.nodes()]
        node_y = [pos[n][1] for n in g.nodes()]
        node_labels = [result.species[n].formula() for n in g.nodes()]
        node_size = [8 + 2 * min(result.counts.get(n, 0), 20) for n in g.nodes()]

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                                  line=dict(width=0.5, color="#888"), hoverinfo="none"))
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y, mode="markers+text", text=node_labels, textposition="top center",
            marker=dict(size=node_size, color="#2a6fdb"),
        ))
        fig.update_layout(showlegend=False, height=650,
                           xaxis=dict(visible=False), yaxis=dict(visible=False))
        st.plotly_chart(fig, use_container_width=True)

# -- Autocatalytic cycles ----------------------------------------------------
with tab_cycles:
    st.caption(
        "Directed cycles in the realized reaction-flow graph, weighted by how often the "
        "limiting (bottleneck) reaction in the cycle actually fired. A short, high-flux "
        "cycle is a candidate chain-propagation / self-amplifying loop -- the closest thing "
        "in this hydrocarbon-only model to a chemical 'reproducing' something. See README "
        "for why this is a diagnostic, not proof of life."
    )
    cycles = find_candidate_cycles(result.reactions, result.reaction_fire_counts, result.species)
    if not cycles:
        st.write("No candidate cycles found in this run.")
    else:
        rows = [{
            "cycle": " -> ".join(c.formulas) + f" -> {c.formulas[0]}",
            "length": c.length,
            "bottleneck flux (fires)": c.bottleneck_flux,
        } for c in cycles]
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# -- Species inventory --------------------------------------------------------
with tab_species:
    rows = [{
        "formula": sp.formula(),
        "count": result.counts.get(sid, 0),
        "carbons": sp.n_carbon,
        "radical": sp.is_radical,
    } for sid, sp in result.species.items()]
    df = pd.DataFrame(rows).sort_values("count", ascending=False).reset_index(drop=True)
    st.dataframe(df, use_container_width=True, hide_index=True)

# -- Parameter sweep: concentration x UV level --------------------------------
st.divider()
st.header("Parameter sweep: does concentration or UV level grow longer chains?")
st.caption(
    "Runs a grid of independent simulations (pure hydrocarbon chemistry, no atmosphere -- "
    "this isolates the question to the elementary methane-photolysis system) varying "
    "starting concentration and UV level, holding the same real simulated exposure time "
    "fixed across every cell. That's the fair comparison: holding *event count* fixed "
    "instead makes higher concentration look like it suppresses growth, but that's an "
    "artifact -- a fixed step budget just covers a smaller fraction of a bigger starting "
    "pool. Uses the k_comb/k_abstr/complexity-ceiling from the sidebar. Larger grids and "
    "higher concentrations take longer (more molecules -> more events needed to cover the "
    "same real time) -- this can take up to a minute or so."
)

sweep_t_max = st.slider("Sweep: fixed exposure time across all cells", 1.0, 20.0, 5.0, step=1.0)
conc_options = st.multiselect(
    "Sweep: starting CH4 concentrations", [50, 100, 300, 800, 1500, 3000], default=[100, 300, 800],
)
kphoto_options = st.multiselect(
    "Sweep: UV levels (k_photo)", [0.005, 0.01, 0.05, 0.1, 0.2, 0.4], default=[0.01, 0.05, 0.2],
)
run_sweep = st.button("Run sweep", use_container_width=True)

if run_sweep:
    if not conc_options or not kphoto_options:
        st.warning("Pick at least one concentration and one UV level.")
    else:
        rows = []
        total = len(conc_options) * len(kphoto_options)
        progress = st.progress(0.0, text="Running sweep...")
        for i, conc in enumerate(sorted(conc_options)):
            for j, kp in enumerate(sorted(kphoto_options)):
                sweep_params = SimParams(
                    k_photo=kp, k_comb=k_comb, k_abstr=k_abstr,
                    max_carbon=max_carbon, t_max=sweep_t_max, max_events=2_000_000,
                    sample_every=1000, seed=int(seed),
                )
                sweep_sim = Simulator(sweep_params)
                sweep_sim.seed_species(methane(), conc)
                sweep_result = sweep_sim.run()
                stats = chain_length_stats(sweep_result)
                ch4_left = sweep_result.counts.get(methane().canonical_id(), 0) / conc
                rows.append(dict(
                    concentration=conc, k_photo=kp, mean_carbon=stats.mean_carbon,
                    max_carbon=stats.max_carbon_present, c3plus_frac=stats.c3plus_carbon_fraction,
                    ch4_remaining=ch4_left, events=len(sweep_result.event_log),
                    stopped=sweep_result.stopped_reason,
                ))
                progress.progress((i * len(kphoto_options) + j + 1) / total, text="Running sweep...")
        progress.empty()
        st.session_state["sweep_df"] = pd.DataFrame(rows)

if "sweep_df" in st.session_state:
    sweep_df = st.session_state["sweep_df"]
    pivot = sweep_df.pivot(index="concentration", columns="k_photo", values="mean_carbon")
    fig = go.Figure(go.Heatmap(
        z=pivot.values, x=[str(c) for c in pivot.columns], y=[str(c) for c in pivot.index],
        colorscale="Viridis", colorbar=dict(title="mean C#"), text=pivot.values,
        texttemplate="%{text:.2f}",
    ))
    fig.update_layout(
        xaxis_title="k_photo (UV level)", yaxis_title="starting CH4 concentration",
        height=400, xaxis=dict(type="category"), yaxis=dict(type="category"),
    )
    st.plotly_chart(fig, use_container_width=True)
    if any(sweep_df["stopped"] != "t_max reached"):
        st.warning(
            "Some sweep cells hit their event cap before reaching the target exposure time "
            "(see the 'stopped' column) -- their numbers understate how far that cell would "
            "actually progress. Try a smaller exposure time or fewer/lower concentrations."
        )
    st.dataframe(sweep_df, use_container_width=True, hide_index=True)

# -- Formose reaction: a path to ribose? --------------------------------------
st.divider()
st.header("Formose reaction: a path to ribose?")
st.caption(
    "A different chemistry regime entirely: formaldehyde self-condenses via base/mineral-"
    "catalyzed aldol addition into sugars, including ribose (C5H10O5) -- a required RNA "
    "building block, and arguably a fairer target than hydrocarbons for 'does life actually "
    "need this molecule have a plausible abiotic source.' The well-known catch (the 'sugar "
    "problem'): aldol addition and its reverse (retro-aldol) constantly scramble formaldehyde "
    "units in and out of sugars of every size, so real formose chemistry produces a mess of "
    "many sugars, not clean ribose -- and this model tracks sugars only by carbon count, so "
    "'C5 sugar' means the whole pentose pool (ribose + its stereoisomers), not ribose "
    f"specifically. See README for the {RIBOSE_FRACTION_ESTIMATE:.0%} literature-style "
    "estimate used to convert that pool size into a ribose-equivalent number."
)

fc1, fc2 = st.columns(2)
with fc1:
    hcho_count = st.slider("Starting formaldehyde (HCHO) count", 100, 5000, 1000, step=100)
    k_init = st.slider(
        "2 HCHO -> glycolaldehyde (k_init)", 0.0, 1.0, 0.01, step=0.01,
        help="Slow, rate-limiting induction step -- real formose famously has a slow start.",
    )
    k_aldol = st.slider("Aldol growth, Cn + HCHO -> C(n+1) (k_aldol)", 0.0, 10.0, 1.0, step=0.1)
    k_retro_aldol = st.slider(
        "Retro-aldol, C(n+1) -> Cn + HCHO (k_retro_aldol)", 0.0, 10.0, 0.3, step=0.1,
        help="Real formose equilibria run both ways -- this reversibility is the heart of the sugar problem.",
    )
    k_cannizzaro = st.slider(
        "Cannizzaro side reaction (k_cannizzaro)", 0.0, 2.0, 0.05, step=0.01,
        help="2 HCHO -> methanol + formate, a real dead end that competes for formaldehyde.",
    )
with fc2:
    max_sugar_carbon = st.slider("Max sugar size (complexity ceiling)", 4, 10, 8)
    stabilize_carbon = st.slider(
        "Which sugar size gets mineral protection", 2, max_sugar_carbon, 5,
        help="5 = ribose's carbon count.",
    )
    k_stabilize = st.slider(
        "Mineral (borate) stabilization rate (k_stabilize)", 0.0, 5.0, 0.0, step=0.05,
        help="Off (0) by default -- the 'no mineral rescue' case. Real borate minerals selectively "
             "bind and protect ribose's ring geometry from further reaction (Ricardo et al., Science "
             "2004). Turn this on to test whether that's actually what it takes.",
    )
    formose_t_max = st.slider("Formose: simulated time", 10.0, 2000.0, 300.0, step=10.0)
    formose_seed = st.number_input("Formose: random seed", value=7, step=1)

st.subheader("Is there a correction factor for ribose specifically?")
st.caption(
    "By default this model can't distinguish ribose from its 3 stereoisomers (arabinose, "
    "xylose, lyxose) -- 'C5 sugar' is the whole pentose pool. But real chemistry isn't fully "
    "uninformative here: aldol addition itself has no stereoselectivity (nothing to tune "
    "there), but borate's binding affinity is NOT uniform across pentose diastereomers -- "
    "ribose's furanose ring binds somewhat preferentially (Ricardo et al., Science 2004). "
    "Turning this on splits the target tier into the 4 real named diastereomers and lets "
    "ribose's stabilization rate be boosted relative to the other 3, so the ribose fraction "
    "becomes something the simulation actually produces instead of an imported guess."
)
track_stereo = st.checkbox(
    "Track pentose stereoisomers (ribose, arabinose, xylose, lyxose)", value=False,
    help="Only takes effect when 'Which sugar size gets mineral protection' is 5.",
)
ribose_selectivity = 1.0
if track_stereo:
    if stabilize_carbon != 5:
        st.warning(
            "Stereoisomer tracking only applies at carbon number 5 (ribose's carbon count) -- "
            "set 'Which sugar size gets mineral protection' to 5 above to use it."
        )
    ribose_selectivity = st.slider(
        "Ribose selectivity (borate binding preference vs. its 3 stereoisomers)",
        1.0, 20.0, 1.0, step=0.5,
        help="1.0 = no preference, matching real unselective aldol chemistry -- expect ribose "
             "to land around 1-in-4 of the stabilized pool. Raise it to see how much borate-"
             "binding selectivity it would take for ribose to dominate.",
    )

run_formose_clicked = st.button("Run formose simulation", use_container_width=True)

if run_formose_clicked:
    formose_params = FormoseParams(
        k_init=k_init, k_aldol=k_aldol, k_retro_aldol=k_retro_aldol, k_cannizzaro=k_cannizzaro,
        k_stabilize=k_stabilize, stabilize_carbon=stabilize_carbon, max_sugar_carbon=max_sugar_carbon,
        track_pentose_stereoisomers=track_stereo, ribose_selectivity=ribose_selectivity,
        t_max=formose_t_max, max_events=300_000, sample_every=1000, seed=int(formose_seed),
    )
    with st.spinner("Running formose simulation..."):
        formose_result = run_formose(formose_params, initial_hcho=hcho_count)
    st.session_state["formose_result"] = formose_result
    st.session_state["formose_params"] = formose_params
    st.session_state["formose_stabilize_carbon"] = stabilize_carbon
    st.session_state["formose_track_stereo"] = track_stereo and stabilize_carbon == 5

if "formose_result" in st.session_state:
    fresult = st.session_state["formose_result"]
    fparams = st.session_state["formose_params"]
    fstab_carbon = st.session_state["formose_stabilize_carbon"]
    ftrack = st.session_state["formose_track_stereo"]

    hcho_now = fresult.counts.get("HCHO", 0)

    if ftrack:
        stab_counts = {v: fresult.counts.get(f"C{fstab_carbon}sugar-stabilized-{v}", 0) for v in PENTOSE_DIASTEREOMERS}
        free_counts = {v: fresult.counts.get(f"C{fstab_carbon}sugar-{v}", 0) for v in PENTOSE_DIASTEREOMERS}
        total_stab = sum(stab_counts.values())
        ribose_share = stab_counts["ribose"] / total_stab if total_stab else 0.0

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("HCHO remaining", hcho_now)
        m2.metric("Stabilized ribose (actual count)", stab_counts["ribose"])
        m3.metric("Stabilized, whole pentose pool", total_stab)
        m4.metric("Ribose's share of stabilized pool", f"{ribose_share:.0%}")

        fig_stereo = go.Figure()
        fig_stereo.add_trace(go.Bar(name="Free (unprotected)", x=list(PENTOSE_DIASTEREOMERS),
                                     y=[free_counts[v] for v in PENTOSE_DIASTEREOMERS]))
        fig_stereo.add_trace(go.Bar(name="Stabilized (protected)", x=list(PENTOSE_DIASTEREOMERS),
                                     y=[stab_counts[v] for v in PENTOSE_DIASTEREOMERS]))
        fig_stereo.update_layout(barmode="group", yaxis_title="molecule count",
                                  xaxis_title="pentose diastereomer", height=400)
        st.plotly_chart(fig_stereo, use_container_width=True)

        if fparams.ribose_selectivity <= 1.0:
            st.info(
                f"**No selectivity applied**: ribose holds {ribose_share:.0%} of the stabilized pool, "
                "close to the 1-in-4 statistical baseline -- consistent with aldol addition having no "
                "stereochemical bias on its own. Raise ribose_selectivity to test whether differential "
                "borate binding is enough to change that."
            )
        else:
            st.success(
                f"**Selective stabilization works**: with ribose_selectivity={fparams.ribose_selectivity:g}, "
                f"ribose holds {ribose_share:.0%} of the stabilized pool ({stab_counts['ribose']} molecules) "
                f"versus {total_stab - stab_counts['ribose']} split across its 3 stereoisomers. This is a real, "
                "adjustable variable -- differential borate-binding affinity -- not a flat guess applied "
                "after the fact."
            )
    else:
        target_sid = f"C{fstab_carbon}sugar"
        stabilized_sid = f"C{fstab_carbon}sugar-stabilized"
        free_target = fresult.counts.get(target_sid, 0)
        stabilized_target = fresult.counts.get(stabilized_sid, 0)
        ribose_equivalent = round((free_target + stabilized_target) * RIBOSE_FRACTION_ESTIMATE, 1)

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("HCHO remaining", hcho_now)
        m2.metric(f"Free C{fstab_carbon} sugar", free_target)
        m3.metric(f"Stabilized C{fstab_carbon} sugar", stabilized_target)
        m4.metric(
            "Ribose-equivalent estimate", ribose_equivalent,
            help=f"({free_target} + {stabilized_target}) x {RIBOSE_FRACTION_ESTIMATE:.0%} -- a "
                 "literature-style heuristic, not something the simulation itself resolves. Turn on "
                 "stereoisomer tracking above to get an actual simulated ribose count instead.",
        )

        if fparams.k_stabilize == 0.0:
            st.info(
                f"**No mineral rescue**: {free_target} free C{fstab_carbon} sugar exists right now, but "
                "nothing protects it -- it stays in the same aldol/retro-aldol equilibrium as every other "
                "sugar size, so it doesn't durably accumulate. Turn on the stabilization rate to test "
                "whether that's what changes the picture."
            )
        elif stabilized_target > 0:
            st.success(
                f"**Mineral stabilization works**: {stabilized_target} C{fstab_carbon} sugar molecules have "
                f"been pulled out of the reactive pool and protected, versus only {free_target} unprotected "
                "free sugar sitting in the same equilibrium as everything else. This is the Ricardo et al. "
                "hypothesis in miniature: a rescue mechanism, not the aldol chemistry alone, is what lets a "
                "persistent ribose-sized pool build up."
            )
        else:
            st.warning("Stabilization is on but nothing has been protected yet -- try a higher rate or longer run.")

    ftab_pop, ftab_dist = st.tabs(["Sugar population over time", "Sugar-size distribution"])

    with ftab_pop:
        peak = {}
        for snap in fresult.history_counts:
            for sid, n in snap.items():
                peak[sid] = max(peak.get(sid, 0), n)
        top_ids = sorted(peak, key=peak.get, reverse=True)[:10]
        rows = []
        for t, snap in zip(fresult.history_t, fresult.history_counts):
            for sid in top_ids:
                rows.append({"t": t, "species": fresult.species[sid].formula(), "count": snap.get(sid, 0)})
        fdf = pd.DataFrame(rows)
        fig = go.Figure()
        for formula, sub in fdf.groupby("species"):
            fig.add_trace(go.Scatter(x=sub["t"], y=sub["count"], mode="lines", name=formula))
        fig.update_layout(
            yaxis_title="molecule count", xaxis_title="time (arbitrary units)",
            yaxis_type="log", height=450, legend_title="species (top 10 by peak count)",
        )
        st.plotly_chart(fig, use_container_width=True)

    with ftab_dist:
        def _sugar_count_at(n):
            if ftrack and n == fstab_carbon:
                return sum(fresult.counts.get(f"C{n}sugar-{v}", 0) for v in PENTOSE_DIASTEREOMERS)
            return fresult.counts.get(f"C{n}sugar", 0)

        sugar_counts = {f"C{n}": _sugar_count_at(n) for n in range(2, fparams.max_sugar_carbon + 1)}
        if fparams.k_stabilize > 0:
            if ftrack:
                sugar_counts[f"C{fstab_carbon} (stabilized)"] = sum(
                    fresult.counts.get(f"C{fstab_carbon}sugar-stabilized-{v}", 0) for v in PENTOSE_DIASTEREOMERS
                )
            else:
                sugar_counts[f"C{fstab_carbon} (stabilized)"] = stabilized_target
        fig2 = go.Figure(go.Bar(x=list(sugar_counts.keys()), y=list(sugar_counts.values())))
        fig2.update_layout(yaxis_title="molecule count", xaxis_title="sugar size", height=400)
        st.plotly_chart(fig2, use_container_width=True)
        if ftrack:
            st.caption(
                f"C{fstab_carbon} bars here sum across all 4 tracked diastereomers -- see the stereoisomer "
                "chart above for the ribose-specific breakdown. Every other sugar size is still a generic "
                "(CH2O)n pool."
            )
        else:
            st.caption(
                "Every sugar size is a generic (CH2O)n -- this doesn't distinguish ribose from its "
                "stereoisomers (arabinose, xylose, lyxose) or from ketopentoses. 'C5 sugar' is the whole "
                "pentose pool, not ribose specifically."
            )
