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
from engine.molecule import ATMOSPHERIC_GASES, Molecule, SEED_MOLECULES, ethane, methane
from engine.simulator import SimParams, Simulator

st.set_page_config(page_title="Abiogenesis: hydrocarbon photochemistry", layout="wide")

st.title("Methane photochemistry -> reaction-network complexity")
st.caption(
    "CH4 + UV -> CH3• + H•, then radical combination and H-abstraction let the "
    "network grow on its own. Add N2/O2/CO2/Ar/H2O to test the central question: does "
    "the C2H6 that forms actually *accumulate*, or does O2 (direct, or self-generated from "
    "water photolysis) outcompete radical self-combination and divert it into oxidized "
    "products instead? This is a toy, non-calibrated model -- see README for what it does "
    "and doesn't capture."
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
        "N2 and Ar never react in this model (accurately: N2's triple bond needs far "
        "harder UV than this sim implies; Ar is a noble gas). O2 and CO2 photolyze; O2 "
        "additionally *scavenges* hydrocarbon radicals directly, and separately builds up "
        "O3 (ozone) with UV-generated O atoms -- O3 is itself highly reactive and scavenges "
        "radicals too. H2O also photolyzes (H2O -> H* + OH*); two OH* can disproportionate "
        "into H2O + O*, which is how a wet planet can bootstrap its *own* O2/O3 from water "
        "alone, no free O2 required. See the callout below the run for how this plays out."
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
