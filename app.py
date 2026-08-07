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

from engine.autocatalysis import build_flow_graph, find_candidate_cycles
from engine.molecule import SEED_MOLECULES
from engine.simulator import SimParams, Simulator

st.set_page_config(page_title="Abiogenesis: hydrocarbon photochemistry", layout="wide")

st.title("Methane photochemistry -> reaction-network complexity")
st.caption(
    "CH4 + UV -> CH3• + H•, then radical combination and H-abstraction let the "
    "network grow on its own. This is a toy, non-calibrated model: the interesting "
    "question isn't absolute rates, it's whether the *shape* of the reaction network "
    "(chain-propagation cycles) gives anything for selection to act on."
)

with st.sidebar:
    st.header("Setup")
    seed_name = st.selectbox("Starting molecule", list(SEED_MOLECULES.keys()), index=0)
    seed_count = st.slider("Starting molecule count", 10, 2000, 300, step=10)

    st.header("Rate constants (relative units)")
    k_photo = st.slider("UV photolysis (k_photo)", 0.001, 1.0, 0.05, step=0.001, format="%.3f")
    k_comb = st.slider("Radical combination (k_comb)", 0.1, 50.0, 5.0, step=0.1)
    k_abstr = st.slider("H-abstraction (k_abstr)", 0.01, 10.0, 0.5, step=0.01)

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
        max_carbon=max_carbon, t_max=t_max, max_events=max_events,
        sample_every=max(1, max_events // 500), seed=int(seed),
    )
    sim = Simulator(params)
    sim.seed_species(SEED_MOLECULES[seed_name](), seed_count)
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

tab_pop, tab_net, tab_cycles, tab_species = st.tabs(
    ["Population dynamics", "Reaction network", "Autocatalytic cycles", "Species inventory"]
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
