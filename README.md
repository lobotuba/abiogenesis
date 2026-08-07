# Abiogenesis: from methane photolysis to autocatalytic networks

A stochastic chemistry simulator exploring a narrow, concrete question:
starting from UV photolysis of methane (CH4 + hv -> CH3• + H•), and letting
the resulting radical chemistry grow in complexity on its own, does the
reaction network ever produce **self-amplifying loops** -- the minimal
structural prerequisite for Darwinian dynamics (differential reproduction)
to have something to act on?

This is a toy model, not a calibrated photochemistry code. It's built to
let you see the *shape* of emergent chemical networks, not to predict real
atmospheric or interstellar methane chemistry quantitatively.

## The model

**Molecules** are represented as carbon-skeleton graphs (a tree of C-C
bonds; no rings in v1). Hydrogen counts are derived, not stored: each
carbon has valence 4, so H count = 4 - (skeleton degree) - (1 if it
currently carries the molecule's unpaired electron). At most one radical
site per molecule is modeled -- this keeps the chemistry to ordinary
mono-radical organic chemistry and avoids carbene/biradical bookkeeping.
Atomic H and H2 are handled as two small fixed species alongside the
carbon-graph molecules.

**Three reaction rules** are applied generically to whatever species exist
in the pool, so the network's complexity grows on its own rather than being
hand-authored:

1. **Photolysis** `M + hv -> M* + H*` -- UV breaks one C-H bond homolytically.
   Rate weighted by the number of chemically-equivalent C-H bonds (e.g.
   methane's 4 equivalent H's).
2. **Combination** `R* + R'* -> R-R'` -- two radicals join to form a new C-C
   (or H-H, or C-H) bond.
3. **Abstraction** `R* + M-H -> R-H + M*` -- a radical pulls an H off a
   closed-shell molecule, regenerating a *different* radical. This is the
   important one: it's how a radical can be **regenerated** a few steps
   after being consumed, which is the real chemical mechanism behind
   radical chain reactions -- the closest thing to autocatalysis available
   in pure hydrocarbon chemistry.

Rate constants (`k_photo`, `k_comb`, `k_abstr`, adjustable in the UI) are
relative, not calibrated to real photon flux or Arrhenius kinetics. A crude
bond-strength proxy biases abstraction/photolysis toward more-substituted
carbons (secondary/tertiary > primary), giving the network a
realistic-flavored preference for branched products instead of a uniform
random walk over isomers.

**Simulation** uses the Gillespie stochastic simulation algorithm (SSA) over
a reaction network that grows on demand: the first time a new molecule
appears, the three rules generate whatever new reactions it enables against
the existing species pool; that reaction list is then reused for the rest of
the run (propensities recompute from current counts every step, but the
candidate list itself doesn't need re-deriving each event). A configurable
`max_carbon` acts as a complexity ceiling so the combinatorial isomer
explosion of organic chemistry doesn't run away.

**Autocatalysis detection** (`engine/autocatalysis.py`) builds a directed
graph from only the reactions that *actually fired* during a run (reactant
-> product edges, weighted by fire count), then searches for short directed
cycles. A short, high-flux cycle -- e.g. `C2H6 -> C2H5• -> C2H6` -- is a
candidate chain-propagation loop: a species that gets regenerated through a
different path a few reactions after being consumed. This is a heuristic
diagnostic, not a rigorous test (true autocatalysis requires showing a
species's presence *causes* its own faster production, which this doesn't
prove) -- see Limitations below.

## Running it

```
pip install -r requirements.txt
streamlit run app.py
```

Or run the non-interactive sanity checks:

```
python tests/test_engine.py
```

## What we've already seen

A default run (methane seed, `max_carbon=6`, `t_max=400`) reliably finds
~20-25 candidate cycles, the strongest of which is typically the simple
ethane/ethyl-radical interconversion (`C2H6 <-> C2H5*`) firing hundreds of
times -- a legitimate, if simple, chain-propagation loop. Larger
`max_carbon` and longer runs surface longer, rarer cycles among
higher-carbon isomers.

## Limitations (read before drawing conclusions)

- **No rings, no heteroatoms.** Only acyclic hydrocarbons and their
  radicals. Real prebiotic chemistry needs O, N (and eventually P) for
  anything resembling nucleotides, amino acids, or lipids.
- **No true catalysis.** Nothing in this model has its *rate* enhanced by
  the presence of another species (which is what "catalysis" actually
  means kinetically). The cycles detected are chain-*propagation* loops
  (a radical carrier gets regenerated), which is a real and relevant
  phenomenon, but is a weaker claim than "autocatalytic set" in the
  Kauffman sense.
- **No heredity.** Darwinian evolution needs variation that is copied with
  error (so selection has something to act on across generations). Nothing
  here stores or transmits sequence information; there's no polymer, no
  template-directed replication. This model is a test of whether the
  *substrate* (a complex, cyclic reaction network) is reachable from simple
  starting chemistry -- not a claim that replication itself emerges.
- **Rates are relative, not physical.** Don't read absolute timescales or
  yields as predictions about real methane photochemistry.
- **Well-mixed, not spatial.** No compartments, gradients, or interfaces
  (e.g. mineral surfaces, lipid vesicles) that real origin-of-life chemistry
  likely depended on for concentrating and protecting reactive intermediates.

## Roadmap toward Darwinian relevance

Roughly in order of how directly each one advances the actual question:

1. **Heteroatoms.** Add O (from H2O/CO2 photochemistry) and N to reach
   carbonyls, alcohols, amines -- unlocks a much richer, more
   Miller-Urey-like reaction space and real catalytic possibilities
   (e.g. a species that lowers another reaction's activation energy).
2. **Rings.** Lift the tree-only restriction; rings enable a different
   class of stable, potentially catalytic structures.
3. **Explicit catalysis.** Let some species appear as a rate multiplier on
   a reaction (not just a reactant/product) so the autocatalysis detector
   can test the real Kauffman condition, not just chain-propagation flux.
4. **RAF (reflexively autocatalytic, food-generated) set analysis.** Once
   catalysis is explicit, implement the actual RAF algorithm instead of the
   current cycle-flux heuristic -- this is the rigorous version of what
   `engine/autocatalysis.py` currently approximates.
5. **Polymers + templating.** The real jump: a backbone chemistry (even a
   toy one) where sequence can be copied with occasional error. This is
   where "Darwinian" stops being a stretch and starts being literal --
   variation + heredity + differential survival.
6. **Spatial structure.** Compartments or surfaces that let useful
   combinations of molecules stay together instead of diluting into a
   well-mixed soup.

## Project layout

```
engine/
  molecule.py       carbon-skeleton graph representation + seed molecules
  reactions.py       the three reaction rules + Reaction/propensity model
  simulator.py       Gillespie SSA with on-demand reaction-network growth
  autocatalysis.py   realized-flow graph + candidate cycle detection
app.py                Streamlit UI
tests/test_engine.py  sanity checks (no pytest dependency)
```
