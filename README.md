# Abiogenesis: from methane photolysis to autocatalytic networks

A stochastic chemistry simulator exploring two linked questions:

1. Starting from UV photolysis of methane (CH4 + hv -> CH3• + H•), and
   letting the resulting radical chemistry grow in complexity on its own,
   does the reaction network ever produce **self-amplifying loops** -- the
   minimal structural prerequisite for Darwinian dynamics (differential
   reproduction) to have something to act on?
2. Add the rest of a plausible atmosphere -- N2, O2, CO2, Ar -- and expose
   it to the same UV: does the C2H6 that photochemistry produces actually
   **accumulate**, or does O2 (and the O3 it generates) outcompete radical
   self-combination and divert the carbon into oxidized products instead?

This is a toy model, not a calibrated photochemistry code. It's built to
let you see the *shape* of emergent chemical networks and the *qualitative*
outcome of chemical competitions, not to predict real atmospheric or
interstellar methane chemistry quantitatively.

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

## Atmosphere: N2, O2, CO2, Ar, and O3

Four more species can be seeded alongside methane, each modeled as a small
fixed species (not full carbon-graph molecules -- see Limitations) with
reactivity picked to match its real photochemistry qualitatively:

- **N2** and **Ar** never react anywhere in the model. This is accurate at
  the level this toy model operates: N2's triple bond (~9.8 eV) needs far
  harder UV than what's driving C-H photolysis or O2's O=O bond (~5.2 eV)
  here, and Ar is a noble gas. Both are pure spectators/diluents.
- **O2** photolyzes to two O atoms (`O2 + hv -> 2 O*`), and -- the
  mechanism the whole atmosphere question hinges on -- directly **scavenges
  hydrocarbon radicals**: `R* + O2 -> ROO*` (a peroxy radical). This
  reaction competes with radical-radical self-combination
  (`R* + R'* -> R-R'`, e.g. `CH3* + CH3* -> C2H6`) for the *same* radical
  pool. Both rate constants default to the same order of magnitude on
  purpose, so which pathway wins is decided by relative **concentration**
  (how much O2 is actually present), not by a thumb on the scale.
- **O3 (ozone)** forms from `O* + O2 -> O3` (a simplified bimolecular
  stand-in for the real termolecular `O + O2 + M -> O3 + M` -- no
  third-body/pressure dependence here), photolyzes back to `O2 + O*`, and
  is itself highly reactive: `R* + O3 -> RO* (alkoxy radical) + O2` is a
  second, independent radical sink alongside direct O2 scavenging, and it
  regenerates O2 as it consumes O3.
- **CO2** photolyzes to CO + O*. CO is treated as an unreactive terminal
  product here, which sidesteps the bond-order reorganization (the
  remaining C=O double bond becoming a C-O triple bond) that a fully
  accurate graph treatment of CO2 photolysis would require.

O and OH radicals produced along the way plug into the *existing* generic
abstraction rule for free (anything with `is_radical=True` can pull an H off
any hydrocarbon), and their combination products (alkoxy radicals, alcohols,
hydroperoxides, water) are all lightweight wrapper species that tag a
hydrocarbon fragment with a fixed O-group rather than full graph
molecules -- no further oxidation is modeled past that first O-addition
step (see Limitations).

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

A default run (methane seed, `max_carbon=6`, `t_max=400`, no atmosphere)
reliably finds ~20-25 candidate cycles, the strongest of which is typically
the simple ethane/ethyl-radical interconversion (`C2H6 <-> C2H5*`) firing
hundreds of times -- a legitimate, if simple, chain-propagation loop. Larger
`max_carbon` and longer runs surface longer, rarer cycles among
higher-carbon isomers.

**Does C2H6 accumulate once O2 is in the mix?** In a controlled test (same
seed, same rates, only O2 abundance varied -- see
`test_o2_competes_with_self_combination_for_c2h6` in `tests/test_engine.py`):
with no O2, ~150 self-combination events fire and C2H6 accumulates freely.
With O2 abundant relative to the trace CH3• pool, self-combination fires
only 5 times against 300 O2-scavenging events, and **final C2H6 count is
zero** -- the carbon instead ends up in peroxy radicals and their
hydroperoxide products. Adding realistic O2 photolysis and ozone chemistry
on top (`test_full_atmosphere_ethane_accumulation`) only reinforces this:
ozone builds up (over 1600 formation events in one run) and adds a second,
independent radical sink. This is a genuine emergent result of the
simulation, not a hand-tuned outcome -- the O2-scavenging and
self-combination rate constants are deliberately set to the *same order of
magnitude*, so the competition is decided by concentration, exactly as it
would be in reality: even if a collision with another CH3• and a collision
with an O2 molecule were equally likely to react, there are vastly more O2
molecules around to collide with. This is consistent with why prebiotic
chemistry (Miller-Urey and successors) is generally thought to require a
*reducing*, low-O2 atmosphere -- an oxidizing one burns the organic
chemistry back down about as fast as photochemistry can make it.

## Limitations (read before drawing conclusions)

- **No rings, no *general* heteroatoms.** The core carbon-skeleton
  representation is still acyclic hydrocarbons only. N2/O2/CO2/O3 and the
  first-generation O-oxidation products (alkoxy/peroxy radicals, alcohols,
  hydroperoxides, water) are bolted on as small fixed/wrapper species (see
  Atmosphere section above), not full graph molecules -- there's no
  representation yet of, say, a carbon chain with *two* different O-groups
  on it, or O incorporated into the carbon skeleton itself (ethers,
  carbonyls-in-a-chain). That's still the fuller "heteroatoms" roadmap item
  below.
- **Oxidation chemistry stops after one O-addition.** Peroxy/alkoxy
  radicals, alcohols, and hydroperoxides don't get further photolyzed,
  don't abstract H themselves, and (mostly) don't combine with each other
  (`combine()` returns `None` for e.g. two peroxy radicals meeting -- see
  `engine/molecule.py`). Real atmospheric oxidation chains go much further
  (aldehydes, further fragmentation, eventually CO2 + H2O). This model only
  needs to show the *first* branch point (self-combination vs. scavenging)
  to answer the accumulation question, so it stops there.
- **Atomic O is simplified to a monoradical.** Real ground-state O has two
  unpaired electrons; modeled here as one reactive site, consistent with
  how H* is modeled. OH is treated as terminal (nothing abstracts its H
  back off).
- **No third-body/pressure-dependent kinetics.** O3 formation
  (`O* + O2 -> O3`) is modeled as simple bimolecular; the real reaction is
  termolecular (`O + O2 + M -> O3 + M`) and pressure-dependent. O2/N2/Ar
  don't create any general pressure or third-body effect in this model.
- **No true catalysis (in the hydrocarbon-only cycle-detection sense).**
  Nothing in the *hydrocarbon* reaction-network engine has its rate
  enhanced by the presence of another species. The cycles detected by
  `engine/autocatalysis.py` are chain-*propagation* loops (a radical
  carrier gets regenerated), which is real and relevant but a weaker claim
  than "autocatalytic set" in the Kauffman sense.
- **No heredity.** Darwinian evolution needs variation that is copied with
  error (so selection has something to act on across generations). Nothing
  here stores or transmits sequence information; there's no polymer, no
  template-directed replication. This model tests whether the *substrate*
  (a complex, cyclic reaction network) is reachable from simple starting
  chemistry, and whether a given product's *concentration* is favored or
  disfavored by the surrounding chemistry -- not a claim that replication
  itself emerges.
- **Rates are relative, not physical.** Don't read absolute timescales or
  yields as predictions about real methane/atmosphere photochemistry.
  Concentration *ratios* and which pathway dominates are the meaningful
  outputs, not absolute numbers.
- **Well-mixed, not spatial.** No compartments, gradients, or interfaces
  (e.g. mineral surfaces, lipid vesicles) that real origin-of-life chemistry
  likely depended on for concentrating and protecting reactive intermediates.
  Also no altitude-dependent UV/pressure profile, so this can't reproduce
  something like a real stratospheric ozone layer.

## Roadmap toward Darwinian relevance

Roughly in order of how directly each one advances the actual question:

1. **General heteroatoms in the graph.** O2/CO2/O3/N2/Ar and first-step
   oxidation products are in now (see Atmosphere section), but as fixed/
   wrapper species bolted onto the carbon-only graph, not as O/N nodes
   *inside* it. Generalizing `Molecule` to carry element type and bond
   order per node/edge (valence table `{C:4, N:3, O:2, H:1}`, H-count still
   implicit) would let carbonyls, ethers, amines, and multi-step oxidation
   chains fall out of the *same* rule engine instead of hand-written wrapper
   classes -- unlocking a much richer, more Miller-Urey-like reaction space
   and real catalytic possibilities (e.g. a species that lowers another
   reaction's activation energy).
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
  molecule.py       carbon-skeleton graph representation, seed molecules,
                     and the fixed/wrapper atmospheric + oxidation species
  reactions.py       the reaction rules + Reaction/propensity model
  simulator.py       Gillespie SSA with on-demand reaction-network growth
  autocatalysis.py   realized-flow graph + candidate cycle detection
app.py                Streamlit UI
tests/test_engine.py  sanity checks (no pytest dependency)
```
