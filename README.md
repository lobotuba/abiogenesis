# Abiogenesis: from methane photolysis to autocatalytic networks

A stochastic chemistry simulator exploring three linked questions:

1. Starting from UV photolysis of methane (CH4 + hv -> CH3• + H•), and
   letting the resulting radical chemistry grow in complexity on its own,
   does the reaction network ever produce **self-amplifying loops** -- the
   minimal structural prerequisite for Darwinian dynamics (differential
   reproduction) to have something to act on?
2. Add the rest of a plausible atmosphere -- N2, O2, CO2, Ar, H2O -- and
   expose it to the same UV: does the C2H6 that photochemistry produces
   actually **accumulate**, or does O2 (direct, or self-generated from
   water) outcompete radical self-combination and divert the carbon into
   oxidized products instead?
3. Add electricity as a distinct energy source: is there a real chemical
   difference between UV and electric discharge, or is "electricity" just
   UV with a different label? (There's a concrete answer: N2's triple bond
   is a wall UV in this model can never get through, but a spark can.)

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

## Water photolysis, hydrogen escape, and abiotic O2

A wet planet doesn't need free O2 seeded to start building its own: H2O
photolyzes too (`H2O + hv -> H* + OH*`). On its own this goes nowhere --
the primary products have no path back to O2/O3 -- so a second real
reaction is included: **`OH* + OH* -> H2O + O*`** (hydroxyl
disproportionation), the secondary step that actually bootstraps a free O
atom out of pure water photochemistry. From there the existing O2/O3
machinery (see above) takes over.

**But a closed system never builds up any free O2/O3 at all**, no matter
how much UV or how much simulated time -- confirmed in testing out to
1000+ time units and nearly a million reaction events: every O atom
liberated from H2O finds its way back into H2O via the same radical
chemistry (abstraction, combination). This matches real planetary
chemistry: the process needs a one-way loss of hydrogen to actually leave
oxygen behind, since hydrogen is the lightest gas and the one real
atmospheres actually lose to space. That's modeled here as **hydrogen
escape** -- `H2 -> (nothing)` and `H* -> (nothing)`, a simple first-order
sink with no altitude/exobase physics, off by default (opt-in via the
sidebar's `k_escape` slider). This is the accepted mechanism proposed for
*abiotic* O2/O3 buildup on a wet, lifeless planet -- studied in the
literature as a potential "false positive" biosignature for exoplanets
(a planet could show O2/O3 in its spectrum with zero biology involved).

**What we found running it:** turning escape on does let free O2/O3
persist, but in this model it reaches a modest, *self-limited* steady
state rather than an unbounded runaway -- O3's own photolysis
(`O3 -> O2 + O*`) provides a reverse flux that grows right along with O3's
forward production, the same balance that maintains a roughly steady
concentration in the real stratospheric ozone layer. Across escape rates
from 0.1 up to 100 (and simulated times up to 3000 units), starting from
300 CH4 + 300 H2O, O3 plateaued around 60-70 molecules -- present, but not
enough of a *standing* population to seriously outcompete hydrocarbon
self-combination for the transient radical pool (scavenging events stayed
near zero in most of these runs, versus hundreds to thousands of
self-combination events). Compare that to the atmosphere experiments
above, where O2 was seeded directly at *thousands* of molecules (a modern
Earth-like ratio) and scavenging dominated overwhelmingly -- the difference
is standing abundance, not mechanism. So: **on this toy model's own terms,
a wet planet exposed only to UV does not obviously drive prebiotic
hydrocarbon chemistry to ~0% -- the oxidant it can self-generate this way
plateaus well short of the seeded-O2-rich scenario**, unless escape and/or
water-photolysis rates are pushed well past what's realistic (both are
adjustable in the UI, so you can push past that line yourself and watch
the crossover happen). This is consistent with the real open question in
astrobiology: exactly how much abiotic O2 a given planet can accumulate
this way depends sensitively on escape rate, UV flux, and geological
timescales most toy models (including this one) don't fully capture --
see Limitations.

## Electricity vs UV as an energy source

Miller-Urey's original experiment used electric spark discharge, not UV,
and that choice was load-bearing, not incidental: **N2's triple bond
(~9.8 eV) is beyond what UV in this model's implied range can break, but a
spark's energy density is high enough to crack it.** That's precisely why
discharge could pull nitrogen into the reducing gas mixture (toward HCN
and, eventually, amino acids) while pure UV photolysis of CH4/NH3/H2/H2O
mostly can't reach N2 at all. UV also tends to be far more *wavelength-
selective* (governed by each bond's absorption cross-section) than a spark,
which is closer to indiscriminately energetic across whatever bonds are
present. On the other hand, UV was almost certainly the *larger total
energy source* on early Earth's surface/upper atmosphere (lightning is
intense but rare and localized; solar UV is diffuse but enormously more
abundant in integrated flux) -- both mechanisms are real and complementary
in the literature, not competing explanations.

**Electricity is now a first-class variable** (`k_discharge_n2`, off by
default): `N2 + spark -> 2 N*`, generated as its own "discharge" reaction
kind, deliberately separate from photolysis. This is the one thing in the
whole model that UV categorically cannot do at any intensity -- confirmed
directly in testing (`test_electricity_cracks_n2_but_uv_alone_never_does`):
with discharge off, N2 sits at its exact seeded count no matter how much UV
or simulated time runs; with discharge on, it visibly depletes. Once N*
exists it plugs into the *existing* combination/abstraction machinery for
free -- no new rules needed there, only new species: N* + N* -> N2
(recombination), N* + H* -> NH* (imidogen), and N* (or NH*) + a hydrocarbon
radical -> a closed-shell amine `R-NH2` in one step. That last one
deliberately skips a valence-correct aminyl-radical intermediate (nitrogen
is trivalent; a bond + a single "radical site" the way this model already
treats O doesn't divide as cleanly for N) -- see `engine/molecule.py` for
the reasoning. In one test run (300 CH4 + 300 N2, discharge on), nitrogen
took over the system: NH3 became the single most abundant product, methane/
ethane/propane were fully converted, and amines (C1H3NH2 through C4H9NH2)
accounted for most of the remaining carbon -- a small-scale, qualitative
echo of what discharge does that UV alone cannot. No N-O cross chemistry
(real NOx chemistry) is modeled -- a deliberate scope cut to keep this
addition focused on the one question it was built to answer.

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

## Chain length: what actually drives longer hydrocarbons?

The app's "Chain-length distribution" tab shows the carbon-number histogram
for a single run, and the "Parameter sweep" section (`engine/analysis.py` +
the sweep block in `app.py`) runs a grid of independent simulations varying
starting CH4 concentration and UV level (`k_photo`) to answer this directly.

**The one thing that matters when comparing runs: hold real simulated time
fixed, not event count.** Comparing at a fixed number of reaction events
makes higher concentration look like it *suppresses* chain growth -- but
that's an artifact: a fixed step budget just covers a smaller fraction of a
bigger starting pool, since propensities (and therefore how much simulated
time each event represents) scale with concentration. The sweep tool holds
`t_max` fixed across every cell for this reason.

With that fair comparison:

- **UV level has a strong, robust effect.** At fixed concentration and fixed
  exposure time, raising `k_photo` reliably produces longer chains and a
  higher C3+ carbon fraction (regression-tested in
  `test_higher_uv_yields_longer_chains_at_fixed_time`). Mechanistically:
  more UV sustains a higher steady-state radical population, and
  radical-radical combination (the *only* chain-growing step) scales with
  `[R*]^2`, so it becomes disproportionately more likely relative to
  non-growing H-abstraction as the radical pool grows.
- **Starting concentration alone, surprisingly, barely matters** when
  compared at the same real exposure time -- mean chain length stayed
  roughly flat across a 50x concentration range (100 to 5000) in testing.
  Higher concentration produces *far more total reaction events* in the
  same time window (bimolecular rates scale with concentration²), but the
  resulting *distribution shape* (average chain length, % converted) ends
  up depending mainly on `k_photo * t` and the rate-constant ratios, not on
  the absolute concentration scale. Caveat: in the tested range, every cell
  reached the `max_carbon` complexity ceiling almost immediately, so this
  can't rule out concentration mattering more once chains are allowed to
  grow past that ceiling -- raise `max_carbon` and re-run the sweep to
  check.

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
- **No full HOx chemistry.** Water photolysis's secondary oxygen-bootstrap
  path stops at one reaction (`OH* + OH* -> H2O + O*`). Real atmospheric
  HOx chemistry also includes H2O2 and HO2, which this model doesn't track.
- **Escape has no altitude/exobase physics.** Hydrogen escape is a flat
  first-order rate applied uniformly to H2/H*, not a function of where in
  an atmosphere they are, temperature, or stellar XUV flux -- real
  atmospheric escape is governed by all of these.
- **Nitrogen chemistry stops at one amine-forming step.** N* + hydrocarbon
  radical goes straight to a closed-shell amine (`R-NH2`), skipping a
  valence-correct aminyl-radical intermediate (see Electricity section
  above). No N-O cross chemistry (NOx) is modeled, and no path toward
  nitriles/HCN -- the real Miller-Urey route to amino acid precursors goes
  further than this model does. NH is treated as terminal, same
  simplification as OH.
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
   oxidation/amination products are in now (see Atmosphere and Electricity
   sections), but as fixed/wrapper species bolted onto the carbon-only
   graph, not as O/N nodes *inside* it. Generalizing `Molecule` to carry
   element type and bond order per node/edge (valence table
   `{C:4, N:3, O:2, H:1}`, H-count still implicit) would let carbonyls,
   ethers, amines with a *correct* aminyl-radical intermediate, and
   multi-step oxidation chains fall out of the *same* rule engine instead
   of hand-written wrapper classes -- unlocking a much richer, more
   Miller-Urey-like reaction space and real catalytic possibilities (e.g. a
   species that lowers another reaction's activation energy).
2. **Nitrile / HCN chemistry.** The real Miller-Urey route to amino acid
   precursors goes through HCN and related nitriles, not just amines --
   this model's nitrogen chemistry stops one step short of that.
3. **Rings.** Lift the tree-only restriction; rings enable a different
   class of stable, potentially catalytic structures.
4. **Explicit catalysis.** Let some species appear as a rate multiplier on
   a reaction (not just a reactant/product) so the autocatalysis detector
   can test the real Kauffman condition, not just chain-propagation flux.
5. **RAF (reflexively autocatalytic, food-generated) set analysis.** Once
   catalysis is explicit, implement the actual RAF algorithm instead of the
   current cycle-flux heuristic -- this is the rigorous version of what
   `engine/autocatalysis.py` currently approximates.
6. **Polymers + templating.** The real jump: a backbone chemistry (even a
   toy one) where sequence can be copied with occasional error. This is
   where "Darwinian" stops being a stretch and starts being literal --
   variation + heredity + differential survival.
7. **Spatial structure.** Compartments or surfaces that let useful
   combinations of molecules stay together instead of diluting into a
   well-mixed soup. Also where altitude-dependent escape/UV physics would
   belong, if the wet-planet question above needs to get more realistic.

## Project layout

```
engine/
  molecule.py       carbon-skeleton graph representation, seed molecules,
                     and the fixed/wrapper atmospheric + oxidation/amination species
  reactions.py       the reaction rules + Reaction/propensity model
  simulator.py       Gillespie SSA with on-demand reaction-network growth
  autocatalysis.py   realized-flow graph + candidate cycle detection
  analysis.py         chain-length distribution stats (single-run tab + sweep tool)
app.py                Streamlit UI
tests/test_engine.py  sanity checks (no pytest dependency)
```
