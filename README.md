# Abiogenesis: from methane photolysis to autocatalytic networks

A stochastic chemistry simulator exploring eight linked questions:

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
4. Set hydrocarbons aside and ask a fairer question: RNA (and so all life
   on Earth) needs ribose. Does the classic prebiotic route to it --
   formaldehyde self-condensing into sugars via the formose reaction --
   actually let a persistent ribose-sized sugar pool build up, or does it
   need something else? (See the Formose section below -- this is a
   genuinely different chemistry regime, and the honest answer involves a
   real, famous problem with formose chemistry.)
5. Push past ribose to an actual nucleotide, using a route that never
   needs free ribose to exist at all: attach the nucleobase precursor to
   the sugar chain *before* the ribose-defining stereocenters even form.
   Does a photochemical selection step actually matter for how much
   nucleotide forms, and if so, *why* -- is it really about rescuing
   material, or about something else? (See the Nucleotide section below --
   the honest answer surprised the model's own author partway through
   building it.)
6. Push past a single monomer to an actual strand, and ask the question
   this whole project started with -- self-amplification -- again, now
   that there's finally something strand-like to copy: does a minimal
   template-directed self-replicator (modeled on von Kiedrowski's real
   1986 self-replicating hexanucleotide) actually amplify itself, or does
   the same chemistry that lets it copy itself also cap how far that can
   go? (See the Polymerization section below -- the real answer is
   famous, and it's not the one "autocatalysis" makes it sound like.)
7. Make the founding question literal: given two heritable replicator
   sequences competing for the same finite food supply, does the
   faster-copying one actually win -- and can that faster variant even
   *arise* from copying error in the first place, rather than being
   seeded in from outside? This needs all three ingredients Darwinian
   selection actually requires (heritable variation, a way for copying to
   produce that variation, and differential survival under a shared
   constraint) at once, for the first time in this project. (See the
   Selection section below.)
8. Ask the opposite question: instead of competing, what if replicators
   **cooperate**? In a closed loop where each one catalyzes replication of
   the *next* one (Eigen & Schuster's hypercycle), can members too weak to
   replicate on their own persist and grow anyway -- purely because the
   loop is closed? (See the Hypercycle section below -- the test is about
   as sharp as a stochastic model can give.)

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

**Does electricity produce O2 in a wet environment, and does that stop
other molecules from existing?** Discharge doesn't touch water or oxygen
chemistry at all -- it only cracks N2. But nitrogen fixation consumes H*
(`N* + H* -> NH*`, `NH* + H* -> NH3`), and that's the *same* H* that
normally recombines with OH* back into H2O. Starve that recombination and
OH* lives long enough to instead disproportionate (`OH* + OH* -> H2O +
O*`), feeding O3. Tested directly in a **closed system with hydrogen
escape off** (so any O2/O3 that appears can only be this side effect, not
escape-driven buildup): with discharge off, O2+O3 stays at exactly zero,
same as any closed wet system; with discharge on, O2+O3 becomes nonzero
purely from this hydrogen-sink effect
(`test_electricity_indirectly_enables_o2_via_hydrogen_sink`). But the
second half of the question has a counterintuitive answer: **the O2/O3
that appears is not what suppresses other molecules.** In the same test
run, O2/O3-scavenging fired zero times -- 13 molecules of O3 is nowhere
near enough standing concentration to compete for the radical pool (it
takes thousands, not tens -- see the Atmosphere section above). What
actually captured 98.7% of all carbon into amines was nitrogen fixation
*directly*, competing for the same hydrocarbon-radical pool that would
otherwise self-combine. The app's "Three-way competition" panel (shown
whenever N2 is present) makes this explicit: it compares self-combination
vs. O2/O3-scavenging vs. amine-forming event counts side by side, plus a
carbon-accounting bar chart (every tracked carbon atom is in exactly one
of: plain hydrocarbons, amines, or O-containing products) so you can see
directly which pathway is actually responsible for what you're seeing,
rather than assuming it's whichever byproduct is easiest to notice.

## Formose reaction: a path to ribose

Everything above is radical chemistry: UV or a spark breaks a bond
homolytically, unpaired electrons pair back up. Ribose (C5H10O5, required
for RNA) doesn't come from that pathway in reality -- it comes from the
**formose reaction**, formaldehyde self-condensing into sugars via
base/mineral-catalyzed **aldol addition**, a closed-shell ionic mechanism
with no radicals anywhere. That's different enough chemistry that it gets
its own module (`engine/formose.py`) and its own section of the app,
rather than being bolted onto the radical engine.

**The reaction network**, all generated upfront (unlike the hydrocarbon
engine, the formose species set is small and fully known in advance, so
there's no need for on-demand discovery):

1. **Induction**: `2 HCHO -> C2 sugar` (glycolaldehyde) -- slow and
   rate-limiting, matching formose's real slow start.
2. **Aldol growth**: `Cn + HCHO -> C(n+1)` -- the autocatalytic
   chain-growing step, for every sugar size up to a `max_sugar_carbon`
   complexity ceiling (same role as `max_carbon` elsewhere).
3. **Retro-aldol**: `C(n+1) -> Cn + HCHO` -- the reverse of (2). Real
   formose equilibria run both ways, and this reversibility is the heart
   of the famous **"sugar problem"**: aldol addition and its reverse
   constantly scramble formaldehyde in and out of sugars of every size, so
   a real formose reaction produces a mess of many different sugars
   ("formose tar"), not clean ribose.
4. **Cannizzaro side reaction**: `2 HCHO -> waste` (methanol + formate) --
   a real, well-known dead end that competes for formaldehyde and limits
   how much of it is ever available for sugar growth at all.
5. **Mineral stabilization** (optional, off by default): `Cn -> stabilized
   Cn` at a chosen carbon number (default 5). This is the proposed
   resolution in the literature (Ricardo, Carrigan, Olcott, Benner,
   ["Borate Minerals Stabilize Ribose"](https://www.science.org/doi/10.1126/science.1102722),
   *Science* 2004): borate ions selectively bind ribose's ring geometry
   and protect it from further reaction, pulling it out of the
   aldol/retro-aldol equilibrium before it's scrambled into hexoses or
   degraded. Modeled as a one-way, terminal reaction -- the same
   "protected from further chemistry" treatment already used elsewhere in
   this project for e.g. H2O or NH3.

**What we found, testing directly**
(`test_formose_ribose_needs_stabilization`): without stabilization, the
free C5 sugar pool stays small (9 molecules, in one representative run) --
consistent with the sugar problem, since nothing stops it from continuing
to react further or reverting. With stabilization turned on, a
**persistent protected pool of 73** builds up instead. The central
hypothesis holds up in this simplified model: **ribose-sized sugar needs a
rescue mechanism to accumulate; the aldol chemistry alone won't durably
produce it.**

### Is ribose's stereoisomer problem fixable, or just estimated away?

The model's biggest limitation was that "C5 sugar" meant the whole pentose
pool -- ribose, arabinose, xylose, and lyxose all lumped together, with a
flat `RIBOSE_FRACTION_ESTIMATE` (1/4) applied after the fact to guess how
much of that pool was "really" ribose. That's not fully honest: real aldol
addition genuinely has no stereoselectivity without a chiral catalyst (so
treating the four diastereomers as equally likely going *into* the pool is
correct), but real borate-binding affinity is *not* uniform across them --
ribose's furanose ring binds somewhat preferentially. So there is a real
variable hiding behind that flat estimate, and it belongs on the
*stabilization* step specifically, not smeared across the whole pool.

`track_pentose_stereoisomers` (off by default) turns this on: the C5 tier
splits into the four real named diastereomers (`PENTOSE_DIASTEREOMERS` in
`engine/formose.py`), aldol addition produces each with equal weight (the
honest, unselective part), and `ribose_selectivity` multiplies *only*
ribose's stabilization rate relative to the other three. Tested directly
(`test_formose_ribose_selectivity_is_a_real_correction_factor`): with
`ribose_selectivity=1.0`, ribose holds ~33% of the stabilized pool (close
to the 1-in-4 baseline, within stochastic noise -- confirming the
unselective-aldol assumption is behaving as intended); with
`ribose_selectivity=5.0`, that jumps to ~57%, with ribose's absolute count
more than doubling while the other three stereoisomers stay flat. The
ribose fraction is now something the simulation actually produces from an
adjustable, chemically-motivated variable, not a number typed in after the
fact.

**A different, more radical real answer, now built as its own module
(see the Nucleotide section below):** Sutherland et al.'s 2009 *Nature*
paper ("Synthesis of activated pyrimidine ribonucleotides in prebiotically
plausible conditions") sidesteps the sugar problem entirely rather than
solving it -- ribose and the nucleobase are synthesized *together* from
glycolaldehyde and cyanamide, with unwanted diastereomers destroyed by UV
along the way, so free ribose never needs to exist and be purified on its
own. That's a fundamentally different pathway, not a parameter tweak.

**Still not modeled**: full 3D stereochemistry beyond the one tracked
tier (real formose stereochemistry starts mattering from C3 onward, not
just C5), crossed aldol reactions between two sugars (only sugar +
formaldehyde growth is included, a scope cut that avoids an M×N
combinatorial explosion of product sizes while still capturing the
essential growth mechanism), and any borate-independent stabilization
pathway.

## Nucleotide pathway: skipping the sugar problem entirely

The formose module's answer to "does ribose accumulate" was "only with a
rescue mechanism." This module asks a sharper version of the same
question by using a route that never produces free ribose to begin with
(`engine/nucleotide.py`, Powner, Gerland & Sutherland, *Nature* 2009):
attach the nucleobase precursor to the sugar chain *before* the
ribose-defining stereocenters even form, so there's no free-sugar mixture
to purify ribose out of in the first place.

**The reaction network**, again fully known upfront:

1. Glycolaldehyde + cyanamide -> 2-aminooxazole. Glycolaldehyde is
   literally `Sugar(2)` imported from `engine/formose.py` -- this pathway
   branches off the *same* earliest formose intermediate, before the sugar
   problem has any chance to develop.
2. 2-aminooxazole + glyceraldehyde (`Sugar(3)`) -> one of 4 diastereomeric
   aminooxazolines (ribo-/arabino-/xylo-/lyxo-configured), formed with
   equal likelihood -- reusing the exact same "unselective aldol-type
   step, real stereoselectivity enters later" logic as the formose
   module's stereoisomer tracking.
3. **Photochemical selection** (off by default): the real, documented
   Sutherland finding is that ribo-configured aminooxazoline is
   photostable while its 3 stereoisomers are destroyed by UV at the
   wavelength used. Modeled as a photolysis reaction touching only the 3
   non-ribo variants.
4. Ribo-aminooxazoline + cyanoacetylene -> anhydronucleoside (productive).
   The 3 non-ribo diastereomers, if not removed by step 3, react with
   cyanoacetylene at the same rate -- just unproductively, into a dead-end
   adduct.
5. Anhydronucleoside + phosphate -> an activated pyrimidine ribonucleotide.

**What we found, testing directly -- including a wrong first guess,
caught by testing it**: the obvious hypothesis is that photochemical
selection helps by "rescuing" material -- destroy the losers, leave more
raw material for the winner. Building it that way first
(no reaction 4 competition, just destruction) produced *no difference at
all* in final nucleotide yield with selection on vs. off, which is what
exposed the actual mechanism: photolysis doesn't run the aldol step
backwards, so destroying a diastereomer doesn't free up any
glyceraldehyde or aminooxazole to make more of the productive one. The
branching ratio into ribo- is fixed the moment the aminooxazoline forms
(rule 2); nothing downstream changes that number.

What selection *actually* does, once reaction 4's competition for
cyanoacetylene is included: it protects a **scarce downstream resource**
from being wasted on unproductive branches. Cyanoacetylene is genuinely
the special, scarce reagent in this pathway (seeded lower than the sugar
precursors in the app's defaults, matching its real prebiotic scarcity
relative to simple sugars). Tested directly
(`test_nucleotide_forms_and_photoselection_protects_scarce_cyanoacetylene`):
with no selection, cyanoacetylene splits ~1-in-4 to the productive branch
like anything else, yielding 121 ribonucleotide in one run; with fast
selection (`k_photo_destroy` large relative to
`k_anhydro x [cyanoacetylene]` -- selection has to *win the kinetic race*,
not just be present), yield triples to 300, with the non-ribo diastereomers
destroyed before they get a chance to consume any cyanoacetylene at all.
**The mechanism is protection of a scarce resource, not recycling of a
diluted one** -- a real distinction this project only found by building
the wrong version first and noticing the result didn't move.

**Not modeled**: the exact real-world numeric yields, the intermediate
anhydronucleoside/hydrolysis chemistry Sutherland's actual synthesis needs
several more steps for, and (same limitation as formose) real 3D
stereochemistry -- these are named diastereomer labels, not geometry.

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

## Polymerization: does copying itself actually amplify?

The nucleotide module stops at one finished, activated monomer. This
module (`engine/polymer.py`) asks the question this whole project started
with -- self-amplification -- one more time, now that there's finally
something strand-like to copy. It's modeled directly on von Kiedrowski,
"A Self-Replicating Hexadeoxynucleotide" (*Nature*, 1986): the first
experimentally demonstrated minimal template-directed self-replicator.

**The reaction network**, same generic-species, no-explicit-sequence
abstraction level as the nucleotide module:

1. Oligomerization: 3 activated ribonucleotides (`nucleotide.RIBONUCLEOTIDE`,
   imported directly -- this pathway starts exactly where that module ends)
   spontaneously join into a FRAGMENT (a trimer). Deliberately slow/weak by
   default: condensation reactions release water, and in bulk solution
   that makes uncatalyzed ligation thermodynamically uphill -- the real
   difficulty behind why biological polymers need help forming at all.
2. Background ligation: two FRAGMENTs join into a TEMPLATE (a hexamer),
   equally slow and uncatalyzed -- the "no templating exists yet" baseline.
3. **Templated ligation, the point of this module**: two FRAGMENTs plus an
   *existing* TEMPLATE ligate into two TEMPLATEs. The old TEMPLATE isn't
   consumed -- it's a genuine catalyst, so this reaction is autocatalytic
   by construction. TEMPLATE is modeled as self-complementary (its own two
   halves are each a copy of FRAGMENT) -- one species acting as substrate,
   catalyst, and product simultaneously, which is the actual trick behind
   von Kiedrowski's real minimal self-replicator.
4. **Duplex formation, the other point of this module**: the same
   self-complementarity that lets one TEMPLATE catalyze step 3 also lets
   two TEMPLATEs hybridize with EACH OTHER into an inert double-stranded
   DUPLEX -- a dead end until it melts back apart. This isn't a penalty
   bolted on from outside; it's mechanistically the same base-pairing
   chemistry as step 3, just acting on two product strands instead of a
   product and two fragments.
5. Duplex melting: DUPLEX -> 2 TEMPLATE, a slow reverse (thermal
   denaturation / wet-dry or day-night cycling in the real system),
   reintroducing catalytically active single strands.

**This needed a real three-reactant reaction** (step 3: two FRAGMENTs and
a TEMPLATE all at once), which the engine didn't support --
`Reaction.propensity` in `engine/reactions.py` special-cased exactly one
or two reactants by hand. It was generalized to a single formula: for each
distinct reactant id appearing with multiplicity *m* in a reaction, the
number of ways to draw *m* molecules of it from the current pool is
`comb(n, m)`; multiply those across all distinct ids. This is a strict
generalization -- it reduces to the exact old `n`, `n_a * n_b`, and
`n_a*(n_a-1)/2` arithmetic for every reaction already in the project
(verified directly in `test_reaction_propensity_generalizes_to_n_body`),
and now supports genuinely termolecular reactions too.

**What we found, tested directly**: turning on templated ligation (step 3)
measurably speeds up TEMPLATE accumulation relative to background ligation
alone -- in one test, background ligation alone reached 26 TEMPLATE in a
fixed time window while the templated pathway reached 102 in the same
window (`test_polymer_templating_accelerates_strand_formation`). That's
real autocatalysis, not just a label for it.

But turning on duplex formation (step 4) reveals von Kiedrowski's actual,
famous, counterintuitive result: in one test run
(`test_polymer_duplex_formation_causes_self_inhibition`), duplex-off and
duplex-on both fully converted the same 1000-FRAGMENT pool into
502 template-equivalent units -- but with duplex off, all 502 stayed
active, single-stranded TEMPLATE; with duplex on, **0** stayed active and
251 ended up locked away as inert duplex. Duplex formation doesn't reduce
how much ligation chemistry happens -- it's mass-conserving, verified
exactly by `test_polymer_mass_conservation` (1 ribonucleotide-unit per
monomer, 3 per fragment, 6 per template, 12 per duplex, constant at every
sampled point of a run) -- it just locks up the catalyst as fast as the
catalyst makes more of itself. This is the real reason minimal
template-directed self-replicators like von Kiedrowski's show *parabolic*
(roughly square-root-of-time) growth instead of exponential/Malthusian
growth: the very self-pairing that lets a strand copy itself is what
increasingly sequesters it out of solution as its own concentration rises.
So the honest answer to this project's founding question, now that there's
finally something strand-like to copy, is: yes, template-directed
self-amplification is real and measurable here -- but it is not
unconditional. The same mechanism that enables it structurally caps it.

**Not modeled here**: explicit sequence/base-pairing (FRAGMENT and TEMPLATE
are tracked only by length, the same abstraction nucleotide.py uses), a
mismatch/error-rate model for copying fidelity, and anything past a single
self-complementary replicator. Those three gaps are exactly what the
Selection section below adds.

## Selection: does a faster replicator actually win?

This is the founding question of the whole project, made literal at last.
Darwinian selection strictly needs three ingredients at once: (1)
**heritable variation** -- more than one distinguishable replicating type,
(2) a copying process that can **produce** that variation from an existing
type (mutation), and (3) **differential survival/reproduction** under a
shared constraint. `engine/selection.py` adds all three on top of
polymer.py's already-tested chemistry, at the same "generic species, no
explicit base sequence" abstraction level formose.py and nucleotide.py use
for their diastereomers -- named variant tags ("A", "B") stand in for
distinct heritable sequences.

**The reaction network**, per named variant, reusing `Oligomer`/`Duplex`
from `engine/polymer.py` (both extended with an optional `variant` tag,
backward-compatible -- `None` preserves polymer.py's original single-
species ids exactly):

1. Oligomerization and background ligation (3 RIBONUCLEOTIDE -> FRAGMENT_v,
   FRAGMENT_v + FRAGMENT_v -> TEMPLATE_v) use the SAME rate constants for
   every variant -- deliberately unbiased, so every lineage has equal, fair
   access to the shared food supply. This is what makes it a genuine
   competition rather than a rigged one.
2. **Templated ligation is where fitness actually lives**: FRAGMENT_v +
   FRAGMENT_v + TEMPLATE_v -> 2 TEMPLATE_v, at each variant's OWN rate.
   This is the one place a "better" sequence is allowed to matter, matching
   the real biophysical claim: fitness differences between self-replicators
   come from how well each one catalyzes its own copying, not from how
   easily its raw material forms.
3. Duplex formation/melting, same self-inhibition mechanism as
   polymer.py, shared across variants. Deliberately no cross-variant
   duplexes (TEMPLATE_A + TEMPLATE_B -> hybrid) -- the simplifying
   assumption that distinct heritable sequences are different enough not to
   cross-hybridize, the same role sequence *specificity* plays in real
   template-directed replication.
4. **Mutation** (off by default): during templated ligation, an existing
   TEMPLATE_v occasionally produces a TEMPLATE_w of a *different* tracked
   variant instead of a faithful copy of itself. This is the actual source
   of heritable variation the Polymerization section's own "not modeled"
   note flagged as missing: with mutation on (and background ligation off,
   for a clean comparison -- otherwise undirected background chemistry
   could also independently stumble into forming any variant from scratch),
   a variant that was never seeded can still appear purely through
   imperfect copying of an existing one.

**What we found, tested directly**:

- *Direct competition* (both variants pre-seeded with identical stock, no
  shared-resource confound): giving one variant 3x the templated-ligation
  rate of the other, from identical starting conditions, produced a
  decisive lead (46 vs 8 in one test,
  `test_selection_fitter_variant_wins_direct_competition`) mid-transient --
  measured before either exhausted its fixed fragment stock, since (same
  lesson as polymer.py's own acceleration test) given infinite time both
  eventually fully convert their own stock regardless of rate.
- *The control*: with variant B never seeded and mutation off, B stayed at
  **exactly 0** across every seed tested -- confirmed there is no hidden
  backdoor for it to appear
  (`test_selection_mutation_is_required_for_unseeded_variant_to_appear`).
- *The capstone result*: seed only variant A, give a never-seeded variant B
  a 10x higher templated-ligation rate, and turn mutation on. In one test
  run, B -- which starts at a strict double disadvantage (doesn't exist yet,
  and can only appear once a mutation event happens to fire) -- not only
  came into existence (via 101 mutation events) but **overtook** the
  already-established, pre-seeded variant A by the time the shared
  ribonucleotide supply ran out (1729 vs 1606,
  `test_selection_fitter_mutant_overtakes_established_lineage`), robust
  across every seed checked during calibration. This is the literal minimal
  signature of Darwinian selection this whole project set out to test for:
  a heritable variant that arose from copying error alone, and then won,
  because it was fitter -- not because it was seeded, not because it was
  favored by an external rule, just because it copied itself faster than
  its competitor given the same shared, finite resource.

**Not modeled here**: real nucleotide-level sequences (still just named
tags, not base-pairing geometry), more than a handful of competing
variants at once, cooperation between distinct replicators (hypercycles --
see the next section, which adds exactly that), and spatial structure
(well-mixed only, so there's no way for a locally-successful variant to
outrun diffusion into a shared global pool) -- see Roadmap.

## Hypercycle: does cooperation let the weak survive?

Selection (above) is replicators *competing*, each relying entirely on its
own copying ability. `engine/hypercycle.py` builds the other classic
arrangement from origin-of-life theory: Eigen & Schuster's **hypercycle**
(Eigen, *Naturwissenschaften*, 1971; Eigen & Schuster, *The Hypercycle: A
Principle of Natural Self-Organization*, 1979) -- a closed loop of
replicators where each one catalyzes replication of the *next* one in the
cycle, not itself. The classic point of a hypercycle is that it lets
replicators too weak to sustain themselves alone persist and grow *only
because the loop is closed*.

**The reaction network**, per named variant at position `i` in a list of
variants, reusing `Oligomer`/`Duplex` from `engine/polymer.py` exactly like
`engine/selection.py` does:

1. Oligomerization and background ligation, shared/unbiased across
   variants -- same as selection.py.
2. **Self-templated ligation, deliberately off for every member by
   default** (`k_self`, per variant, defaults to 0.0): the classic
   hypercycle setup studies replicators that are individually *too weak to
   be self-sufficient*, so whatever growth happens has to come from
   cooperation, not solo copying.
3. **Cross-catalyzed ligation, the actual hypercycle step**: variant `v`'s
   TEMPLATE catalyzes ligation of the *next* variant's own fragments into a
   new copy of that neighbor (FRAGMENT_next + FRAGMENT_next + TEMPLATE_v ->
   TEMPLATE_v + TEMPLATE_next) -- structurally the same three-reactant
   reaction shape as `templated_ligation` (self-copying, polymer.py) and
   `mutation` (error-copying, selection.py), just pointed at a neighbor
   instead of at itself or at random.
4. Duplex formation/melting, same self-inhibition mechanism as the other
   replicator modules, shared across variants.

A `closed` flag controls whether the *last* variant's cross-catalyzed
reaction wraps back around to the *first* variant (a true closed loop,
A -> B -> C -> A) or is simply omitted (an open chain, A -> B -> C, where
nothing on Earth ever catalyzes A). This one flag is the entire experiment:
same species, same rates, same everything else -- does closing that one
edge change what the system can do?

**What we found, tested directly**: with every member's self-templating
rate at 0 (nobody can replicate alone) and background ligation also at 0
(so a TEMPLATE's count can only ever change by being produced as a
reaction *product*), variant A in an open 3-member chain is not just
"slower to grow" -- it is **mathematically incapable of changing at all**,
since no reaction in the whole open-chain network has it as a product.
Tested directly across every seed tried
(`test_hypercycle_closing_the_loop_lets_an_unhelped_member_grow`): open
chain, A stays at *exactly* its seeded value (5) every single time.
Closing the loop is the only structural change, and in one test run it let
A grow to 1135 -- right alongside B (1107) and C (1105), all three rising
together, none of them capable of doing it alone. This is about as sharp a
result as a stochastic model can produce: not "closing the loop helps," but
"closing the loop is the *only* way this particular growth can happen at
all." Mass conservation holds throughout
(`test_hypercycle_mass_conservation`), confirming the loop is redistributing
a shared, finite ribonucleotide supply among all three members, not
creating material from nowhere.

**Not modeled**: more than a simple ring topology (real hypercycle theory
also studies branched and higher-connectivity networks), parasites (a
"cheater" species that gets catalyzed by the cycle but doesn't catalyze
anyone in return -- the classic hypercycle *vulnerability*, not just its
strength), and, same as Selection above, real base-level sequences and
spatial structure.

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
- **Heredity, competition, and cooperation -- but still not real
  sequences.** `engine/polymer.py` has genuine template-directed
  replication, `engine/selection.py` adds heritable variation, mutation,
  and demonstrated differential survival under a shared resource
  constraint, and `engine/hypercycle.py` adds cooperative cross-catalysis
  (replicators too weak to sustain themselves alone, persisting only
  because a loop is closed) -- see the Selection and Hypercycle sections
  above. What's still missing: variant tags aren't real base sequences (no
  explicit A/U/G/C, no mismatch/binding geometry), only a handful of
  variants are tracked at once (not an open-ended population), and there's
  no "cheater" species that exploits a hypercycle without contributing to
  it -- the classic hypercycle vulnerability, not just its strength. See
  Roadmap item 3.
- **Rates are relative, not physical.** Don't read absolute timescales or
  yields as predictions about real methane/atmosphere photochemistry.
  Concentration *ratios* and which pathway dominates are the meaningful
  outputs, not absolute numbers.
- **Well-mixed, not spatial.** No compartments, gradients, or interfaces
  (e.g. mineral surfaces, lipid vesicles) that real origin-of-life chemistry
  likely depended on for concentrating and protecting reactive intermediates.
  Also no altitude-dependent UV/pressure profile, so this can't reproduce
  something like a real stratospheric ozone layer.
- **Stereochemistry is only tracked at one tier, opt-in.** With
  `track_pentose_stereoisomers` on, the formose module distinguishes
  ribose from its 3 stereoisomers at C5 specifically (see the Formose
  section above) -- but that's the only place in the whole project with
  any stereochemistry at all, and it stops at naming 4 diastereomers with
  adjustable rates, not real 3D geometry. Off by default, "C5 sugar" is
  still the whole pentose pool and `RIBOSE_FRACTION_ESTIMATE` is a carried-
  in heuristic. No crossed sugar-sugar aldol reactions (only sugar +
  formaldehyde growth), and mineral stabilization is a single one-way
  reaction, not a real borate-binding equilibrium.

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
2. **Deeper stereochemistry for the formose module.** Named-diastereomer
   tracking with a selective stabilization rate now exists at the C5 tier
   (`track_pentose_stereoisomers`), which was the most direct fix for the
   ribose question specifically -- but it's still a label with an
   adjustable rate, not real 3D/chirality geometry, and it only applies at
   one carbon number. Real geometry (so the model could in principle
   *derive* rather than assume which sugars borate should prefer) and
   tracking through more of the carbon ladder (real stereochemistry starts
   at C3) are the natural next steps. The bigger structural alternative,
   sidestepping free-ribose synthesis entirely, is now built -- see the
   Nucleotide pathway section above.
3. **Beyond the activated ribonucleotide, competition, and cooperation --
   all three now built.** `engine/polymer.py` polymerizes the activated
   ribonucleotide into a self-complementary, template-copying strand (von
   Kiedrowski's 1986 minimal self-replicator); `engine/selection.py` adds
   heritable variation, mutation, and tested differential survival under a
   shared resource constraint; `engine/hypercycle.py` adds cooperative
   cross-catalysis, showing replicators too weak to sustain themselves
   alone can persist and grow together purely because a loop is closed
   (see Polymerization, Selection, and Hypercycle sections above). This is
   where this project's two founding questions -- self-amplifying networks
   (item 1's original goal) and a plausible route to ribose -- finally met,
   and where "Darwinian" stopped being a stretch and became literal: a
   fitter, never-seeded variant was shown to emerge from copying error
   alone and overtake an established competitor, and a closed cooperative
   loop was shown to be the *only* way a particular member could grow at
   all. What's still missing, and is now the most direct remaining path
   forward: real base-level sequences (not just named tags), an open-ended
   population of variants instead of a handful, and a "cheater" species
   that exploits a hypercycle without contributing to it -- the classic
   vulnerability of cooperation, and a natural place for this project's two
   modes (competition and cooperation) to finally interact directly.
4. **Nitrile / HCN chemistry.** The real Miller-Urey route to amino acid
   precursors goes through HCN and related nitriles, not just amines --
   this model's nitrogen chemistry stops one step short of that.
5. **Rings.** Lift the tree-only restriction; rings enable a different
   class of stable, potentially catalytic structures.
6. **Explicit catalysis.** Let some species appear as a rate multiplier on
   a reaction (not just a reactant/product) so the autocatalysis detector
   can test the real Kauffman condition, not just chain-propagation flux.
7. **RAF (reflexively autocatalytic, food-generated) set analysis.** Once
   catalysis is explicit, implement the actual RAF algorithm instead of the
   current cycle-flux heuristic -- this is the rigorous version of what
   `engine/autocatalysis.py` currently approximates.
8. **Spatial structure.** Compartments or surfaces that let useful
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
  formose.py           separate aldol-chemistry module: formaldehyde -> sugars -> ribose
  nucleotide.py         separate module: glycolaldehyde/glyceraldehyde + cyanamide/
                         cyanoacetylene/phosphate -> activated ribonucleotide, skipping
                         free ribose entirely (imports Sugar from formose.py)
  polymer.py             separate module: ribonucleotide -> fragment -> a
                         self-complementary, template-copying strand (von Kiedrowski
                         1986 minimal self-replicator; imports RIBONUCLEOTIDE from
                         nucleotide.py)
  selection.py            separate module: adds heritable variant tags, mutation
                         (copying error), and tested differential survival on top of
                         polymer.py's chemistry (imports Oligomer/Duplex from
                         polymer.py, RIBONUCLEOTIDE from nucleotide.py) -- the literal
                         minimal Darwinian-selection test
  hypercycle.py           separate module: cross-catalytic closed-loop cooperation
                         between distinct replicators (Eigen and Schuster's hypercycle;
                         imports Oligomer/Duplex from polymer.py, RIBONUCLEOTIDE from
                         nucleotide.py) -- tests whether a closed loop lets replicators
                         too weak to self-sustain persist and grow together
app.py                Streamlit UI
tests/test_engine.py  sanity checks (no pytest dependency)
```
