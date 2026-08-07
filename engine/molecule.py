"""
Graph-based representation of saturated hydrocarbons and hydrocarbon radicals.

A molecule's carbon skeleton is a tree (no rings in v1): nodes are carbon
atoms, edges are C-C bonds. Each carbon has valence 4; its hydrogen count is
derived (not stored) as 4 minus its skeleton degree minus 1 if it currently
carries the molecule's unpaired electron ("radical site"). At most one
radical site per molecule is supported in v1 -- this keeps the model to
mono-radical organic chemistry (the dominant regime in low-pressure gas-phase
photolysis) and avoids carbene/biradical bookkeeping.

Atomic hydrogen (H*) and H2 are handled as small fixed species since they
have no carbon skeleton.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import networkx as nx


@runtime_checkable
class Species(Protocol):
    """Common interface shared by Molecule and the atomic hydrogen species."""

    def canonical_id(self) -> str: ...
    def formula(self) -> str: ...
    @property
    def n_carbon(self) -> int: ...
    @property
    def is_radical(self) -> bool: ...
    def radicalizable_sites(self) -> list["Site"]: ...


@dataclass(frozen=True)
class Site:
    """A candidate place to strip one H and place a radical there.

    multiplicity: number of chemically-equivalent H atoms at this position
        (e.g. methane's carbon has multiplicity 4).
    product: the Species that results from removing one H at this site.
    degree: skeleton degree of the carbon losing the H (0 for atomic/H2);
        used as a crude proxy for C-H bond strength (primary < secondary <
        tertiary reactivity) when weighting abstraction rates.
    """

    multiplicity: int
    product: "Species"
    degree: int


class Molecule:
    __slots__ = ("graph", "_formula", "_hash", "_radical_node")

    def __init__(self, graph: nx.Graph):
        self.graph = graph
        self._formula: str | None = None
        self._hash: str | None = None
        rad_nodes = [n for n, d in graph.nodes(data=True) if d.get("radical")]
        if len(rad_nodes) > 1:
            raise ValueError("v1 model supports at most one radical site per molecule")
        self._radical_node = rad_nodes[0] if rad_nodes else None

    # -- basic properties -------------------------------------------------

    @property
    def n_carbon(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def is_radical(self) -> bool:
        return self._radical_node is not None

    def h_count(self, node) -> int:
        deg = self.graph.degree[node]
        rad = 1 if node == self._radical_node else 0
        return 4 - deg - rad

    @property
    def total_h(self) -> int:
        return sum(self.h_count(n) for n in self.graph.nodes)

    def formula(self) -> str:
        if self._formula is None:
            c, h = self.n_carbon, self.total_h
            marker = "•" if self.is_radical else ""
            self._formula = f"C{c}H{h}{marker}"
        return self._formula

    def canonical_id(self) -> str:
        if self._hash is None:
            self._hash = nx.weisfeiler_lehman_graph_hash(
                self.graph, node_attr="radical", iterations=3
            )
        return self._hash

    def __eq__(self, other):
        return isinstance(other, Molecule) and self.canonical_id() == other.canonical_id()

    def __hash__(self):
        return hash(self.canonical_id())

    def __repr__(self):
        return f"Molecule({self.formula()})"

    # -- reaction-relevant queries -----------------------------------------

    def radicalizable_sites(self) -> list[Site]:
        """All (multiplicity, product, degree) for stripping one H from a carbon.

        Skips the carbon that already carries the radical (no biradicals in
        v1) and skips carbons with no remaining H.
        """
        sites = []
        for node in self.graph.nodes:
            if node == self._radical_node:
                continue
            mult = self.h_count(node)
            if mult <= 0:
                continue
            new_graph = self.graph.copy()
            for n in new_graph.nodes:
                new_graph.nodes[n]["radical"] = n == node
            sites.append(Site(multiplicity=mult, product=Molecule(new_graph), degree=self.graph.degree[node]))
        return sites

    def saturate_radical(self) -> "Molecule":
        """Return the closed-shell molecule formed by adding H at the radical site."""
        if self._radical_node is None:
            raise ValueError("molecule has no radical site to saturate")
        new_graph = self.graph.copy()
        new_graph.nodes[self._radical_node]["radical"] = False
        return Molecule(new_graph)

    @staticmethod
    def join(a: "Molecule", b: "Molecule") -> "Molecule":
        """Form a new C-C bond between a's and b's radical sites."""
        if a._radical_node is None or b._radical_node is None:
            raise ValueError("both molecules must have a radical site to combine")
        g = nx.disjoint_union(a.graph, b.graph)
        offset = a.graph.number_of_nodes()
        a_node = a._radical_node
        b_node = offset + b._radical_node
        g.add_edge(a_node, b_node)
        g.nodes[a_node]["radical"] = False
        g.nodes[b_node]["radical"] = False
        return Molecule(g)


class _AtomicH:
    """Monatomic hydrogen radical, H*."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "H"

    def formula(self) -> str:
        return "H•"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return True

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _AtomicH)

    def __hash__(self):
        return hash("H")

    def __repr__(self):
        return "Molecule(H•)"


class _H2:
    """Molecular hydrogen, H2. Not directly UV-photolyzed in this model;
    only participates via radical H-abstraction."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "H2"

    def formula(self) -> str:
        return "H2"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        # Breaking the H-H bond leaves a lone H atom; multiplicity 2 because
        # either H is an equivalent abstraction target.
        return [Site(multiplicity=2, product=ATOMIC_H, degree=0)]

    def __eq__(self, other):
        return isinstance(other, _H2)

    def __hash__(self):
        return hash("H2")

    def __repr__(self):
        return "Molecule(H2)"


ATOMIC_H = _AtomicH()
MOLECULAR_H2 = _H2()


# -- atmospheric gases -------------------------------------------------
#
# These are modeled as small fixed species rather than extending the
# carbon-skeleton graph to arbitrary heteroatoms/bond orders (that fuller
# generalization is still on the roadmap). Each one's *reactivity* is
# hand-picked to match its real photochemistry at the qualitative level
# this toy model operates at:
#
#   Ar   -- noble gas, zero reaction rules attached: pure spectator/diluent.
#   N2   -- triple bond, dissociation energy ~9.8 eV, far higher than the UV
#           that drives C-H photolysis or O2's ~5.2 eV O=O bond in this
#           model. No photolysis rule attached: N2 is treated as inert
#           under the "actinic" UV this sim implies, same as it is in
#           reality below the far-UV (<100 nm).
#   O2   -- photolyzes to two O atoms, and is the one new species that
#           *actively competes* with hydrocarbon-radical self-combination
#           (see the "o2_scavenging" rule in reactions.py): R* + O2 -> ROO*.
#           This competition is the mechanism that decides whether C2H6
#           accumulates or gets diverted into oxidized products. O2 also
#           combines with O atoms to form O3 (ozone), which photolyzes back
#           to O2 + O* and separately scavenges hydrocarbon radicals of its
#           own (R* + O3 -> RO* + O2, rule "ozone_scavenging") -- a second,
#           independent radical sink alongside direct O2 scavenging.
#   CO2  -- photolyzes to CO + O. CO is treated as an unreactive terminal
#           product in this model (real CO is comparatively inert at these
#           conditions); this sidesteps the bond-order reorganization
#           (C=O double bond -> C-O triple bond) that a fully accurate
#           graph-based treatment of CO2 photolysis would require.
#
# O atoms and OH radicals produced along the way plug into the *existing*
# generic abstraction rule for free (anything with is_radical=True can pull
# an H off any hydrocarbon with radicalizable_sites()) -- no extra code
# needed there. What needed new code is what they *become*: O2, atomic H,
# and hydroxyl combine into water; O atoms and OH combine with hydrocarbon
# radicals into alkoxy radicals / alcohols; peroxy radicals combine with H
# into hydroperoxides. See combine() below.


class _AtomicO:
    """Atomic oxygen, O*. Real ground-state O is a diradical; simplified to
    a single reactive site here, consistent with how H* is modeled."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "O"

    def formula(self) -> str:
        return "O•"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return True

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _AtomicO)

    def __hash__(self):
        return hash("O")

    def __repr__(self):
        return "Molecule(O•)"


class _Hydroxyl:
    """Hydroxyl radical, OH*."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "OH"

    def formula(self) -> str:
        return "OH•"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return True

    def radicalizable_sites(self) -> list[Site]:
        return []  # nothing abstracts the H back off OH in this model

    def __eq__(self, other):
        return isinstance(other, _Hydroxyl)

    def __hash__(self):
        return hash("OH")

    def __repr__(self):
        return "Molecule(OH•)"


class _O2:
    """Dioxygen. Splits into 2 O* under UV (see reactions.py); also the
    radical scavenger that competes with hydrocarbon self-combination."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "O2"

    def formula(self) -> str:
        return "O2"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _O2)

    def __hash__(self):
        return hash("O2")

    def __repr__(self):
        return "Molecule(O2)"


class _O3:
    """Ozone. Forms from O* + O2 (see reactions.py), photolyzes back to
    O2 + O* under UV, and -- the reason it matters here -- reacts fast with
    hydrocarbon radicals (R* + O3 -> RO* + O2), giving a second pathway
    (alongside direct O2 scavenging) that pulls radicals away from forming
    C2H6. Modeled as a simple bimolecular O* + O2 -> O3 step rather than the
    real termolecular O + O2 + M -> O3 + M (no third-body/pressure
    dependence in this toy model)."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "O3"

    def formula(self) -> str:
        return "O3"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _O3)

    def __hash__(self):
        return hash("O3")

    def __repr__(self):
        return "Molecule(O3)"


class _N2:
    """Dinitrogen. No reaction rule attached anywhere -- inert spectator."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "N2"

    def formula(self) -> str:
        return "N2"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _N2)

    def __hash__(self):
        return hash("N2")

    def __repr__(self):
        return "Molecule(N2)"


class _CO2:
    """Carbon dioxide. Photolyzes to CO + O* (see reactions.py)."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "CO2"

    def formula(self) -> str:
        return "CO2"

    @property
    def n_carbon(self) -> int:
        return 0  # deliberately excluded from the hydrocarbon max_carbon accounting

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _CO2)

    def __hash__(self):
        return hash("CO2")

    def __repr__(self):
        return "Molecule(CO2)"


class _CO:
    """Carbon monoxide. Treated as an unreactive terminal product here."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "CO"

    def formula(self) -> str:
        return "CO"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _CO)

    def __hash__(self):
        return hash("CO")

    def __repr__(self):
        return "Molecule(CO)"


class _Argon:
    """Argon. Pure spectator: no photolysis, no bonding, ever."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "Ar"

    def formula(self) -> str:
        return "Ar"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _Argon)

    def __hash__(self):
        return hash("Ar")

    def __repr__(self):
        return "Molecule(Ar)"


class _H2O:
    """Water. Terminal/unreactive in this model."""

    __slots__ = ()

    def canonical_id(self) -> str:
        return "H2O"

    def formula(self) -> str:
        return "H2O"

    @property
    def n_carbon(self) -> int:
        return 0

    @property
    def is_radical(self) -> bool:
        return False

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return isinstance(other, _H2O)

    def __hash__(self):
        return hash("H2O")

    def __repr__(self):
        return "Molecule(H2O)"


ATOMIC_O = _AtomicO()
HYDROXYL_OH = _Hydroxyl()
MOLECULAR_O2 = _O2()
MOLECULAR_O3 = _O3()
MOLECULAR_N2 = _N2()
MOLECULAR_CO2 = _CO2()
MOLECULAR_CO = _CO()
ATOMIC_AR = _Argon()
MOLECULAR_H2O = _H2O()


class _HydrocarbonOxideWrapper:
    """Shared machinery for the four "hydrocarbon fragment + O-group"
    species (peroxy/alkoxy radicals, hydroperoxides, alcohols): each just
    tags a hydrocarbon Molecule fragment with a fixed O-containing group and
    otherwise delegates to it. No further oxidation chemistry is modeled
    past this first step (radicalizable_sites is always empty), which caps
    this model's oxidation chemistry at one O-addition per radical."""

    __slots__ = ("parent",)
    _tag: str = ""
    _suffix: str = ""
    _is_radical: bool = False

    def __init__(self, parent: "Molecule"):
        self.parent = parent

    def canonical_id(self) -> str:
        return f"{self._tag}:{self.parent.canonical_id()}"

    def formula(self) -> str:
        return f"C{self.parent.n_carbon}H{self.parent.total_h}{self._suffix}"

    @property
    def n_carbon(self) -> int:
        return self.parent.n_carbon

    @property
    def is_radical(self) -> bool:
        return self._is_radical

    def radicalizable_sites(self) -> list[Site]:
        return []

    def __eq__(self, other):
        return type(other) is type(self) and self.canonical_id() == other.canonical_id()

    def __hash__(self):
        return hash(self.canonical_id())

    def __repr__(self):
        return f"Molecule({self.formula()})"


class PeroxyRadical(_HydrocarbonOxideWrapper):
    """R* + O2 -> ROO*"""
    _tag, _suffix, _is_radical = "ROO", "OO•", True


class AlkoxyRadical(_HydrocarbonOxideWrapper):
    """R* + O* -> RO*"""
    _tag, _suffix, _is_radical = "RO", "O•", True


class Hydroperoxide(_HydrocarbonOxideWrapper):
    """ROO* + H* -> ROOH"""
    _tag, _suffix, _is_radical = "ROOH", "OOH", False


class Alcohol(_HydrocarbonOxideWrapper):
    """RO* + H* (or R* + OH*) -> ROH"""
    _tag, _suffix, _is_radical = "ROH", "OH", False


def _pair(a, b, ta, tb) -> bool:
    return isinstance(a, ta) and isinstance(b, tb)


def combine(a: Species, b: Species) -> Species | None:
    """Form the product of a radical-radical combination a + b -> product.

    Returns None for radical pairs this model deliberately doesn't track
    further (e.g. two peroxy radicals meeting) -- callers must skip
    generating a reaction in that case rather than treat it as an error.
    """
    a_is_h, b_is_h = isinstance(a, _AtomicH), isinstance(b, _AtomicH)
    if a_is_h and b_is_h:
        return MOLECULAR_H2
    if a_is_h and isinstance(b, Molecule):
        return b.saturate_radical()
    if isinstance(a, Molecule) and b_is_h:
        return a.saturate_radical()
    if isinstance(a, Molecule) and isinstance(b, Molecule):
        return Molecule.join(a, b)

    if _pair(a, b, _AtomicO, _AtomicO):
        return MOLECULAR_O2
    if _pair(a, b, _AtomicO, _AtomicH) or _pair(b, a, _AtomicO, _AtomicH):
        return HYDROXYL_OH
    if _pair(a, b, _Hydroxyl, _AtomicH) or _pair(b, a, _Hydroxyl, _AtomicH):
        return MOLECULAR_H2O
    if _pair(a, b, _AtomicO, Molecule) or _pair(b, a, _AtomicO, Molecule):
        radical = a if isinstance(a, Molecule) else b
        return AlkoxyRadical(radical)
    if _pair(a, b, _Hydroxyl, Molecule) or _pair(b, a, _Hydroxyl, Molecule):
        radical = a if isinstance(a, Molecule) else b
        return Alcohol(radical)
    if _pair(a, b, PeroxyRadical, _AtomicH) or _pair(b, a, PeroxyRadical, _AtomicH):
        wrapper = a if isinstance(a, PeroxyRadical) else b
        return Hydroperoxide(wrapper.parent)
    if _pair(a, b, AlkoxyRadical, _AtomicH) or _pair(b, a, AlkoxyRadical, _AtomicH):
        wrapper = a if isinstance(a, AlkoxyRadical) else b
        return Alcohol(wrapper.parent)

    return None  # not modeled: e.g. two peroxy/alkoxy radicals meeting, OH+O, etc.


# -- seed molecule factories, for choosing the simulation's starting species --

def methane() -> Molecule:
    g = nx.Graph()
    g.add_node(0, radical=False)
    return Molecule(g)


def ethane() -> Molecule:
    g = nx.Graph()
    g.add_node(0, radical=False)
    g.add_node(1, radical=False)
    g.add_edge(0, 1)
    return Molecule(g)


def propane() -> Molecule:
    g = nx.Graph()
    g.add_nodes_from([(0, {"radical": False}), (1, {"radical": False}), (2, {"radical": False})])
    g.add_edge(0, 1)
    g.add_edge(1, 2)
    return Molecule(g)


SEED_MOLECULES = {
    "Methane (CH4)": methane,
    "Ethane (C2H6)": ethane,
    "Propane (C3H8)": propane,
}

# Fixed singleton atmospheric gases, exposed as zero-arg factories for the
# same seeding interface as SEED_MOLECULES (Simulator.seed_species expects
# a fresh-looking Species instance; these are stateless so returning the
# shared singleton each time is safe).
ATMOSPHERIC_GASES = {
    "N2": lambda: MOLECULAR_N2,
    "O2": lambda: MOLECULAR_O2,
    "CO2": lambda: MOLECULAR_CO2,
    "Ar": lambda: ATOMIC_AR,
}
