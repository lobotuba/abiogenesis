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


def combine(a: Species, b: Species) -> Species:
    """Form the product of a radical-radical combination a + b -> product."""
    a_is_h, b_is_h = isinstance(a, _AtomicH), isinstance(b, _AtomicH)
    if a_is_h and b_is_h:
        return MOLECULAR_H2
    if a_is_h and isinstance(b, Molecule):
        return b.saturate_radical()
    if isinstance(a, Molecule) and b_is_h:
        return a.saturate_radical()
    if isinstance(a, Molecule) and isinstance(b, Molecule):
        return Molecule.join(a, b)
    raise ValueError(f"cannot combine {a!r} and {b!r}")


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
