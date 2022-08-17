import pytest
import ase.io
import networkx as nx
from tmc_tools.graphs.constructors import graph_from_ase_atoms, graph_from_mol_file


def test_water(resource_path_root):
    g_ref = nx.Graph()
    g_ref.add_nodes_from(
        [(0, {"symbol": "O"}), (1, {"symbol": "H"}), (2, {"symbol": "H"})]
    )
    g_ref.add_edges_from([(0, 1), (0, 2)])

    atoms = ase.io.read(resource_path_root / "water.mol")
    g = graph_from_ase_atoms(atoms)
    assert g.nodes == g_ref.nodes
    assert g.edges == g_ref.edges

    g = graph_from_mol_file(resource_path_root / "water.mol")
    assert g.nodes == g_ref.nodes
    assert g.edges == g_ref.edges


@pytest.fixture
def furan_graph():
    g = nx.Graph()
    g.add_nodes_from(
        [
            (0, {"symbol": "O"}),
            (1, {"symbol": "C"}),
            (2, {"symbol": "C"}),
            (3, {"symbol": "C"}),
            (4, {"symbol": "C"}),
            (5, {"symbol": "H"}),
            (6, {"symbol": "H"}),
            (7, {"symbol": "H"}),
            (8, {"symbol": "H"}),
        ]
    )
    g.add_edges_from(
        [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (1, 5), (2, 6), (4, 7), (3, 8)]
    )
    return g


def test_furan(resource_path_root, furan_graph):

    atoms = ase.io.read(resource_path_root / "furan.mol")
    g = graph_from_ase_atoms(atoms)
    assert g.nodes == furan_graph.nodes
    assert g.edges == furan_graph.edges

    g = graph_from_mol_file(resource_path_root / "furan.mol")
    assert g.nodes == furan_graph.nodes
    assert g.edges == furan_graph.edges
