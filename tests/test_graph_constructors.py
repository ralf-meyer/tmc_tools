import ase.io
import networkx as nx
from tmc_tools.graphs.constructors import graph_from_atoms


def test_water(resource_path_root):
    atoms = ase.io.read(resource_path_root / 'water.mol')
    g = graph_from_atoms(atoms)

    g_ref = nx.Graph()
    g_ref.add_nodes_from([(0, {'symbol': 'O'}),
                          (1, {'symbol': 'H'}),
                          (2, {'symbol': 'H'})])
    g_ref.add_edges_from([(0, 1), (0, 2)])

    assert g.nodes == g_ref.nodes
    assert g.edges == g_ref.edges


def test_furan(resource_path_root):
    atoms = ase.io.read(resource_path_root / 'furan.mol')
    g = graph_from_atoms(atoms)

    g_ref = nx.Graph()
    g_ref.add_nodes_from([(0, {'symbol': 'O'}),
                          (1, {'symbol': 'C'}),
                          (2, {'symbol': 'C'}),
                          (3, {'symbol': 'C'}),
                          (4, {'symbol': 'C'}),
                          (5, {'symbol': 'H'}),
                          (6, {'symbol': 'H'}),
                          (7, {'symbol': 'H'}),
                          (8, {'symbol': 'H'})])
    g_ref.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (1, 5),
                          (2, 6), (4, 7), (3, 8)])

    assert g.nodes == g_ref.nodes
    assert g.edges == g_ref.edges
