import networkx as nx
import numpy as np
from tmc_tools.constants import atomic_numbers, atomic_masses


def molecular_graph_determinant(graph, return_string=True):
    weights = np.array(
        [
            atomic_masses[atomic_numbers[atom[1]]]
            for atom in graph.nodes(data="symbol", default="X")
        ]
    )
    adjacency = nx.adjacency_matrix(graph)
    mat = np.outer(weights, weights) * adjacency.toarray()
    np.fill_diagonal(mat, weights)
    det = np.linalg.det(mat)
    if return_string:
        det = str(det)
        if "e+" in det:
            sp = det.split("e+")
            det = sp[0][0:12] + "e+" + sp[1]
        else:
            det = det[0:12]
    return det
