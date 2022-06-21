import numpy as np
import networkx as nx
from tmc_tools.constants import (atomic_numbers, electronegativity,
                                 covalent_radii)


def property_vector(graph, node):
    output = np.zeros(5)
    symbol = graph.nodes[node]['symbol']
    Z = atomic_numbers[symbol]
    # property (i): nuclear charge Z
    output[0] = Z
    # property (ii): Pauling electronegativity chi
    output[1] = electronegativity[Z]
    # property (iii): topology T, coordination number
    output[2] = len(list(graph.neighbors(node)))
    # property (iv): identity
    output[3] = 1.
    # property (v): covalent radius S
    output[4] = covalent_radii[Z]
    return output


def atom_centered_AC(graph, starting_node, depth: int = 3):
    output = np.zeros((depth + 1, 5))
    # Generate all paths from the starting node to all possible nodes
    paths = nx.shortest_path(graph, source=starting_node)
    p_i = property_vector(graph, starting_node)
    for node, path in paths.items():
        # Depth for this atom is given by the length of the path
        d_ij = len(path) - 1
        if d_ij <= depth:
            p_j = property_vector(graph, node)
            output[d_ij] += p_i * p_j
    return output.flatten()
