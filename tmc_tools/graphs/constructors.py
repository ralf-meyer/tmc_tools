import numpy as np
import networkx as nx
from tmc_tools.constants import covalent_radii


def graph_from_atoms(atoms, threshold=1.2):
    g = nx.Graph()
    for i, atom in enumerate(atoms):
        g.add_node(i, symbol=atom.symbol)

    for i, ai in enumerate(atoms):
        for j, aj in enumerate(atoms[i+1:]):
            r = np.linalg.norm(ai.position - aj.position)
            if (r < threshold * (covalent_radii[ai.number]
                                 + covalent_radii[aj.number])):
                g.add_edge(i, j+i+1)
    return g
