import numpy as np
import networkx as nx
from tmc_tools.constants import covalent_radii
from tmc_tools.utils import deprecated


def graph_from_ase_atoms(atoms, threshold=1.2, covalent_radii=covalent_radii):
    g = nx.Graph()
    for i, atom in enumerate(atoms):
        g.add_node(i, symbol=atom.symbol)

    for i, ai in enumerate(atoms):
        for j, aj in enumerate(atoms[i + 1 :]):
            r = np.linalg.norm(ai.position - aj.position)
            if r < threshold * (covalent_radii[ai.number] + covalent_radii[aj.number]):
                g.add_edge(i, j + i + 1)
    return g


@deprecated("Use graph_from_ase_atoms() instead")
def graph_from_atoms(atoms, **kwargs):
    return graph_from_ase_atoms(atoms, **kwargs)


def graph_from_xyz_file(file, **kwargs):
    import ase.io

    atoms = ase.io.read(file)
    return graph_from_ase_atoms(atoms, **kwargs)


def graph_from_mol_file(file):
    with open(file, "r") as fin:
        lines = fin.readlines()

    # Read counts line:
    sp = lines[3].split()
    n_atoms = int(sp[0])
    n_bonds = int(sp[1])

    g = nx.Graph()

    # Add atoms (offset of 4 for the header lines):
    for i, line in enumerate(lines[4 : 4 + n_atoms]):
        sp = line.split()
        g.add_node(i, symbol=sp[3])

    # Add bonds:
    for line in lines[4 + n_atoms : 4 + n_atoms + n_bonds]:
        sp = line.split()
        # Subtract 1 because of zero indexing vs. one indexing
        g.add_edge(int(sp[0]) - 1, int(sp[1]) - 1)
    return g


def graph_from_mol2_file(file, positions=False):
    with open(file, "r") as fin:
        lines = fin.readlines()

    # Read counts line:
    sp = lines[2].split()
    n_atoms = int(sp[0])
    n_bonds = int(sp[1])

    g = nx.Graph()

    atom_start = lines.index("@<TRIPOS>ATOM\n") + 1
    for i, line in enumerate(lines[atom_start : atom_start + n_atoms]):
        sp = line.split()
        sym = sp[5].split(".")[0]
        data = dict(symbol=sym)
        if positions:
            pos = [float(p) for p in sp[2:5]]
            data.update(position=pos)
        g.add_node(i, **data)

    bond_start = lines.index("@<TRIPOS>BOND\n") + 1
    for line in lines[bond_start : bond_start + n_bonds]:
        sp = line.split()
        # Subtract 1 because of zero indexing vs. one indexing
        g.add_edge(int(sp[1]) - 1, int(sp[2]) - 1)

    return g
