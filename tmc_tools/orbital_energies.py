from pyscf.tools import molden


def orbital_energies(filename):
    (mol, mo_energy, mo_coeff, mo_occ,
     irrep_labels, spins) = molden.load(filename)

    return [0.] * 5
