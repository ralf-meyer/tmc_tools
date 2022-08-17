from pyscf.tools import molden
from pyscf.lo import Boys


def get_d_orbital_energies(filename):
    (mol, mo_energy, mo_coeff, mo_occ, irrep_labels, spins) = molden.load(filename)

    localizer = Boys(mol)
    # First axis AO coefficient, second MO number
    Ca = localizer.kernel(mo_coeff[0])
    return [0.0] * 5, [0.0] * 5
