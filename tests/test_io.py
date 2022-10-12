import numpy as np
from tmc_tools.io import read_mol, fix_terachem_molden


def test_read_mol(resource_path_root):
    with open(resource_path_root / "cyanide.mol", "r") as fin:
        atoms = read_mol(fin)
    np.testing.assert_equal(atoms.get_initial_charges(), [-1.0, 0.0])

    with open(resource_path_root / "oxalate.mol", "r") as fin:
        atoms = read_mol(fin)
    np.testing.assert_equal(
        atoms.get_initial_charges(), [0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
    )


def test_fix_terachem_molden(resource_path_root):

    patched_lines = fix_terachem_molden(resource_path_root / "crf6.molden")

    with open(resource_path_root / "crf6_fixed.molden", "r") as fin:
        ref_lines = fin.readlines()

    assert patched_lines == ref_lines

    # Assert that patched files do not get patched again

    patched_lines = fix_terachem_molden(resource_path_root / "crf6_fixed.molden")
    assert patched_lines == ref_lines
