from tmc_tools.orbital_energies import orbital_energies


def test_orbital_diagram(resource_path_root):

    energies = orbital_energies(resource_path_root / 'crf6_fixed.molden')

    assert len(energies) == 5
