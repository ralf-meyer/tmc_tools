from tmc_tools.orbital_energies import orbital_energies


def test_orbital_diagram(resource_path_root):

    energ_a, energ_b = orbital_energies(
        resource_path_root / 'crf6_fixed.molden')

    assert len(energ_a) == 5
    assert len(energ_b) == 5
