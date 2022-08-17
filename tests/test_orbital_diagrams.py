from tmc_tools.orbital_energies import get_d_orbital_energies


def test_orbital_diagram(resource_path_root):

    energ_a, energ_b = get_d_orbital_energies(resource_path_root / "crf6_fixed.molden")

    assert len(energ_a) == 5
    assert len(energ_b) == 5
