from tmc_tools.orbital_energies import orbital_energies
from pkg_resources import resource_filename, Requirement


def test_orbital_diagram(tmpdir):
    molden_file = resource_filename(
        Requirement.parse('tmc_tools'),
        'tests/resources/crf6_fixed.molden')

    energies = orbital_energies(molden_file)

    assert len(energies) == 5
