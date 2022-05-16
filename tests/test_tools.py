from tmc_tools.tools import fix_terachem_molden
from pkg_resources import resource_filename, Requirement


def test_fix_terachem_molden():
    prepatch = resource_filename(
        Requirement.parse('tmc_tools'),
        'tests/resources/crf6.molden')

    patched_lines = fix_terachem_molden(prepatch)

    postpatch = resource_filename(
        Requirement.parse('tmc_tools'),
        'tests/resources/crf6_fixed.molden')

    with open(postpatch, 'r') as fin:
        ref_lines = fin.readlines()

    assert patched_lines == ref_lines

    # Assert that patched files do not get patched again

    patched_lines = fix_terachem_molden(postpatch)
    assert patched_lines == ref_lines
