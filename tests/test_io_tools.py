from tmc_tools.io_tools import fix_terachem_molden


def test_fix_terachem_molden(resource_path_root):

    patched_lines = fix_terachem_molden(resource_path_root / 'crf6.molden')

    with open(resource_path_root / 'crf6_fixed.molden', 'r') as fin:
        ref_lines = fin.readlines()

    assert patched_lines == ref_lines

    # Assert that patched files do not get patched again

    patched_lines = fix_terachem_molden(
        resource_path_root / 'crf6_fixed.molden')
    assert patched_lines == ref_lines
