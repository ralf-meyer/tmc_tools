import tmc_tools.constants as constants


def d_orbital_occupation(metal, oxidation, multiplicity):
    """Calculates total number of electrons and number of unpaired electrons
    in the t2g and eg subspaces.

    Parameters
    ----------
    metal : str
        _description_
    oxidation : int | str
        _description_
    multiplicity : int
        _description_

    Returns
    -------
    t2g_total : int
        Total number of electrons in the t2g subspace
    eg_total : int
        Total number of electrons in the eg subspace
    t2g_unpaired : int
        Number of unpaired electrons in the t2g subspace
    eg_unpaired : int
        Number of unpaired electrons in the eg subspace
    """
    if isinstance(oxidation, str):
        oxidation = constants.roman_numerals[oxidation]
    # Check that metal is in the 3d block
    assert 21 <= constants.atomic_numbers[metal] <= 30
    # Subtract 18 electrons (Argon configuration) and oxidation state.
    n_valence = constants.atomic_numbers[metal] - 18 - oxidation
    n_unpaired = multiplicity - 1
    n_paired = n_valence - n_unpaired
    assert n_paired % 2 == 0

    # Total number of spin up electrons (half of the pairs + unpaired)
    n_alpha = n_paired // 2 + n_unpaired
    # Distinguish 3 possible cases:
    # eg empty
    if n_alpha <= 3:
        return n_valence, 0, n_unpaired, 0
    # t2g partially filled
    if n_paired < 6:
        t2g_unpaired = 3 - n_paired // 2
        return (n_paired + t2g_unpaired,
                n_unpaired - t2g_unpaired,
                t2g_unpaired,
                n_unpaired - t2g_unpaired)
    # t2g fully filled
    return 6, n_valence - 6, 0, n_unpaired


def dbloc_coefficients(metal, oxidation, gs_mult, es_mult):

    (gs_t2g_total, gs_eg_total, gs_t2g_unpaired,
     gs_eg_unpaired) = d_orbital_occupation(metal, oxidation, gs_mult)
    (es_t2g_total, es_eg_total, es_t2g_unpaired,
     es_eg_unpaired) = d_orbital_occupation(metal, oxidation, es_mult)

    p = 0.5 * (gs_t2g_unpaired - es_t2g_unpaired
               + gs_eg_unpaired - es_eg_unpaired)
    ss = 0.5 * (es_t2g_unpaired * (es_t2g_unpaired - 1)
                + es_eg_unpaired * (es_eg_unpaired - 1)
                - gs_t2g_unpaired * (gs_t2g_unpaired - 1)
                - gs_eg_unpaired * (gs_eg_unpaired - 1))
    return p, ss
