from tmc_tools.models.dbloc import d_orbital_occupation, dbloc_coefficients


def test_d_orbital_occuptation():
    """Tests the occupation function on all cases presented in
    Hughes and Friesner J.C.T.C. 7, 19-32, 2011.
    Additional higher or intermediate spin cases are added where
    applicable."""

    # Test both integer / roman numeral call of oxidation and upper and lower
    # case versions of the metal string.
    assert d_orbital_occupation("Fe", "II", 5) == d_orbital_occupation("fe", 2, 5)

    # Extended Hughes and Friesner cases
    # 2 valence electrons
    assert d_orbital_occupation("V", "III", 3) == (2, 0, 2, 0)
    assert d_orbital_occupation("V", "III", 1) == (2, 0, 0, 0)

    # 3 valence electrons
    assert d_orbital_occupation("Mn", "IV", 4) == (3, 0, 3, 0)
    assert d_orbital_occupation("Mn", "IV", 2) == (3, 0, 1, 0)

    # 4 valence electrons
    assert d_orbital_occupation("Mn", "III", 5) == (3, 1, 3, 1)
    assert d_orbital_occupation("Mn", "III", 3) == (4, 0, 2, 0)
    assert d_orbital_occupation("Mn", "III", 1) == (4, 0, 0, 0)

    # 5 valence electrons
    assert d_orbital_occupation("Fe", "III", 6) == (3, 2, 3, 2)
    assert d_orbital_occupation("Fe", "III", 4) == (4, 1, 2, 1)
    assert d_orbital_occupation("Fe", "III", 2) == (5, 0, 1, 0)

    # 6 valence electrons
    assert d_orbital_occupation("Fe", "II", 5) == (4, 2, 2, 2)
    assert d_orbital_occupation("Fe", "II", 3) == (5, 1, 1, 1)
    assert d_orbital_occupation("Fe", "II", 1) == (6, 0, 0, 0)

    # 7 valence electrons
    assert d_orbital_occupation("Co", "II", 4) == (5, 2, 1, 2)
    assert d_orbital_occupation("Co", "II", 2) == (6, 1, 0, 1)

    # 8 valence electrons
    assert d_orbital_occupation("Ni", "II", 3) == (6, 2, 0, 2)
    assert d_orbital_occupation("Ni", "II", 1) == (6, 2, 0, 0)


def test_dbloc_coefficients():

    # 2 valence electrons: V(III) F6
    assert dbloc_coefficients("V", "III", 3, 1) == (1.0, -1.0, 0.0, 0.0, 0.0)
    # 3 valence electrons: Cr(III) 223tetcl2
    assert dbloc_coefficients("Cr", "III", 4, 2) == (1.0, -3.0, 0.0, 0.0, 0.0)
    # 4 valence electrons: Mn(III) CN6
    assert dbloc_coefficients("Mn", "III", 3, 1) == (1.0, -1.0, 0.0, 0.0, 0.0)
    assert dbloc_coefficients("Mn", "III", 5, 1) == (2.0, -3.0, 0.0, 0.0, 0.0)
    # 5 valence electrons: Fe(III) CN6
    # Not same as Table 8 because here lambda=0
    assert dbloc_coefficients("Fe", "III", 6, 4) == (1.0, -3.0, 0.0, 0.0, 0.0)
    assert dbloc_coefficients("Fe", "III", 6, 2) == (2.0, -4.0, 0.0, 0.0, 0.0)
    assert dbloc_coefficients("Fe", "III", 2, 4) == (-1.0, 1.0, 0.0, 0.0, 0.0)
    assert dbloc_coefficients("Fe", "III", 2, 6) == (-2.0, 4.0, 0.0, 0.0, 0.0)
    # 6 valence electrons: Fe(II) CN6
    assert dbloc_coefficients("Fe", "II", 1, 3) == (-1.0, 0.0, 0.0, 0.0, 0.0)
    # 7 valence electrons: Co(II) F6
    # Not same as Table 8 because here lambda=0
    assert dbloc_coefficients("Co", "II", 4, 2) == (1.0, -1.0, 0.0, 0.0, 0.0)
    # 8 valence electrons: Ni(II) F6
    assert dbloc_coefficients("Ni", "II", 3, 1) == (1.0, -1.0, 0.0, 0.0, 0.0)

    # Fe(II) CN6
    assert dbloc_coefficients("Fe", "III", 2, 4, ["carbonyl"] * 6) == (
        -1.0,
        1.0,
        0.0,
        0.0,
        6.0,
    )
