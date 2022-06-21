import numpy as np
roman_numerals = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}

# Copied from the atomic simulation environment
chemical_symbols = [
    # 0
    'X',
    # 1
    'H', 'He',
    # 2
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    # 3
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    # 4
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    # 5
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    # 6
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    # 7
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']

atomic_numbers = {}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z


# Covalent radii from ase:
#
#  Covalent radii revisited,
#  Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
#  Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez,
#  Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J
missing = 0.2
covalent_radii = np.array([
    missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    missing,  # Bk
    missing,  # Cf
    missing,  # Es
    missing,  # Fm
    missing,  # Md
    missing,  # No
    missing,  # Lr
    missing,  # Rf
    missing,  # Db
    missing,  # Sg
    missing,  # Bh
    missing,  # Hs
    missing,  # Mt
    missing,  # Ds
    missing,  # Rg
    missing,  # Cn
    missing,  # Nh
    missing,  # Fl
    missing,  # Mc
    missing,  # Lv
    missing,  # Ts
    missing,  # Og
])


# From:
# https://en.wikipedia.org/wiki/Electronegativities_of_the_elements_(data_page)
no_data = np.NaN
electronegativity = np.array([
    2.20,  # H
    no_data,  # He
    0.98,  # Li
    1.57,  # Be
    2.04,  # B
    2.55,  # C
    3.04,  # N
    3.44,  # O
    3.98,  # F
    no_data,  # Ne
    0.93,  # Na
    1.31,  # Mg
    1.61,  # Al
    1.90,  # Si
    2.19,  # P
    2.58,  # S
    3.16,  # Cl
    no_data,  # Ar
    0.82,  # K
    1.00,  # Ca
    1.36,  # Sc
    1.54,  # Ti
    1.63,  # V
    1.66,  # Cr
    1.55,  # Mn
    1.83,  # Fe
    1.88,  # Co
    1.91,  # Ni
    1.90,  # Cu
    1.65,  # Zn
    1.81,  # Ga
    2.01,  # Ge
    2.18,  # As
    2.55,  # Se
    2.96,  # Br
    3.00,  # Kr
    0.82,  # Rb
    0.95,  # Sr
    1.22,  # Y
    1.33,  # Zr
    1.60,  # Nb
    2.16,  # Mo
    1.90,  # Tc
    2.20,  # Ru
    2.28,  # Rh
    2.20,  # Pd
    1.93,  # Ag
    1.69,  # Cd
    1.78,  # In
    1.96,  # Sn
    2.05,  # Sb
    2.10,  # Te
    2.66,  # I
    2.60,  # Xe
    0.79,  # Cs
    0.89,  # Ba
    1.10,  # La
    1.12,  # Ce
    1.13,  # Pr
    1.14,  # Nd
    no_data,  # Pm
    1.17,  # Sm
    no_data,  # Eu
    1.20,  # Gd
    no_data,  # Tb
    1.22,  # Dy
    1.23,  # Ho
    1.24,  # Er
    1.25,  # Tm
    no_data,  # Yb
    1.27,  # Lu
    1.30,  # Hf
    1.50,  # Ta
    2.36,  # W
    1.90,  # Re
    2.20,  # Os
    2.20,  # Ir
    2.28,  # Pt
    2.54,  # Au
    2.00,  # Hg
    1.62,  # Tl
    2.33,  # Pb
    2.02,  # Bi
    2.00,  # Po
    2.20,  # At
    no_data,  # Rn
    no_data,  # Fr
    0.90,  # Ra
    1.10,  # Ac
    1.30,  # Th
    1.50,  # Pa
    1.38,  # U
    1.36,  # Np
    1.28,  # Pu
    1.30,  # Am
    1.30,  # Cm
    1.30,  # Bk
    1.30,  # Cf
    1.30,  # Es
    1.30,  # Fm
    1.30,  # Md
    1.30   # No
])
