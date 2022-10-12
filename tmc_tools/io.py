from ase.atoms import Atoms


def read_mol(fileobj):
    # Based on the ase function read_mol() and extended to support
    # reading of charges
    lines = fileobj.readlines()
    L1 = lines[3]

    # The V2000 dialect uses a fixed field length of 3, which means there
    # won't be space between the numbers if there are 100+ atoms, and
    # the format doesn't support 1000+ atoms at all.
    if L1.rstrip().endswith("V2000"):
        natoms = int(L1[:3].strip())
    else:
        natoms = int(L1.split()[0])
    positions = []
    symbols = []
    for line in lines[4 : 4 + natoms]:
        x, y, z, symbol = line.split()[:4]
        symbols.append(symbol)
        positions.append([float(x), float(y), float(z)])
    charges = [0] * natoms
    for line in lines[4 + natoms :]:
        if line.startswith("M  CHG"):
            split = line.split()[3:]
            for i in range(0, len(split), 2):
                charges[int(split[i]) - 1] = float(split[i + 1])
            break
    return Atoms(symbols=symbols, positions=positions, charges=charges)


def fix_terachem_molden(moldenfile):
    """The current implementation in terachem is missing the symmetry
    labels for the molecular orbitals. This function iterates through the
    molecular orbitals and adds the missing symmetry labels."""

    with open(moldenfile, "r") as fin:
        lines = fin.readlines()

    output_lines = []
    prev_line = ""
    for line in lines:
        if "ene=" in line.lower() and "sym=" not in prev_line.lower():
            output_lines.append(" Sym= C1\n")
        output_lines.append(line)
        prev_line = line
    return output_lines
