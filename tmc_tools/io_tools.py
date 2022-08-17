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
