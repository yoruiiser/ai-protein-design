from Bio.PDB import PDBParser, DSSP
def assign_secondary_structure(pdb_file):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]

    dssp = DSSP(
        model,
        pdb_file,
        dssp="mkdssp"
    )

    results = []

    for key in dssp.keys():

        residue_data = dssp[key]

        chain_id = key[0]
        residue_id = key[1][1]

        amino_acid = residue_data[1]
        secondary_structure = residue_data[2]
        relative_accessibility = residue_data[3]

        results.append({
            "chain": chain_id,
            "residue_number": residue_id,
            "amino_acid": amino_acid,
            "secondary_structure": secondary_structure,
            "relative_accessibility": relative_accessibility
        })

    return results