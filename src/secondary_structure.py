from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from collections import Counter


def assign_secondary_structure(file_path):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        file_path
    )

    model = structure[0]

    dssp = DSSP(
        model,
        file_path,
        dssp="mkdssp"
    )

    assignments = {}

    for key in dssp.keys():

        residue_data = dssp[key]

        amino_acid = residue_data[1]
        secondary_structure = residue_data[2]

        assignments[key] = {
            "amino_acid": amino_acid,
            "secondary_structure": secondary_structure
        }

    return assignments

def calculate_secondary_composition(assignments):

    counts = Counter()

    for residue in assignments.values():

        code = residue["secondary_structure"]

        if code in ["H", "G", "I"]:
            counts["Helix"] += 1

        elif code in ["E", "B"]:
            counts["Sheet"] += 1

        else:
            counts["Other"] += 1

    total = sum(counts.values())

    percentages = {
        category: (count / total) * 100
        for category, count in counts.items()
    }

    return counts, percentages