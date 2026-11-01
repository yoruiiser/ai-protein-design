from collections import Counter


def calculate_secondary_structure_composition(results):

    # Extract secondary structure codes
    structures = [
        residue["secondary_structure"]
        for residue in results
    ]

    counts = Counter(structures)

    # Group DSSP assignments
    helix = sum(counts[code] for code in "HGI")
    sheet = sum(counts[code] for code in "EB")
    coil = sum(counts[code] for code in "TS-")

    total = len(structures)

    if total == 0:
        raise ValueError("No secondary structure assignments found.")

    composition = {
        "Helix": {
            "count": helix,
            "percentage": helix / total * 100
        },
        "Sheet": {
            "count": sheet,
            "percentage": sheet / total * 100
        },
        "Coil": {
            "count": coil,
            "percentage": coil / total * 100
        }
    }

    return composition