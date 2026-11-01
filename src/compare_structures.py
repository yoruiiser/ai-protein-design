from structural_alignment import (
    extract_ca_by_residue,
    match_residues,
    kabsch_alignment,
    calculate_rmsd
)
from plot_alignment import plot_superposition

def main():

    reference = extract_ca_by_residue(
        "data/raw/1LYZ.pdb"
    )

    mobile = extract_ca_by_residue(
        "data/raw/2LYZ.pdb"
    )

    print("Reference residues:", len(reference))
    print("Mobile residues:", len(mobile))

    coords1, coords2 = match_residues(
        reference,
        mobile
    )

    print("Matched residues:", len(coords1))

    raw_rmsd = calculate_rmsd(coords1, coords2)

    aligned_coords, rotation, aligned_rmsd = (
        kabsch_alignment(coords1, coords2)
    )
    plot_superposition(
    coords1,
    aligned_coords,
    "results/figures/structure_alignment.png"
)

    print("\nSTRUCTURAL COMPARISON")
    print("Raw RMSD:", raw_rmsd)
    print("Aligned RMSD:", aligned_rmsd)

    print("\nRotation matrix:")
    print(rotation)


if __name__ == "__main__":
    main()