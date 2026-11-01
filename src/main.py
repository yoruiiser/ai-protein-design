from pdb_parser import parse_pdb
from sequence_analysis import (
    extract_sequence,
    calculate_composition,
    calculate_frequency,
    calculate_statistics
)
from visualization import plot_composition
from visualization import plot_distance_matrix
from visualization import plot_secondary_structure
from structure_analysis import (
    extract_ca_coordinates,
    calculate_distance_matrix,
    calculate_contact_map
)
from backbone_analysis import extract_backbone_angles
from visualization import plot_ramachandran
from secondary_structure import assign_secondary_structure
from secondary_structure_analysis import calculate_secondary_structure_composition


def main():

    structure = parse_pdb("data/raw/1LYZ.pdb")

    sequences = extract_sequence(structure)

    for chain_id, sequence in sequences.items():

        print("\nCHAIN:", chain_id)
        print("SEQUENCE:", sequence)

        composition = calculate_composition(sequence)
        plot_composition(composition,f"results/figures/chain_{chain_id}_composition.png")
        frequencies = calculate_frequency(sequence)
        statistics = calculate_statistics(sequence)
        

        print("\nAMINO ACID COMPOSITION:")

        for amino_acid, count in sorted(composition.items()):
            print(amino_acid, ":", count)

        print("\nAMINO ACID FREQUENCIES:")

        for amino_acid, frequency in sorted(frequencies.items()):
            print(amino_acid, ":", round(frequency, 4))

        print("\nSEQUENCE STATISTICS:")

        for key, value in statistics.items():
            print(key, ":", value)

    
        coordinates, residue_ids = extract_ca_coordinates(structure)

        distance_matrix = calculate_distance_matrix(coordinates)

        contact_map = calculate_contact_map(distance_matrix)
        plot_distance_matrix(
    distance_matrix,
    "results/figures/distance_matrix.png")

        print("\nSTRUCTURAL ANALYSIS")

        print("Coordinate shape:", coordinates.shape)

        print("Distance matrix shape:", distance_matrix.shape)

        print("Contact map shape:", contact_map.shape)

        print("Distance between first two CA atoms:",
                distance_matrix[0, 1])

    
    angles = extract_backbone_angles("data/raw/1LYZ.pdb")
    for item in angles[:10]:
        print(
            item["amino_acid"],
            "Phi:", item["phi"],
            "Psi:", item["psi"]
        )

    plot_ramachandran(
        angles,
        "results/figures/ramachandran.png"
    )
    results = assign_secondary_structure("data/raw/1LYZ.pdb")

    print("\nSECONDARY STRUCTURE RESULTS")
    composition = calculate_secondary_structure_composition(results)

    print("\nSECONDARY STRUCTURE COMPOSITION")

    for category, values in composition.items():
        print(
        f"{category}: {values['count']} residues "
        f"({values['percentage']:.2f}%)"
    )

    for residue in results[:10]:
        print(residue)
    plot_secondary_structure(
    composition,
    "results/figures/secondary_structure.png"
)


if __name__ == "__main__":
 main()