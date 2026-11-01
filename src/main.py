from pdb_parser import parse_pdb
from sequence_analysis import (
    extract_sequence,
    calculate_composition,
    calculate_frequency,
    calculate_statistics
)
from visualization import plot_composition
from visualization import plot_distance_matrix
from structure_analysis import (
    extract_ca_coordinates,
    calculate_distance_matrix,
    calculate_contact_map
)


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
      


if __name__ == "__main__":
    main()