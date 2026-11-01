from pdb_parser import parse_pdb
from sequence_analysis import (
    extract_sequence,
    calculate_composition,
    calculate_frequency,
    calculate_statistics
)
from visualization import plot_composition


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


if __name__ == "__main__":
    main()