from secondary_structure import (
    assign_secondary_structure,
    calculate_secondary_composition
)

from visualization import plot_secondary_structure


def main():

    assignments = assign_secondary_structure(
        "data/raw/1LYZ.pdb"
    )

    counts, percentages = calculate_secondary_composition(
        assignments
    )

    print("\nSECONDARY STRUCTURE COUNTS")

    for category, count in counts.items():
        print(category, ":", count)

    print("\nSECONDARY STRUCTURE PERCENTAGES")

    for category, percentage in percentages.items():
        print(category, ":", round(percentage, 2), "%")

    plot_secondary_structure(
        percentages,
        "results/figures/secondary_structure.png"
    )


if __name__ == "__main__":
    main()