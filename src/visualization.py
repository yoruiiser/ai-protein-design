import matplotlib.pyplot as plt


def plot_composition(composition, output_path):

    amino_acids = list(composition.keys())
    counts = list(composition.values())

    plt.figure(figsize=(10, 5))

    plt.bar(amino_acids, counts)

    plt.xlabel("Amino Acid")
    plt.ylabel("Count")
    plt.title("Protein Amino Acid Composition")

    plt.tight_layout()

    plt.savefig(output_path)

    plt.show()

import matplotlib.pyplot as plt


def plot_distance_matrix(distance_matrix, output_path):

    plt.figure(figsize=(8, 7))

    plt.imshow(distance_matrix, cmap="viridis")

    plt.colorbar(label="Distance (Å)")

    plt.xlabel("Residue index")
    plt.ylabel("Residue index")

    plt.title("Protein Residue Distance Matrix")

    plt.tight_layout()

    plt.savefig(output_path, dpi=300)

    print("Distance matrix saved to:", output_path)

    plt.show()