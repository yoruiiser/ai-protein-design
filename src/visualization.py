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