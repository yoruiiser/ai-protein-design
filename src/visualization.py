import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

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



def plot_ramachandran(results, output_path):

    phi_values = [
        item["phi"] for item in results
        if item["phi"] is not None and item["psi"] is not None
    ]

    psi_values = [
        item["psi"] for item in results
        if item["phi"] is not None and item["psi"] is not None
    ]

    plt.figure(figsize=(8, 8))

    plt.scatter(phi_values, psi_values, alpha=0.7)

    plt.xlim(-180, 180)
    plt.ylim(-180, 180)

    plt.xlabel("Phi (degrees)")
    plt.ylabel("Psi (degrees)")

    plt.title("Ramachandran Plot")

    plt.grid(alpha=0.2)

    plt.tight_layout()

    plt.savefig(output_path, dpi=300)

    plt.show()



def plot_secondary_structure(percentages, output_path):
    labels = list(percentages.keys())
    values = list(percentages.values())

    plt.figure(figsize=(7, 5))

    plt.bar(labels, values)

    plt.xlabel("Secondary Structure")
    plt.ylabel("Percentage (%)")

    plt.title("Protein Secondary Structure Composition")

    plt.tight_layout()

    plt.savefig(output_path, dpi=300)

    plt.show()

def plot_secondary_structure(composition, output_path):

    categories = list(composition.keys())

    percentages = [
        composition[category]["percentage"]
        for category in categories
    ]

    counts = [
        composition[category]["count"]
        for category in categories
    ]

    plt.figure(figsize=(8, 5))

    bars = plt.bar(categories, percentages)

    plt.xlabel("Secondary Structure")
    plt.ylabel("Percentage of Residues")
    plt.title("Secondary Structure Composition of Lysozyme (1LYZ)")

    plt.ylim(0, 100)
    
    for bar, percentage, count in zip(bars, percentages, counts):

        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1,
            f"{percentage:.1f}% ({count})",
            ha="center"
        )

    plt.tight_layout()

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")

    plt.show()

    plt.close()