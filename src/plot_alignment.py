import matplotlib.pyplot as plt


def plot_superposition(reference, mobile, output_path):

    fig = plt.figure(figsize=(9, 7))

    ax = fig.add_subplot(111, projection="3d")

    ax.plot(
        reference[:, 0],
        reference[:, 1],
        reference[:, 2],
        label="Reference",
        linewidth=2
    )

    ax.plot(
        mobile[:, 0],
        mobile[:, 1],
        mobile[:, 2],
        label="Aligned structure",
        linewidth=2
    )

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")

    ax.set_title("Protein Structural Superposition")

    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()