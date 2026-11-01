import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np

from structural_alignment import (
    calculate_rmsd,
    kabsch_alignment
)



reference = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 0.0, 3.0],
    [1.0, 2.0, 3.0]
])



angle = np.radians(45)

rotation = np.array([
    [np.cos(angle), -np.sin(angle), 0],
    [np.sin(angle), np.cos(angle), 0],
    [0, 0, 1]
])

mobile = reference @ rotation


mobile += np.array([5.0, 8.0, 3.0])



raw_rmsd = calculate_rmsd(reference, mobile)


aligned, calculated_rotation, aligned_rmsd = (
    kabsch_alignment(reference, mobile)
)


print("Raw RMSD:", raw_rmsd)
print("Aligned RMSD:", aligned_rmsd)

print("\nCalculated rotation matrix:")
print(calculated_rotation)

print("\nAligned coordinates:")
print(aligned)