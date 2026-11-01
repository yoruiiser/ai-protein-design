import numpy as np


def calculate_rmsd(coords1, coords2):

    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have equal shapes.")

    if coords1.ndim != 2 or coords1.shape[1] != 3:
        raise ValueError("Coordinates must have shape (N, 3).")

    if len(coords1) == 0:
        raise ValueError("Coordinate arrays cannot be empty.")

    difference = coords1 - coords2

    squared_distances = np.sum(difference ** 2, axis=1)

    rmsd = np.sqrt(np.mean(squared_distances))

    return rmsd

def kabsch_alignment(reference, mobile):

    if reference.shape != mobile.shape:
        raise ValueError("Coordinate arrays must have equal shapes.")

    if reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("Coordinates must have shape (N, 3).")

    if len(reference) < 3:
        raise ValueError("At least three points are required.")

    centroid_reference = np.mean(reference, axis=0)
    centroid_mobile = np.mean(mobile, axis=0)


    centered_reference = reference - centroid_reference
    centered_mobile = mobile - centroid_mobile

    covariance = centered_mobile.T @ centered_reference

    U, S, Vt = np.linalg.svd(covariance)

    correction = np.eye(3)

    if np.linalg.det(U @ Vt) < 0:
        correction[2, 2] = -1

    rotation = U @ correction @ Vt
    aligned_mobile = (
        centered_mobile @ rotation + centroid_reference
    )

    rmsd = calculate_rmsd(reference, aligned_mobile)

    return aligned_mobile, rotation, rmsd