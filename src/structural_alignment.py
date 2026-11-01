from Bio.PDB import PDBParser, is_aa
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



def extract_ca_by_residue(file_path, chain_id="A"):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)

    model = structure[0]
    chain = model[chain_id]

    residues = {}

    for residue in chain:

        if not is_aa(residue, standard=True):
            continue

        if not residue.has_id("CA"):
            continue

        residue_key = (
            residue.id[1],
            residue.id[2]
        )

        residues[residue_key] = {
            "name": residue.resname,
            "coordinate": residue["CA"].coord
        }

    return residues



def match_residues(reference, mobile):

    common_keys = sorted(
        set(reference.keys()) & set(mobile.keys())
    )

    matched_reference = []
    matched_mobile = []

    for key in common_keys:

        ref_residue = reference[key]
        mob_residue = mobile[key]

        if ref_residue["name"] != mob_residue["name"]:
            continue

        matched_reference.append(
            ref_residue["coordinate"]
        )

        matched_mobile.append(
            mob_residue["coordinate"]
        )

    return (
        np.array(matched_reference),
        np.array(matched_mobile)
    )