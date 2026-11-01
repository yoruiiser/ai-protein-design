import numpy as np
from Bio.PDB import is_aa
def extract_ca_coordinates(structure):

    coordinates = []
    residue_ids = []

    for model in structure:

        for chain in model:

            for residue in chain:

                if is_aa(residue, standard=True):

                    if residue.has_id("CA"):

                        ca_atom = residue["CA"]

                        coordinates.append(ca_atom.coord)
                        residue_ids.append(
                            (chain.id, residue.id[1])
                        )

    coordinates = np.array(coordinates)

    return coordinates, residue_ids


def calculate_distance_matrix(coordinates):

    difference = (
        coordinates[:, np.newaxis, :]
        - coordinates[np.newaxis, :, :]
    )

    squared_difference = difference ** 2

    distance_matrix = np.sqrt(
        np.sum(squared_difference, axis=2)
    )

    return distance_matrix

def calculate_contact_map(distance_matrix, cutoff=8.0):

    contact_map = distance_matrix < cutoff

    return contact_map.astype(int)