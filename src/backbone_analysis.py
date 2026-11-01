from Bio.PDB import PDBParser, PPBuilder
import numpy as np


def extract_backbone_angles(file_path):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)

    polypeptides = PPBuilder().build_peptides(structure)

    results = []

    for peptide in polypeptides:

        sequence = str(peptide.get_sequence())

        phi_psi_angles = peptide.get_phi_psi_list()

        for index, (phi, psi) in enumerate(phi_psi_angles):

            phi_degrees = (
                np.degrees(phi) if phi is not None else None
            )

            psi_degrees = (
                np.degrees(psi) if psi is not None else None
            )

            results.append({
                "residue_index": index + 1,
                "amino_acid": sequence[index],
                "phi": phi_degrees,
                "psi": psi_degrees
            })

    return results