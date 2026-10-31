from Bio.SeqUtils import seq1
from Bio.PDB import is_aa

def extract_sequence(structure):
    sequences={}
    for model in structure:
        for chain in model:
            sequence=""
            for residue in chain:
                if is_aa(residue, standard=True):
                    amino_acid=seq1(residue.resname)
                    sequence= sequence + amino_acid
        sequences[chain.id]=sequence
    return sequences
                