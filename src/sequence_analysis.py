from Bio.SeqUtils import seq1
from Bio.PDB import is_aa
from collections import Counter

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


def calculate_composition(sequence):
    counts=Counter(sequence)
    return counts
def calculate_frequency(sequence):
    counts=calculate_composition(sequence)
    length=len(sequence)
    if length==0:
        return {}
    frequencies={amino_acid: count/length 
                 for amino_acid, count in counts.items()}
    return frequencies
def calculate_statistics(sequence):
    return{
        "length": len(sequence),
        "unique_amino_acids": len(set(sequence)),
        "glycine_count": sequence.count('G'),
        "proline_count": sequence.count('P'),
        "aromatic_count": sum(sequence.count(aa) for aa in "FWY")}
   
    