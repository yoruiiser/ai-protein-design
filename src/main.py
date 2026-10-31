from pdb_parser import parse_pdb
from sequence_analysis import extract_sequence

def main():
    structure=parse_pdb("data/raw/1LYZ.pdb")
    sequences=extract_sequence(structure)
    for chain_id, sequence in sequences.items():
        print(f"Chain {chain_id}")
        print(f"Sequence: {sequence}")
        print("Length:", len(sequence))
if __name__ == "__main__":
    main()