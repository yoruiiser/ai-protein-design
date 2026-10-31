from Bio.PDB import PDBParser
def parse_pdb(path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('Protein', path)
    return structure
def explore_structure(structure):
    print(f"Structure ID: {structure.id}")
    models = list(structure.get_models())
    print(f"Number of models: {len(models)}")
    for model in models:
        print("\nModel ID:", model.id)
        for chain in model:
            print("Chain:", chain.id)
def main():
    path="data/raw/1LYZ.pdb"
    structure=parse_pdb(path)
    explore_structure(structure)
if __name__ == "__main__":
    main()