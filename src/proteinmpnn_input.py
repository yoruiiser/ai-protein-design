from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


class ProteinChainSelect(Select):

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):

        if not is_aa(residue, standard=True):
            return False

        required_atoms = ["N", "CA", "C"]

        for atom in required_atoms:
            if atom not in residue:
                return False

        return True


def prepare_protein_structure(
    input_file,
    output_file,
    chain_id="A"
):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", input_file)

    model = structure[0]

    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found.")

    chain = model[chain_id]

    valid_residues = 0
    rejected_residues = 0

    for residue in chain:

        if is_aa(residue, standard=True):

            required_atoms = ["N", "CA", "C"]

            if all(atom in residue for atom in required_atoms):
                valid_residues += 1
            else:
                rejected_residues += 1

    io = PDBIO()
    io.set_structure(structure)

    io.save(
        output_file,
        ProteinChainSelect(chain_id)
    )

    print("\nPROTEIN PREPARATION REPORT")
    print("Chain:", chain_id)
    print("Valid residues:", valid_residues)
    print("Rejected residues:", rejected_residues)
    print("Output file:", output_file)