from Bio import SeqIO
from collections import Counter


def calculate_conservation(fasta_file):

    records = list(SeqIO.parse(fasta_file, "fasta"))

    sequences = [str(record.seq) for record in records]

    if not sequences:
        raise ValueError("No sequences found.")

    lengths = {len(sequence) for sequence in sequences}

    if len(lengths) != 1:
        raise ValueError("Sequences must have equal lengths.")

    sequence_length = len(sequences[0])

    conservation = []

    for position in range(sequence_length):

        amino_acids = [
            sequence[position]
            for sequence in sequences
        ]

        counts = Counter(amino_acids)

        most_common_aa, count = counts.most_common(1)[0]

        frequency = count / len(sequences)

        conservation.append({
            "position": position + 1,
            "most_common_aa": most_common_aa,
            "conservation": frequency
        })

    return conservation