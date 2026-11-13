from Bio import SeqIO


def analyze_generated_sequences(fasta_file):

    records = list(SeqIO.parse(fasta_file, "fasta"))

    if len(records) < 2:
        raise ValueError("Expected a reference and generated sequences.")

    reference = str(records[0].seq)

    print("\nPROTEINMPNN SEQUENCE ANALYSIS")
    print("Reference length:", len(reference))

    for record in records[1:]:

        sequence = str(record.seq)

        if len(sequence) != len(reference):
            print(f"{record.id}: Length mismatch")
            continue

        matches = sum(
            a == b
            for a, b in zip(reference, sequence)
        )

        identity = matches / len(reference) * 100

        mutations = len(reference) - matches

        print(f"\n{record.description}")
        print("Length:", len(sequence))
        print("Identical positions:", matches)
        print("Sequence identity:", round(identity, 2), "%")
        print("Different positions:", mutations)