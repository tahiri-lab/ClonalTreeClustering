def read_fasta(file_path):
    """
    Reads a FASTA format file and returns a list of tuples, each containing the header and the sequence.
    """
    sequences = []
    current_sequence = ""
    header = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_sequence:
                    sequences.append((header, current_sequence))
                header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_sequence:
            sequences.append((header, current_sequence))

    return sequences

def count_sequence_repetitions(fasta_sequences):
    """
    Takes a list of tuples (header, sequence) from a FASTA file and returns a dictionary
    where each key is a header and the value is the number of times that sequence is repeated.
    """
    sequence_counts = {}

    for header, sequence in fasta_sequences:
        key = (header, sequence)
        sequence_counts[key] = sequence_counts.get(key, 0) + 1

    return sequence_counts

# The path to FASTA file
file_path = '40_8.txt'
fasta_sequences = read_fasta(file_path)
repetition_counts = count_sequence_repetitions(fasta_sequences)

# Print the header and count of each sequence
for (header, sequence), count in repetition_counts.items():
    print(f"Header: {header}\nSequence:\n{sequence}\nCount: {count}\n")