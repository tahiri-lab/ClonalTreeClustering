
import os

def read_fasta(file_path):
    """
    Reads a FASTA format file and returns a list of tuples containing header and sequence.
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

def read_newick(file_path):
    """
    Read Newick tree from file
    """
    with open(file_path, 'r') as file:
        return file.read().strip()

def count_sequence_repetitions(fasta_sequences):
    """
    Count how many times each sequence appears in the FASTA file.
    """
    sequence_counts = {}
    sequence_to_headers = {}
    
    # Group identical sequences
    for header, sequence in fasta_sequences:
        if sequence not in sequence_to_headers:
            sequence_to_headers[sequence] = []
        sequence_to_headers[sequence].append(header)
    
    # Count repetitions for each header
    for sequence, headers in sequence_to_headers.items():
        count = len(headers)
        for header in headers:
            sequence_counts[header] = count
    
    # Get weight for naive (use the first sequence's count)
    naive_count = len(sequence_to_headers[list(sequence_to_headers.keys())[0]])
    sequence_counts['naive'] = naive_count
    
    return sequence_counts

def main():
    # Set up paths
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal/60")
    fasta_file = os.path.join(base_path, "60_10.fasta")
    newick_file = os.path.join(base_path, "60_10.GT.nk")
    
    # Read files
    fasta_sequences = read_fasta(fasta_file)
    newick_string = read_newick(newick_file)
    
    # Get sequence counts
    repetition_counts = count_sequence_repetitions(fasta_sequences)
    
    # Initialize modified Newick string
    modified_newick = newick_string
    
    # Add weights after sequence names
    for header, count in repetition_counts.items():
        # Handle all sequences including naive
        modified_newick = modified_newick.replace(f"{header}:", f"{header}@{count}:")
    
    # If the Newick string doesn't end with naive node, add it
    if not modified_newick.rstrip(');').endswith('naive'):
        modified_newick = modified_newick.rstrip(';') + ")naive@" + str(repetition_counts['naive']) + ":1;"
    
    # Print only the modified Newick string
    print(modified_newick)

if __name__ == "__main__":
    main()
