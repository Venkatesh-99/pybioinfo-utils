def align_sequences(input_file, output_file):
    """
    Align sequences by padding or truncating each sequence to the maximum length.

    :param input_file: Path to the input file containing sequences.
    :param output_file: Path to the output file where aligned sequences will be saved.
    """
    sequences = []
    headers = []

    # Read sequences and headers from input file
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                headers.append(line.strip())
            else:
                sequences.append(line.strip())

    # Find the maximum length among all sequences
    max_length = max(len(seq) for seq in sequences)

    # Pad or truncate each sequence to the maximum length
    aligned_sequences = [seq.ljust(max_length, '-')[:max_length] for seq in sequences]

    # Write aligned sequences to the output file
    with open(output_file, 'w') as f:
        for header, seq in zip(headers, aligned_sequences):
            f.write(header + '\n')
            f.write(seq + '\n')

print("Made sequences the same length")
