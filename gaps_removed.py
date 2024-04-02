from Bio import SeqIO

def fasta_to_nested_list(input_file):
    """
    Convert a FASTA file to a nested list format.

    :param input_file: Path to the input FASTA file.
    :return: Nested list containing (sequence_id, sequence_list) tuples.
    """
    sequences = []

    with open(input_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence_list = list(str(record.seq))
            sequences.append((record.id, sequence_list))
    return sequences

def remove_columns_with_high_gap_frequency(nested_list, threshold):
    """
    Remove columns with high gap frequency from a nested list of sequences.

    :param nested_list: Nested list containing (sequence_id, sequence_list) tuples.
    :param threshold: Threshold for gap frequency (default is 0.5).
    :return: Modified nested list with high-gap columns removed.
    """
    num_sequences = len(nested_list)
    sequence_length = len(nested_list[0][1])
    columns_to_remove = []
    
    for i in range(sequence_length):
        gap_count = sum(1 for _, sequence in nested_list if sequence[i] == '-')
        if gap_count / num_sequences > threshold:
            columns_to_remove.append(i)
    
    for i in reversed(columns_to_remove):
        for j, (seq_id, sequence) in enumerate(nested_list):
            nested_list[j] = (seq_id, sequence[:i] + sequence[i+1:])
    
    return nested_list

def write_to_fasta(output_file, data):
    """
    Write sequences in FASTA format to an output file.

    :param output_file: Path to the output FASTA file.
    :param data: Nested list containing (sequence_id, sequence_list) tuples.
    """
    with open(output_file, 'w') as out:
        for seq_id, sequence in data:
            out.write(f">{seq_id}\n{''.join(sequence)}\n")

print("Removed columns with gap frequency above the specified threshold.")
