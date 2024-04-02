def count_sequences_in_fasta(fasta_file):
    """
    Counts the number of sequences in the given FASTA file and returns the count.
    
    Parameters:
    fasta_file (str): The path to the FASTA file to be processed.
    
    Returns:
    int: The number of sequences in the FASTA file.
    """
    sequence_count = 0

    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                sequence_count += 1
    
    
    return f"The number of sequences in {fasta_file} is: {sequence_count}"


