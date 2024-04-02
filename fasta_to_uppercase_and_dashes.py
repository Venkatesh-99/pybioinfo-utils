def fasta_to_uppercase_and_dashes(input_file, output_file):
    """
    Convert a FASTA file to uppercase and replace dots with dashes.

    :param input_file: Path to the input FASTA file.
    :param output_file: Path to the output file where modified sequences will be saved.
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):  # Header line
                f_out.write(line)
            else:
                sequence = line.strip().upper().replace('.', '-')
                f_out.write(sequence + '\n')

print("Converted FASTA file to uppercase and replaced '.' with '-'")
