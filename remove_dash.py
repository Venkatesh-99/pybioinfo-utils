def remove_dash(input_file, output_file):
    """
    Remove dash characters '-' from protein sequences in a FASTA file.

    :param input_file: Path to the input FASTA file.
    :param output_file: Path to the output file where cleaned sequences will be saved.
    """
    with open(input_file, 'r') as f_in:
        with open(output_file, 'w') as f_out:
            sequence = ''
            for line in f_in:
                if line.startswith('>'):
                    if sequence:
                        clean_sequence = sequence.replace('-', '')
                        f_out.write(clean_sequence + '\n')
                        sequence = ''
                    f_out.write(line)
                else:
                    sequence += line.strip()
            # Write the last sequence
            if sequence:
                clean_sequence = sequence.replace('-', '')
                f_out.write(clean_sequence + '\n')

print("Dash characters removed successfully!")
