from Bio import SeqIO
import os


def convert_seq_to_fasta(input_file, output_file):
    """
    Convert a .seq file to FASTA format.

    :param input_file: Path to the input .seq file.
    :param output_file: Path to the output FASTA file.
    """
    with open(output_file, "w") as fasta_output:
        SeqIO.convert(input_file, "seq", fasta_output, "fasta")

def batch_convert_seq_to_fasta(input_dir, output_dir):
    """
    Batch convert .seq files in a directory to FASTA format.

    :param input_dir: Path to the input directory containing .seq files.
    :param output_dir: Path to the output directory for FASTA files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".seq"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, f"{filename[:-4]}.fasta")  # Change file extension to .fasta

            # Convert .seq file to FASTA
            convert_seq_to_fasta(input_file, output_file)

    print("Conversion completed.")

