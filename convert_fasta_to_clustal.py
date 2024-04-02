from Bio import SeqIO
import os



def fasta_to_clustal(input_fasta, output_clustal):
    """
    Converts a FASTA file to Clustal format and saves the result to a file.

    Parameters:
        input_fasta (str): The path to the input FASTA file.
        output_clustal (str): The path to the output Clustal file.

    Returns:
        None
    """
    # Read the FASTA file and convert it to Clustal format
    records = SeqIO.parse(input_fasta, "fasta")
    SeqIO.write(records, output_clustal, "clustal")

def convert_directory(input_dir, output_dir):
    # Iterate over FASTA files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_fasta = os.path.join(input_dir, filename)
            output_clustal = os.path.join(output_dir, filename.replace(".fasta", ".clustal"))

            # Convert the FASTA file to Clustal format
            fasta_to_clustal(input_fasta, output_clustal)



print("Conversion of FASTA files to Clustal format completed.")
