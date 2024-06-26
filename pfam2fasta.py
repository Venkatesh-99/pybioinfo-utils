from Bio import AlignIO
from Bio import SeqIO

def pfam2fasta(input_file, output_file):
    """
    Convert a PFAM alignment file to FASTA format.

    :param input_file: Path to the input PFAM alignment file.
    :param output_file: Path to the output FASTA file.
    """
    alignment = AlignIO.read(input_file, 'stockholm')
    for rec in alignment:
        print(rec.id[:-10], end=" ")
    with open(output_file, 'w') as outfile:
        SeqIO.write(alignment, outfile, 'fasta')
