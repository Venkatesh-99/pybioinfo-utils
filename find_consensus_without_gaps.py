from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
import os




def remove_gaps(sequence):
    """
    Remove gaps ('-') from the sequence
    """
    # Remove gaps ('-') from the sequence
    return sequence.replace('-', '')

def find_consensus(input_dir, output_file):
    """
    Function to find consensus sequences from alignment files in the input directory and write them to the output file.

    :param input_dir: The directory containing alignment files
    :param output_file: The file to write the consensus sequences to
    :return: None
    """
    consensus_sequences = []

    # Iterate over each alignment file in the input directory
    for alignment_file in os.listdir(input_dir):
        if alignment_file.endswith(".fasta"):
            input_file = os.path.join(input_dir, alignment_file)

            # Read the alignment
            alignment = AlignIO.read(input_file, "fasta")

            # Calculate consensus sequence
            summary = SummaryInfo(MultipleSeqAlignment(alignment))
            consensus = summary.dumb_consensus(0.7, "-")

            # Remove gaps from the consensus sequence
            gap_removed_consensus = remove_gaps(str(consensus))

            # Extract the filename without extension
            filename_noext = os.path.splitext(alignment_file)[0][:2]

            # Define the consensus sequence with its identifier
            consensus_seq = f">{filename_noext}\n{gap_removed_consensus}"

            # Append to the list of consensus sequences
            consensus_sequences.append(consensus_seq)

    # Write all gap-removed consensus sequences to the output file
    with open(output_file, "a") as f:
        f.write("\n".join(sorted(consensus_sequences)))



print("Consensus sequence calculation and concatenation completed.")
