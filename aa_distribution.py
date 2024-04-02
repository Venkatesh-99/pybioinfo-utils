from Bio import SeqIO
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

def plot_amino_acid_distribution(fasta_file):
    """
    Plot the amino acid distribution for each position in the alignment.

    :param fasta_file: Path to the input FASTA file containing aligned sequences.
    """
    # Read the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))  # Convert iterator to list

    # Initialize a counter to store amino acid frequencies for each position
    amino_acid_counts = []

    # Iterate through each position in the alignment
    for i in range(len(sequences[0].seq)):  # Use the length of the first sequence for alignment length
        # Extract the amino acids at position i for all sequences
        aligned_amino_acids = [record.seq[i] for record in sequences]
        # Count the occurrence of each amino acid
        counts = Counter(aligned_amino_acids)
        # Append the counts to the amino_acid_counts list
        amino_acid_counts.append(counts)

    # Calculate the frequency of each amino acid at each position
    num_sequences = len(sequences)
    amino_acid_distribution = [{aa: count / num_sequences for aa, count in counts.items()} for counts in amino_acid_counts]

    # Define a color palette for amino acids
    color_palette = plt.cm.get_cmap('tab20', len(amino_acid_distribution[0]))

    # Optionally, visualize the amino acid distribution
    positions = np.arange(len(amino_acid_distribution))
    amino_acids = sorted(amino_acid_distribution[0].keys())  # Assuming all positions have the same amino acids
    for idx, aa in enumerate(amino_acids):
        frequencies = [dist.get(aa, 0) for dist in amino_acid_distribution]
        plt.plot(positions, frequencies, label=aa, color=color_palette(idx))

    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.title("Amino Acid Distribution")
    plt.legend(ncol=2)
    plt.show()

