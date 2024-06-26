# pybioinfo-utils

---
pybioinfo-utils is a collection of small Python functions designed for bioinformatics applications, particularly focused on protein sequence processing. These functions aim to simplify common tasks such as sequence manipulation, file format conversion, and sequence analysis.

## Functions

### `remove_gaps.py`

This function removes dash characters ("-") from protein sequences in a FASTA file and writes the cleaned sequences to a new file.

```python
remove_gaps(input_file, output_file)
```

### `fasta_to_uppercase_and_dashes.py`

This function converts lowercase characters to uppercase and replaces dots (".") with dashes ("-") in a protein sequence.

```python
fasta_to_uppercase_and_dashes(input_file, output_file)
```

### `aa_distribution.py`

This function plots the amino acid distribution for each position in the Multiple Sequence Aignment of protein sequences from the given input FASTA file.

```python
aa_distribution(fasta_file)
```

### `remove_columns_with_high_gap_frequency`

This function removes columns with a high gap frequency from a nested list of protein sequences.

```python
new_nested_list = remove_columns_with_high_gap_frequency(nested_list, threshold=0.5)
```

### `write_to_fasta`

This function writes sequences from a nested list to a FASTA file.

```python
write_to_fasta(output_file, data)
```

### `convert_fasta_to_clustal.py`

This Python script contains two functions. One to convert a FASTA file to Clustal format and the other to convert a directory of FASTA files to Clustal format files.

```python
fasta_to_clustal(input_fasta, output_clustal)
```
```python
convert_directory(input_dir, output_dir)
```

### `convert_seq_to_fasta.py`

This Python script contains two functions. One to convert a file in seq format to FASTA format and the other to convert a directory of seq files to FASTA files.

```python
convert_seq_to_fasta(input_file, output_file)
```
```python
batch_convert_seq_to_fasta(input_dir, output_dir)
```

### `count_sequences_in_fasta.py`

This function counts the number of sequences in a FASTA file.

```python
sequence_count = count_sequences_in_fasta(fasta_file)
```

### `find_consensus_without_gaps.py`

This script contains a function that calculates consensus sequences from multiple sequence alignments and writes them to an output file.

```python
find_consensus(input_dir, output_file)
```

### `make_seq_same_length.py`

This script contains a function that 

### `pfam2fasta`

This function converts PFAM alignment files to FASTA format.

```python
pfam2fasta(input_file, output_file)
```

## Usage

To use these functions, simply import them into your Python script or interactive session and provide the required input parameters. See individual function descriptions for usage examples.
