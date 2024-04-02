# pybioinfo-utils

---

# pybioinfo-utils

pybioinfo-utils is a collection of small Python functions designed for bioinformatics applications, particularly focused on protein sequence processing. These functions aim to simplify common tasks such as sequence manipulation, file format conversion, and sequence analysis.

## Functions

### `remove_dash`

This function removes dash characters ("-") from protein sequences in a FASTA file and writes the cleaned sequences to a new file.

```python
remove_dash(input_file, output_file)
```

### `fasta_to_uppercase_and_dashes`

This function converts lowercase characters to uppercase and replaces dots (".") with dashes ("-") in a protein sequence.

```python
fasta_to_uppercase_and_dashes(input_file, output_file)
```

### `fasta_to_nested_list`

This function reads sequences from a FASTA file and converts them into a nested list format, where each sequence is represented as a list of amino acids.

```python
nested_list = fasta_to_nested_list(input_file)
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

### `fasta_to_clustal`

This function converts FASTA files to Clustal format.

```python
fasta_to_clustal(input_fasta, output_clustal)
```

### `align_sequences`

This function aligns sequences in a FASTA file by padding or truncating them to the maximum length.

```python
align_sequences(input_file, output_file)
```

### `convert_directory`

This function converts multiple FASTA files in a directory to Clustal format.

```python
convert_directory(input_dir, output_dir)
```

### `count_sequences_in_fasta`

This function counts the number of sequences in a FASTA file.

```python
sequence_count = count_sequences_in_fasta(fasta_file)
```

### `find_consensus`

This function calculates consensus sequences from multiple sequence alignments and writes them to an output file.

```python
find_consensus(input_dir, output_file)
```

### `pfam2fasta`

This function converts PFAM alignment files to FASTA format.

```python
pfam2fasta(input_file, output_file)
```

## Usage

To use these functions, simply import them into your Python script or interactive session and provide the required input parameters. See individual function descriptions for usage examples.
