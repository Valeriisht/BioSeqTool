- [**English version**](#NASeqTool)

  
# NASeqTool

<img align=right src="https://github.com/user-attachments/assets/70983c21-48f4-41d0-a7dc-9d5b63982eff" alt="#NASeqTool" width="100"/>

**NASeqTool** is a bioinformatics tool  for analyzing Biological sequences.

The repository consists of a main script are used for DNA/RNA sequence processing and sequencing data processing. The scripts are written in the Python programming language.

## Installation

Clone the repository to a local device.

```sh
git clone git@github.com:Valeriisht/NASeqTool.git
```

The utility was designed to run on Python version 3.9 or newer.

- ### NASeqTool
    - [**Biological sequences**](#biological_sequences) 
    - [**filter_fastq_tool**](#filter_fastq_tool)
    - [**bio_files_processor**](#bio_files_processor)

A utility for manipulating DNA/RNA molecules. The script is written using object-oriented programming.

## Biological sequences

### Usage
You will need:

- DNA/RNA sequences are list (or string when performing an operation on a single sequence), which are required parameters for function calls.

### The utility has five Classes: 

- Polymer Sequence
  - Nuclein Acid Sequence
    - DNA Sequence
    - RNA Sequence 
  - Amino Acid Sequence 

### The utility can perform the following operations: 

- ### Function ```reverse``` in Nuclein Acid Sequence class

Returns the complementary DNA strand in the 3' to 5' direction. For an RNA molecule, cDNA is returned.

- ### Function ```reverse_complement``` in Nuclein Acid Sequence class

Returns the complementary DNA/cDNA strand in the 5' to 3' direction.

- ### Function ```gc_content``` in Nuclein Acid Sequence class

Counts the GC composition of the sequence and returns the percentage of GC nucleotides.

- ### Function ```Motif``` in Amino Sequence class.

 Finds a specific motif in a protein molecule 

- ### Function ```transcribe``` in DNA Sequence class.

The meaningful DNA sequence is translated into an mRNA molecule. 

- ### Function ```is_na_sequence``` in RNA Sequence class.

Checks whether the sequence is DNA or RNA.

- ### Function ```length_cds``` in RNA Sequence class

Determines the length of the CDS coding sequence.


## filter_fastq_tool

The utility for analyzing sequencing data in fastq format. The function is written using the functionality from BioPython

Using:

- SeqIO is a Biopython module that allows to read and write sequence files in various formats such as FASTA, FASTQ, GenBank
- SeqUtils provides various functions for working with biological sequences (includes tools for analyzing and manipulating DNA, RNA, and protein sequences - GC-content, weight)
- SeqRecord is an object that contains a sequence and information about it. It allows the storage of biological sequences

### The utility performs the steps below:

- Filtering the sequence with regard to gc_composition. The gc- composition is counted and determines if the value is within the thresholds.
By default, the threshold values are defined from 0 to 100.

- Filtering the sequence with regard to length. Checks whether the length of the sequence is within the specified threshold values.
By default, the threshold values are defined from 0 to 2**32.

- Sequence filtering with respect to sequence quality. Checks whether the quality on the phred33 scale is within the specified threshold values. The sequence quality is translated from phred33 scale to ASCII.
By default, the threshold value is 0.

The result of the program is an output_file containing only those reads that have passed each of the three filters.

## bio_files_processor

This utility is located in the root of the repository and is designed to work with large .gbk, .fasta, .txt files.

#### Usage

You will need:

- A .fasta file for the convert_multiline_fasta_to_oneline function.
- A .txt file for the parse_blast_output function.
- A .gbk file for the select_genes_from_gbk_to_fasta function.

### The utility performs the following actions:

- ### Function ```convert_multiline_fasta_to_oneline```

Converts fasta format from multiline to single line format

The function accepts 2 arguments (input_fasta and output_fasta - the default parameter). Reads the input fasta file, saves it to a new fasta file where each sequence is written to a single line. 

- ### Function ```parse_blast_output```

Parses the file resulting after alignment in the Blast program. Pulling the name of the best match with the database and saving the results in a single file.

The function takes 2 arguments (input_file, output_file) as input. The program reads the given txt file, for each QUERY best match against the base on the first line from the Description column. The set of obtained proteins is saved to a new file with one column sorted alphabetically.

- ### Function ```select_genes_from_gbk_to_fasta```

Detect a certain number of genes before and after each gene of interest and save their protein sequence into a new fasta file. The fasta-file can be sent to the input of BLAST explicitly. 

The function takes as input input_gbk (path to the input GBK file), genes (genes of interest near which neighbors are searched), n_before, n_after (number of genes before and after (>=1)), output_fasta (name of the output file).

# NASeqTool_

<img align=right src="https://github.com/user-attachments/assets/612014fa-fa0d-4102-a23c-0ef233ae3b48" alt="#NASeqTool" width="100"/>



