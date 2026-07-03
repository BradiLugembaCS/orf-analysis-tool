# ORF Analysis Tool

A Python tool for analysing DNA sequences from FASTA files, translating all six reading frames, and identifying potential open reading frames.

## Overview

This project analyses DNA sequences and identifies realistic open reading frames based on a configurable minimum amino acid length.

The tool reads a FASTA file, translates the DNA sequence into amino acid sequences across all three forward reading frames and all three reverse-complement reading frames, then reports potential ORFs that begin with a start codon and end with a stop codon.

This project demonstrates algorithmic problem solving, biological sequence processing, file parsing, string manipulation, and structured output generation using Python.

## Features

* Reads DNA sequences from FASTA files
* Parses descriptor lines and sequence data
* Removes whitespace and normalises DNA input
* Translates DNA codons into amino acids
* Analyses all three forward reading frames
* Generates the reverse-complement DNA strand
* Analyses all three reverse-complement reading frames
* Identifies ORFs beginning with a start codon and ending with a stop codon
* Allows the user to set a minimum amino acid length threshold
* Reports nucleotide positions, amino acid positions, ORF lengths, DNA sequences, and amino acid sequences

## Technologies Used

* Python
* Tkinter file selection
* FASTA file parsing
* Dictionary-based codon translation
* String processing
* Algorithmic sequence analysis

## How It Works

The program uses a codon translation table to convert DNA triplets into amino acids. It then searches the translated amino acid sequences for open reading frames that begin with methionine, represented by `M`, and end with a stop codon, represented by `*`.

The analysis is performed across six reading frames:

* `+1`
* `+2`
* `+3`
* `-1`
* `-2`
* `-3`

For each detected ORF, the program reports:

* Frame number
* Nucleotide start and stop positions
* Amino acid start and stop positions
* ORF length in nucleotides
* ORF length in amino acids
* DNA sequence
* Amino acid sequence

## Installation

Clone the repository:

```bash
git clone https://github.com/BradiLugembaCS/orf-analysis-tool.git
```

Navigate into the project folder:

```bash
cd orf-analysis-tool
```

Run the program:

```bash
python orf_analysis_tool.py
```

The program will open a file selection window where you can choose a FASTA file.

## Example FASTA Format

```text
>Example DNA sequence
ATGAAATTTGGGCCCTAA
```

## Project Structure

```text
orf-analysis-tool/
├── examples/
│   └── sample.fasta
├── orf_analysis_tool.py
├── README.md
└── .gitignore
```

## Key Learning Outcomes

Through this project, I strengthened my understanding of:

* Python programming fundamentals
* File input and parsing
* String manipulation
* Dictionary-based lookup tables
* DNA-to-amino-acid translation
* Reading frame analysis
* Algorithm design
* Structured console output


