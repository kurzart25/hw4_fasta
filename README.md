# ib_nucleotide_utils

A lightweight educational bioinformatics toolkit written in pure Python.  
It provides two main utilities:

1. **DNA/RNA sequence processing** (`run_dna_rna_tools`)  
2. **FASTQ record filtering** (`filter_fastq`)
3. **Utilities for input BLAST/FASTA/GenBank file ('bio_files_processor.py')


The project is designed for training purposes and uses no external dependencies.

---

## Overview

This package contains tools for simple nucleotide operations and FASTQ read filtering.  
It demonstrates modular design, documentation practices (docstrings, typing), and data validation.

### Features
Basic operations with nucleotide sequences, FASTQ quality metrics, and on-the-fly streaming filtering are supported. A separate module has been added for working with BLAST output, multiline FASTA, and extraction of neighboring genes from GenBank.

-Reverse/Complement/Transcription for DNA & RNA

-GC-percentage

-Phred+33 decoding and average read quality

-Filtering by GC%, length and average quality

-Structural validation of FASTQ records

-Converting a multiline FASTA → "one line per sequence"

-Parsing text BLAST output and extracting top hits

-Extracting neighboring genes from GenBank in FASTA

---

## Project structure
```
hw4_fasta/
│
├── example_data/
│   ├── example_blast_results.txt
│   ├── example_fastq.fastq
│   ├── example_gbk.gbk
│   └── example_multiline_fasta.fasta
│
├── modules/
│   ├── fastq_utils.py         
│   └── rna_dna_tools.py       
│
├── bio_files_processor.py     
├── filter_fastq.py           
└── README.md
```

---

## Installation:
```bash
git clone https://github.com/kurzart25/hw4_fasta/tree/hw5_fastq
cd hw4_fasta
```
---

## Usage examples
### DNA/RNA tools
```python
from modules.rna_dna_tools import run_dna_rna_tools
print(run_dna_rna_tools("ATGC", "reverse"))          
print(run_dna_rna_tools("ATGC", proc="transcribe")) 
print(run_dna_rna_tools("AUGC", proc="complement"))  
```

### FASTQ filtering

1) filter_fastq
```python
from filter_fastq import filter_fastq
from modules.example_data import EXAMPLE_FASTQ

input_fastq="example_data/example_fastq.fastq",
    output_fastq="filtered_example.fastq",
    gc_bounds=(40, 60),
    length_bounds=(10, 100),
    quality_threshold=20.0,
    on_duplicate="skip",     # 3 options: "skip" | "rename" | "collect"
    dup_dir="duplicates"     # for duplicates
)
print("Written to:", out_path)

print("Kept reads:", list(filtered.keys()))
```
2) filter_fastq_stream

```python
from filter_fastq import filter_fastq_stream
out_path = filter_fastq_stream(
    input_fastq="example_data/example_fastq.fastq",
    output_fastq="filtered_stream.fastq",
    gc_bounds=(40, 60),
    length_bounds=(10, 100),
    quality_threshold=20.0,
    strict_headers=True  # True allows to exclude corrupted headers
)
print("Written to:", out_path)
```
### Input file utilities (bio_files_processor.py)
```python
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
    select_genes_from_gbk_to_fasta,
)

# Multiline FASTA
one_line_fasta = convert_multiline_fasta_to_oneline(
    "example_data/example_multiline_fasta.fasta"
)
print(one_line_fasta)  

# BLAST - top hits extraction from .txt
top_hits_txt = parse_blast_output(
    input_file="example_data/example_blast_results.txt",
    output_file="blast_top_hits.txt",
)
print(top_hits_txt)

# GenBank: neighbouring gene extraction to FASTA file
gene_selection = select_genes_from_gbk_to_fasta(
    input_gbk="example_data/example_gbk.gbk",
    genes=["rpoB", "gene1234"],  # имя гена или locus_tag
    n_before=1,
    n_after=1,
    output_fasta="selected_genes.fasta",
)
print(gene_selection)
```
---

## Example output
Demo run_dna_rna_tools('ATGC', 'reverse'): CGTA
Kept reads: ['@SRX079802', '@SRX079803', '@SRX079804', ...]

---

## Notes
Phred+33 scoring scheme: score = ord(char) - 33

All bounds are inclusive (e.g., GC between 40% and 60%)

Invalid FASTQ entries are skipped silently

filter_fastq_stream is focused on large files and low memory consumption.

---

## Author
Artem Stetoi (kurzart25)
GitHub: https://github.com/kurzart25