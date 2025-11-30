# ib_nucleotide_utils

A lightweight educational bioinformatics toolkit written in pure Python.  
It provides two main utilities:

1. **DNA/RNA sequence processing** (`run_dna_rna_tools`)  
2. **FASTQ record filtering** (`filter_fastq`)

The project is designed for training purposes and uses no external dependencies.

---

## Overview

This package contains tools for simple nucleotide operations and FASTQ read filtering.  
It demonstrates modular design, documentation practices (docstrings, typing), and data validation.

### Features
- Reverse, complement, and transcription for DNA & RNA
- GC content calculation (%)
- Phred+33 quality decoding and mean score computation
- Read filtering by GC%, length, and mean quality thresholds
- Structural FASTQ record validation

---

## Project structure
```
ib_nucleo_utils/
│
├── modules/
│ ├── rna_dna_tools.py # DNA/RNA utilities
│ ├── fastq_utils.py # GC%, quality, and validation helpers
│ ├── example_data.py # Example dataset for testing
│
├── filter_fastq.py # Main script (entry point)
└── README.md
```

---

## Installation:
```bash
git clone https://github.com/<your-username>/ib_nucleo_utils.git
cd ib_nucleo_utils
python filter_fastq.py
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
```python
from filter_fastq import filter_fastq
from modules.example_data import EXAMPLE_FASTQ

filtered = filter_fastq(
    EXAMPLE_FASTQ,
    gc_bounds=(40, 60),
    length_bounds=(10, 100),
    quality_threshold=20.0,
)

print("Kept reads:", list(filtered.keys()))
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

---

## Author
Artem Stetoi (kurzart25)
GitHub: https://github.com/kurzart25