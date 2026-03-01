# ib_nucleotide_utils

A lightweight educational bioinformatics toolkit written in Python using OOP principles.  
It provides three main utilities:

1. **Biological sequence classes** (`main.py`) — DNA, RNA, and protein sequence objects
2. **FASTQ filtering** (`filter_fastq.py`) — read filtering via Biopython
3. **Bioinformatics file processing** (`bio_files_processor_oop.py`) — FASTA, BLAST, GenBank

---

## Installation

```bash
git clone https://github.com/kurzart25/hw4_fasta
cd hw4_fasta
pip install -r requirements.txt
```

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
├── test_scripts/
│   └── test.py
│
├── main.py                      # Biological sequence classes  + FilterFASTQ         
├── bio_files_processor_oop.py   # FASTA / BLAST / GenBank utilities
├── requirements.txt
└── README.md
```

---

## Usage

### 1. Biological sequences (`main.py`)

Four classes are provided: `DNASequence`, `RNASequence`, `AminoAcidSequence`,  
all inheriting from the abstract base class `BiologicalSequence`.

```python
from main import DNASequence, RNASequence, AminoAcidSequence

# DNA
dna = DNASequence("ATGCNATGC")
print(dna)                    
print(dna.complement())       
print(dna.reverse())          
print(dna.reverse_complement())
print(dna.transcribe())      
print(dna[1:4])              
print(len(dna))               

# RNA
rna = RNASequence("AUGCN")
print(rna.complement())       

# Protein
prot = AminoAcidSequence("ACDEFGHIKLM")
print(prot.aa_composition())  # {'A': 1, 'C': 1, 'D': 1, ...}
```

All sequences validate their alphabet on creation and raise `ValueError` on invalid symbols.

---

### 2. FASTQ filtering (`filter_fastq.py`)

Filters reads by GC content, length, and mean Phred quality using Biopython.  
Returns the number of reads written.

```python
from filter_fastq import filter_fastq

n = filter_fastq(
    input_fastq="example_data/example_fastq.fastq",
    output_fastq="filtered.fastq",
    gc_bounds=(40, 60),      
    length_bounds=(50, 200),  
    quality_threshold=20.0,   
)
print(f"Reads saved: {n}")
```

---

### 3. File processing utilities (`bio_files_processor_oop.py`)

#### 3a. Multiline FASTA → one-line FASTA

```python
from bio_files_processor_oop import FastaOnelineConverter

conv = FastaOnelineConverter()
out = conv.convert("example_data/example_multiline_fasta.fasta")
print(out)  
```

#### 3b. BLAST top hits extraction

```python
from bio_files_processor_oop import BlastTopHitExtractor

blast = BlastTopHitExtractor()
hits = blast.parse("example_data/example_blast_results.txt")
for h in hits:
    print(h.query_index, h.description)

blast.write("top_hits.txt")
# or in one call:
blast.parse_and_write("example_data/example_blast_results.txt", "top_hits.txt")
```

#### 3c. GenBank neighboring CDS → FASTA

Extracts protein translations of CDS entries neighboring the target gene(s).

```python
from bio_files_processor_oop import GenBankNeighborsExtractor, GenBankNeighborsConfig

gbk = GenBankNeighborsExtractor(GenBankNeighborsConfig(n_before=1, n_after=1))
out = gbk.extract(
    input_gbk="example_data/example_gbk.gbk",
    genes=["phrB", "IFLAKNEJ_00001"],   # gene name or locus_tag
    output_fasta="neighbors.fasta",
)
print(out)  
```

---

## Notes

- All bounds are inclusive
- `gc_bounds` accepts either a tuple `(lo, hi)` or a single float (treated as upper bound from 0)
- `DNASequence`, `RNASequence`, `AminoAcidSequence` are immutable (`frozen=True`)
- `NucleicAcidSequence` is abstract and cannot be instantiated directly

---

## Author

Artem Stetoi (kurzart25)  
GitHub: https://github.com/kurzart25