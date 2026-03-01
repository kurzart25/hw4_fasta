"""
Тестовый скрипт для проверки всех трёх модулей на example_data.
Запускать из корня проекта: python test_all.py
"""

from pathlib import Path

# --- пути ---
DATA = Path("example_data")
INPUT_FASTA   = DATA / "example_multiline_fasta.fasta"
INPUT_FASTQ   = DATA / "example_fastq.fastq"
INPUT_BLAST   = DATA / "example_blast_results.txt"
INPUT_GBK     = DATA / "example_gbk.gbk"

OUTPUT_FASTA  = DATA / "out_oneline.fasta"
OUTPUT_FASTQ  = DATA / "out_filtered.fastq"
OUTPUT_BLAST  = DATA / "out_blast_hits.txt"
OUTPUT_GBK    = DATA / "out_neighbors.fasta"


def separator(title: str) -> None:
    print(f"\n{'='*50}")
    print(f"  {title}")
    print('='*50)


# ============================================================
# 1. main.py — классы последовательностей
# ============================================================
separator("1. Biological sequences (main.py)")

from main import DNASequence, RNASequence, AminoAcidSequence

# DNASequence
dna = DNASequence("ATGCNATGC")
print(f"DNA:               {dna}")
print(f"complement:        {dna.complement()}")
print(f"reverse:           {dna.reverse()}")
print(f"reverse_complement:{dna.reverse_complement()}")
print(f"transcribe:        {dna.transcribe()}")
print(f"type(transcribe):  {type(dna.transcribe())}")  # должно быть RNASequence
print(f"slice [1:4]:       {dna[1:4]}")
print(f"type(slice):       {type(dna[1:4])}")          # должно быть DNASequence
print(f"index [0]:         {dna[0]}")
print(f"len:               {len(dna)}")

# RNASequence
rna = RNASequence("AUGCN")
print(f"\nRNA:               {rna}")
print(f"complement:        {rna.complement()}")
print(f"type(complement):  {type(rna.complement())}")  # должно быть RNASequence

# AminoAcidSequence
prot = AminoAcidSequence("ACDEFGHIKLM")
print(f"\nProtein:           {prot}")
print(f"aa_composition:    {prot.aa_composition()}")

# Проверка что NucleicAcidSequence нельзя создать напрямую
from main import NucleicAcidSequence
try:
    NucleicAcidSequence("ATGC")
    print("\n[FAIL] NucleicAcidSequence должен кидать NotImplementedError!")
except (NotImplementedError, TypeError) as e:
    print(f"\n[OK] NucleicAcidSequence напрямую: {type(e).__name__}: {e}")

# Проверка невалидного алфавита
try:
    DNASequence("ATGCX")
    print("[FAIL] DNASequence('ATGCX') должен кидать ValueError!")
except ValueError as e:
    print(f"[OK] Невалидный алфавит ДНК: {e}")

# Проверка пустого среза
try:
    dna[5:5]
    print("[FAIL] Пустой срез должен кидать ValueError!")
except ValueError as e:
    print(f"[OK] Пустой срез: {e}")


# ============================================================
# 2. filter_fastq.py
# ============================================================
separator("2. FastQ filter (filter_fastq.py)")

from filter_fastq import filter_fastq

# широкие фильтры — должны пропустить всё
n = filter_fastq(
    str(INPUT_FASTQ), str(OUTPUT_FASTQ),
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
)
print(f"[широкие фильтры]  записей сохранено: {n}")

# строгие фильтры
n_strict = filter_fastq(
    str(INPUT_FASTQ), str(OUTPUT_FASTQ),
    gc_bounds=(40, 60),
    length_bounds=(50, 200),
    quality_threshold=20,
)
print(f"[строгие фильтры]  записей сохранено: {n_strict}")

# gc_bounds как одно число
n_single = filter_fastq(
    str(INPUT_FASTQ), str(OUTPUT_FASTQ),
    gc_bounds=50,
    quality_threshold=10,
)
print(f"[gc_bounds=50]     записей сохранено: {n_single}")

print(f"Выходной файл: {OUTPUT_FASTQ} (exists={OUTPUT_FASTQ.exists()})")


# ============================================================
# 3. bio_files_processor_oop.py
# ============================================================
separator("3a. FASTA oneline (bio_files_processor_oop.py)")

from bio_files_processor_oop import (
    FastaOnelineConverter, FastaOnelineConfig,
    BlastTopHitExtractor, BlastParseConfig,
    GenBankNeighborsExtractor, GenBankNeighborsConfig,
)

conv = FastaOnelineConverter()
out_fa = conv.convert(INPUT_FASTA, OUTPUT_FASTA)
print(f"Выходной файл: {out_fa} (exists={out_fa.exists()})")

# проверим что в выходном файле нет многострочных последовательностей
lines = out_fa.read_text().splitlines()
seq_lines = [l for l in lines if not l.startswith(">")]
multiline = any(len(l) < 10 for l in seq_lines if l)  # очень короткие — признак разбивки
print(f"Строк заголовков: {sum(1 for l in lines if l.startswith('>'))}")
print(f"Строк последовательностей: {len(seq_lines)}")


separator("3b. BLAST top hits")

blast = BlastTopHitExtractor()
try:
    hits = blast.parse(INPUT_BLAST)
    print(f"Найдено хитов: {len(hits)}")
    for h in hits:
        print(f"  query {h.query_index}: {h.description[:80]}")
    out_bl = blast.write(OUTPUT_BLAST)
    print(f"Выходной файл: {out_bl} (exists={out_bl.exists()})")
except ValueError as e:
    print(f"[WARN] {e}")


separator("3c. GenBank neighbors")

gbk = GenBankNeighborsExtractor(GenBankNeighborsConfig(n_before=1, n_after=1))
try:
    out_gbk = gbk.extract(
        INPUT_GBK,
        genes=["phrB", "IFLAKNEJ_00001"],
        output_fasta=OUTPUT_GBK,
    )
    print(f"Выходной файл: {out_gbk} (exists={out_gbk.exists()})")
    print(f"Всего CDS в файле: {len(gbk.last_all_genes)}")
    print(f"Выбранные индексы: {gbk.last_selected_indices}")
    print("Соседи:")
    for i in gbk.last_selected_indices:
        g = gbk.last_all_genes[i]
        print(f"  [{i}] gene={g.gene}, locus_tag={g.locus_tag}, product={g.product}")
except ValueError as e:
    print(f"[WARN] {e}")


separator("Готово!")