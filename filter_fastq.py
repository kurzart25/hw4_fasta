from typing import Dict, Tuple, Union
from modules.rna_dna_tools import run_dna_rna_tools
from modules import fastq_utils

FastqDict = Dict[str, Tuple[str, str]]


def filter_fastq(
    seqs: FastqDict,
    gc_bounds: Union[Tuple[float, float], float] = (0.0, 100.0),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0.0,
) -> FastqDict:
    """
    Filters FASTQ reads by GC content (percent), read length, and mean Phred+33 quality.

    Intervals are inclusive. If a single number is provided for `gc_bounds` or
    `length_bounds`, it is interpreted as an upper bound (i.e., (0, X)).

    Args:
        seqs:
            Dictionary of FASTQ-like records:
            { read_name: (sequence, quality_string) }.
        gc_bounds:                                
            GC% interval as (lo, hi) or a single float X interpreted as (0.0, X).
            Defaults to (0.0, 100.0).
        length_bounds:
            Length interval as (lo, hi) or a single int N interpreted as (0, N).
            Defaults to (0, 2**32).
        quality_threshold:
            Minimum allowed mean Phred+33 quality. Reads with mean quality
            strictly below this threshold are discarded. Defaults to 0.0.

    Returns:
        A new dictionary with only the reads that satisfy all criteria.
    """
    gc_lo, gc_hi = fastq_utils.normalize_bounds_float(gc_bounds)
    len_lo, len_hi = fastq_utils.normalize_bounds_int(length_bounds)

    kept: FastqDict = {}

    for name, (seq, qual) in seqs.items():
        try:
            fastq_utils.validate_fastq_entry(name, seq, qual)
        except fastq_utils.FastqRecordError:
            continue

        length_seq = len(seq)
        if not (len_lo <= length_seq <= len_hi):
            continue

        gc = fastq_utils.gc_content_percent(seq)
        if not (gc_lo <= gc <= gc_hi):
            continue

        mean_q = fastq_utils.mean_phred33(qual)
        if mean_q < quality_threshold:
            continue

        kept[name] = (seq, qual)

    return kept


if __name__ == "__main__":
    from modules.example_data import EXAMPLE_FASTQ

    print("Demo run_dna_rna_tools('ATGC','reverse'):", run_dna_rna_tools("ATGC", "reverse"))
    filtered = filter_fastq(
        EXAMPLE_FASTQ,
        gc_bounds=(40, 60),
        length_bounds=(10, 100),
        quality_threshold=20.0,
    )
    print("Kept reads:", list(filtered.keys()))
