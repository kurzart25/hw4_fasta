"""
fastq_utils.py

Auxiliary utilities for working with FASTQ data:
- normalization of filtering intervals (GC%, read length),
- calculation of GC content (%),
- calculation of average base quality (Phred+33).
"""

from typing import Tuple, Union


class FastqRecordError(ValueError):
    """Indicates an invalid FASTQ record structure """
    pass


def normalize_bounds_float(
    bounds: Union[Tuple[float, float], float],
    default_hi: float = 100.0
) -> Tuple[float, float]:
    """
    Normalizes a floating-point interval for GC content filtering.
    If a single number X is provided, it is interpreted as (0.0, X).
    Ensures that the lower bound does not exceed the upper bound.

    Args:
        bounds: Either (lo, hi) or a single number X (interpreted as the upper bound).
        default_hi: Default upper bound (100.0), corresponding to the maximum possible GC%.

    Returns:
        (lo, hi): Tuple of ordered bounds (lo <= hi).
    """
    if isinstance(bounds, (int, float)):
        lo, hi = 0.0, float(bounds)
    else:
        lo, hi = float(bounds[0]), float(bounds[1])
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


def normalize_bounds_int(
    bounds: Union[Tuple[int, int], int],
    default_hi: int = 2 ** 32
) -> Tuple[int, int]:
    """
    Normalizes an integer interval for read length filtering.
    If a single number N is provided, it is interpreted as (0, N).
    Ensures that the lower bound does not exceed the upper bound.

    Args:
        bounds: Either (lo, hi) or a single number N (interpreted as the upper bound).
        default_hi: Default upper bound (2**32), representing an effectively unlimited read length.

    Returns:
        (lo, hi): Tuple of ordered bounds (lo <= hi).
    """
    if isinstance(bounds, int):
        lo, hi = 0, int(bounds)
    else:
        lo, hi = int(bounds[0]), int(bounds[1])
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


def gc_content_percent(seq: str) -> float:
    """
    Calculates the GC content of a nucleotide sequence in percentage.

    Args:
        seq: A DNA/RNA sequence. Alphabet validity is not checked.

    Returns:
        GC percentage (0.0 - 100.0). Returns 0.0 for an empty string.
    """
    if not seq:
        return 0.0
    gc = sum(c in "GgCc" for c in seq)
    return 100.0 * gc / len(seq)


def phred33_score(ch: str) -> int:
    """
    Converts a single ASCII character from a quality string to its Phred+33 score.

    Args:
        ch: A single character from a quality string.

    Returns:
        int: ord(ch) - 33
    """
    return ord(ch) - 33


def mean_phred33(quality: str) -> float:
    """
    Calculates the average base quality from a Phred+33 encoded quality string.

    Each character in a FASTQ quality string represents an integer quality score
    according to the Phred+33 scheme, where:
        score = ord(char) - 33

    Args:
        quality: A quality string whose length equals the corresponding read length.

    Returns:
        float: The mean Phred score. Returns 0.0 if the string is empty.
    """
    if not quality:
        return 0.0
    base = 33  
    total = 0  
    for ch in quality:
        total += ord(ch) - base
    return total / len(quality)


def validate_fastq_entry(name: str, seq: str, qual: str) -> None:
    """
    Performs a basic structural validation of a FASTQ record:
    - name must be a string,
    - seq and qual must be strings,
    - len(seq) == len(qual).

    Args:
        name: Record identifier (read name).
        seq: Nucleotide sequence.
        qual: Quality string.

    Raises:
        FastqRecordError: If the record structure is invalid or the lengths do not match.
    """
    if not isinstance(name, str):
        raise FastqRecordError("FASTQ key (name) must be a string.")
    if not isinstance(seq, str) or not isinstance(qual, str):
        raise FastqRecordError(f"{name}: sequence and quality must be strings.")
    if len(seq) != len(qual):
        raise FastqRecordError(
            f"{name}: length mismatch: len(seq)={len(seq)} != len(qual)={len(qual)}"
        )
