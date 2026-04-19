from __future__ import annotations

import logging
from pathlib import Path
from typing import Tuple, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

Bounds = Union[Tuple[float, float], float, Tuple[int, int], int]


def setup_logger(log_file: str = "fastq_filter.log") -> logging.Logger:
    logger = logging.getLogger("fastq_filter")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)

    handler = logging.FileHandler(log_file, encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def normalize_bounds(
    bounds: Union[Tuple[float, float], float, Tuple[int, int], int],
    cast: type = float,
) -> tuple:
    if isinstance(bounds, (int, float)):
        lo, hi = cast(0), cast(bounds)
    else:
        if len(bounds) != 2:
            raise ValueError("Bounds must be a single number or a tuple of two values.")
        lo, hi = cast(bounds[0]), cast(bounds[1])

    return (lo, hi) if lo <= hi else (hi, lo)


def passes_filter(
    rec: SeqRecord,
    gc_lo: float,
    gc_hi: float,
    len_lo: int,
    len_hi: int,
    quality_threshold: float,
) -> bool:
    seq_len = len(rec.seq)
    if not (len_lo <= seq_len <= len_hi):
        return False

    gc_pct = 100.0 * gc_fraction(rec.seq)
    if not (gc_lo <= gc_pct <= gc_hi):
        return False

    phred = rec.letter_annotations.get("phred_quality")
    if phred is None or len(phred) == 0:
        return False

    mean_q = sum(phred) / len(phred)
    return mean_q >= quality_threshold


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0.0, 100.0),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0.0,
    log_file: str | None = None,
) -> int:
    logger = setup_logger(log_file)

    input_path = Path(input_fastq)
    output_path = Path(output_fastq)

    if not input_path.exists():
        logger.error(f"Input FASTQ file not found: {input_path}")
        raise FileNotFoundError(f"Input FASTQ file not found: {input_path}")

    gc_lo, gc_hi = normalize_bounds(gc_bounds, float)
    len_lo, len_hi = normalize_bounds(length_bounds, int)

    logger.info(
        f"Starting FASTQ filtering: input={input_path}, output={output_path}, "
        f"gc_bounds=({gc_lo}, {gc_hi}), length_bounds=({len_lo}, {len_hi}), "
        f"quality_threshold={quality_threshold}"
    )

    records = SeqIO.parse(str(input_path), "fastq")
    filtered = (
        rec for rec in records
        if passes_filter(rec, gc_lo, gc_hi, len_lo, len_hi, quality_threshold)
    )

    written = SeqIO.write(filtered, str(output_path), "fastq")
    logger.info(f"Filtering finished successfully. Reads written: {written}")

    return written