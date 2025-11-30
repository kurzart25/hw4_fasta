from typing import Tuple, Union, Optional
from modules.fastq_utils import (
    FastqDict,
    FastqRecordError,
    read_fastq_to_dict,
    write_fastq_from_dict,
    normalize_bounds_float,
    normalize_bounds_int,
    validate_fastq_entry,
    gc_content_percent,
    mean_phred33,
    iter_fastq,            
    write_fastq_record,   
    _unique_filename_in,   
)

def filter_fastq_dict(
    seqs: FastqDict,
    gc_bounds: Union[Tuple[float, float], float] = (0.0, 100.0),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0.0,
) -> FastqDict:
    """Filter reads in memory by GC content, length, and mean Phred+33 quality."""
    gc_lo, gc_hi = normalize_bounds_float(gc_bounds)
    len_lo, len_hi = normalize_bounds_int(length_bounds)

    kept: FastqDict = {}
    for name, (seq, qual) in seqs.items():
        try:
            validate_fastq_entry(name, seq, qual)
        except FastqRecordError:
            continue

        if not (len_lo <= len(seq) <= len_hi):
            continue

        gc = gc_content_percent(seq)
        if not (gc_lo <= gc <= gc_hi):
            continue

        if mean_phred33(qual) < quality_threshold:
            continue

        kept[name] = (seq, qual)
    return kept


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0.0, 100.0),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0.0,
    *,
    on_duplicate: str = "skip",
    dup_dir: Optional[str] = None,
) -> str:
    """
    Main function according to the assignment:
      - Reads a FASTQ file (`input_fastq`) into a dictionary (using the selected duplicate-handling policy),
      - Filters reads by GC content, sequence length, and mean Phred+33 quality,
      - Writes the filtered results to 'filtered/<output_fastq>' (with overwrite protection),
      - Returns the full path to the generated file.
    """
    seqs = read_fastq_to_dict(
        input_fastq,
        on_duplicate=on_duplicate,
        dup_dir=dup_dir,
    )
    kept = filter_fastq_dict(
        seqs,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )
    return write_fastq_from_dict(kept, output_fastq)
# я надеюсь это имелось ввиду под "на лету".
"""  
      1. Read a FASTQ record (4 lines: name, sequence, '+', quality)
      2. Validate the structure and compute GC%, length, and mean quality
      3. If the record passes all thresholds → write it to output
      4. Continue until end-of-file
      I hope its memory efficient enough
      Drawback: cannot handle duplicate names"""
def filter_fastq_stream(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0.0, 100.0),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0.0,
    *,
    strict_headers: bool = True,
) -> str:
    gc_lo, gc_hi = normalize_bounds_float(gc_bounds)
    len_lo, len_hi = normalize_bounds_int(length_bounds)

    out_path = _unique_filename_in("filtered", output_fastq)
    total, kept = 0, 0

    with open(out_path, "w", encoding="utf-8") as out_fh:
        for name, seq, qual in iter_fastq(input_fastq):
            total += 1
            try:
                validate_fastq_entry(name, seq, qual)
            except FastqRecordError as e:
                if strict_headers:
                    raise
                continue

            if not (len_lo <= len(seq) <= len_hi):
                continue
            gc = gc_content_percent(seq)
            if not (gc_lo <= gc <= gc_hi):
                continue
            if mean_phred33(qual) < quality_threshold:
                continue

            write_fastq_record(out_fh, name, seq, qual)
            kept += 1

    return out_path