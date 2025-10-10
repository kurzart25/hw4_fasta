"""
fastq_utils.py

Auxiliary utilities for working with FASTQ data:
- normalization of filtering intervals (GC%, read length),
- calculation of GC content (%),
- calculation of average base quality (Phred+33),
- reading/writing FASTQ files (dict form),
- safe output path generation.
"""

import os
from typing import Dict, Tuple, Union


class FastqRecordError(ValueError):
    """Indicates an invalid FASTQ record structure."""
    pass


def normalize_bounds_float(
    bounds: Union[Tuple[float, float], float],
    default_hi: float = 100.0
) -> Tuple[float, float]:
    if isinstance(bounds, bool):
        raise TypeError("bounds must not be a boolean")

    if isinstance(bounds, (int, float)):
        lo, hi = 0.0, float(bounds)
    else:
        lo, hi = float(bounds[0]), float(bounds[1])

    if lo > hi:
        lo, hi = hi, lo

    lo = max(0.0, min(lo, default_hi))
    hi = max(0.0, min(hi, default_hi))
    return lo, hi


def normalize_bounds_int(
    bounds: Union[Tuple[int, int], int],
    default_hi: int = 2 ** 32
) -> Tuple[int, int]:
    if isinstance(bounds, bool):
        raise TypeError("bounds must not be a boolean")

    if isinstance(bounds, int):
        lo, hi = 0, int(bounds)
    else:
        lo, hi = int(bounds[0]), int(bounds[1])

    if lo > hi:
        lo, hi = hi, lo

    lo = max(0, min(lo, default_hi))
    hi = max(0, min(hi, default_hi))
    return lo, hi


def gc_content_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(c in "GgCc" for c in seq)
    return 100.0 * gc / len(seq)


def phred33_score(ch: str) -> int:
    return ord(ch) - 33


def mean_phred33(quality: str) -> float:
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
    - len(seq) == len(qual),
    - quality chars must be in ASCII range [33,126].
    """
    if not isinstance(name, str):
        raise FastqRecordError("FASTQ key (name) must be a string.")
    if not isinstance(seq, str) or not isinstance(qual, str):
        raise FastqRecordError(f"{name}: sequence and quality must be strings.")
    if len(seq) != len(qual):
        raise FastqRecordError(
            f"{name}: length mismatch: len(seq)={len(seq)} != len(qual)={len(qual)}"
        )

    for i, ch in enumerate(qual):
        code = ord(ch)
        if code < 33 or code > 126:
            raise FastqRecordError(
                f"{name}: invalid quality character at pos {i}: {repr(ch)} (ord={code}), "
                "expected ASCII 33..126 (Phred+33)."
            )


FastqDict = Dict[str, Tuple[str, str]]


def _split_base_ext(filename: str) -> Tuple[str, str]:
    """Split filename into (base, ext) using rsplit(., 1); ext includes the leading dot or is empty."""
    filename = os.path.basename(filename).strip()
    if "." in filename:
        base, ext = filename.rsplit(".", 1)
        return base, "." + ext
    return filename, ""


def _duplicates_filename_for(input_filename: str) -> str:
    """Make a duplicates filename from the input file name."""
    base, ext = _split_base_ext(input_filename)
    if not ext:
        ext = ".fq"
    return f"{base}.dups{ext}"


def _unique_filename_in(dirpath: str, filename: str) -> str:
    """Return a path in `dirpath` with `filename`, ensuring no overwrite."""
    os.makedirs(dirpath, exist_ok=True)
    filename = os.path.basename(filename).strip()

    base, ext = _split_base_ext(filename)
    candidate = os.path.join(dirpath, base + ext)
    if not os.path.exists(candidate):
        return candidate

    i = 1
    while True:
        cand = os.path.join(dirpath, f"{base}_{i}{ext}")
        if not os.path.exists(cand):
            return cand
        i += 1


def read_fastq_to_dict(
    path: str,
    on_duplicate: str = "skip",
    dup_dir: Union[str, None] = None
) -> FastqDict:
    """
    Reads a plain-text FASTQ file into a dict: {name: (seq, qual)}.
      - Skips blank lines.
      - on_duplicate: 'skip' (default), 'overwrite', or 'error'.
      - dup_dir: if provided, duplicates are written to a separate FASTQ file in this directory.
    """
    if on_duplicate not in {"skip", "overwrite", "error"}:
        raise ValueError("on_duplicate must be one of {'skip','overwrite','error'}")

    records: FastqDict = {}
    dup_fh = None

    try:
        if dup_dir is not None:
            dup_name = _duplicates_filename_for(path)
            dup_path = _unique_filename_in(dup_dir, dup_name)
            dup_fh = open(dup_path, "a", encoding="utf-8")

        with open(path, "r", encoding="utf-8") as fh:
            while True:
                line1 = fh.readline()
                while line1 and not line1.strip():
                    line1 = fh.readline()
                if not line1:
                    break

                line2 = fh.readline()
                line3 = fh.readline()
                line4 = fh.readline()
                if not line4:
                    raise FastqRecordError("Incomplete FASTQ record at file end.")

                name_line = line1.rstrip("\r\n")
                seq = line2.rstrip("\r\n")
                plus = line3.rstrip("\r\n")
                qual = line4.rstrip("\r\n")

                if not name_line.startswith("@") or not plus.startswith("+"):
                    raise FastqRecordError("Malformed FASTQ headers (@/+) encountered.")

                name = name_line[1:]  # remove '@'
                validate_fastq_entry(name, seq, qual)

                if name in records:
                    if on_duplicate == "error":
                        raise FastqRecordError(f"Duplicate read name encountered: {name}")
                    if dup_fh is not None:
                        dup_fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    if on_duplicate == "skip":
                        continue

                records[name] = (seq, qual)
    finally:
        if dup_fh is not None:
            dup_fh.close()

    return records


def write_fastq_from_dict(records: FastqDict, output_name: str) -> str:
    """
    Writes records {name:(seq,qual)} to filtered/<output_name>,
    with overwrite protection. Returns full path.
    """
    out_name_only = os.path.basename(output_name).strip()
    out_path = _unique_filename_in("filtered", out_name_only)

    with open(out_path, "w", encoding="utf-8") as fh:
        for name, (seq, qual) in records.items():
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

    return out_path

def iter_fastq(path: str):
    """
    Yield FASTQ records (name, seq, qual) one-by-one.
    Skips blank lines; performs basic header checks (@ and +).
    Raises FastqRecordError on incomplete 4-line record at EOF.
    """
    with open(path, "r", encoding="utf-8") as fh:
        while True:
            line1 = fh.readline()
            while line1 and not line1.strip():
                line1 = fh.readline()
            if not line1:
                break  # EOF

            line2 = fh.readline()
            line3 = fh.readline()
            line4 = fh.readline()
            if not line4:
                raise FastqRecordError("Incomplete FASTQ record at file end.")

            name_line = line1.rstrip("\r\n")
            seq = line2.rstrip("\r\n")
            plus = line3.rstrip("\r\n")
            qual = line4.rstrip("\r\n")

            if not name_line.startswith("@") or not plus.startswith("+"):
                raise FastqRecordError("Malformed FASTQ headers (@/+) encountered.")

            name = name_line[1:]  # drop '@'
            yield name, seq, qual


def write_fastq_record(fh, name: str, seq: str, qual: str) -> None:
    """Write a single FASTQ record to an open text file handle."""
    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")