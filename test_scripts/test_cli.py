from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from filter_fastq import (
    filter_fastq,
    normalize_bounds,
    passes_filter,
)


def make_fastq_record(seq: str, qualities: list[int], record_id: str) -> SeqRecord:
    rec = SeqRecord(Seq(seq), id=record_id, description="")
    rec.letter_annotations["phred_quality"] = qualities
    return rec


@pytest.fixture
def sample_fastq(tmp_path: Path) -> Path:
    records = [
        make_fastq_record("ATGC", [40, 40, 40, 40], "read1"),   # GC=50, len=4, high Q
        make_fastq_record("AAAA", [40, 40, 40, 40], "read2"),   # GC=0, len=4, high Q
        make_fastq_record("GGGGGG", [30] * 6, "read3"),         # GC=100, len=6, medium Q
        make_fastq_record("ATATAT", [10] * 6, "read4"),         # GC=0, len=6, low Q
    ]
    path = tmp_path / "input.fastq"
    SeqIO.write(records, path, "fastq")
    return path


def test_normalize_bounds_tuple_order():
    assert normalize_bounds((60, 40), float) == (40.0, 60.0)


def test_normalize_bounds_single_number():
    assert normalize_bounds(50, float) == (0.0, 50.0)


def test_passes_filter_true():
    rec = make_fastq_record("ATGC", [30, 30, 30, 30], "ok")
    assert passes_filter(rec, 40.0, 60.0, 4, 10, 20.0) is True


def test_passes_filter_false_due_to_quality():
    rec = make_fastq_record("ATGC", [5, 5, 5, 5], "badq")
    assert passes_filter(rec, 40.0, 60.0, 4, 10, 20.0) is False


def test_filter_fastq_by_gc(sample_fastq: Path, tmp_path: Path):
    output = tmp_path / "out_gc.fastq"
    written = filter_fastq(
        str(sample_fastq),
        str(output),
        gc_bounds=(40, 60),
        length_bounds=(0, 100),
        quality_threshold=0,
        log_file=str(tmp_path / "gc.log"),
    )
    assert written == 1


def test_filter_fastq_by_length(sample_fastq: Path, tmp_path: Path):
    output = tmp_path / "out_len.fastq"
    written = filter_fastq(
        str(sample_fastq),
        str(output),
        gc_bounds=(0, 100),
        length_bounds=(5, 10),
        quality_threshold=0,
        log_file=str(tmp_path / "len.log"),
    )
    assert written == 2


def test_filter_fastq_output_file_created(sample_fastq: Path, tmp_path: Path):
    output = tmp_path / "out.fastq"
    written = filter_fastq(
        str(sample_fastq),
        str(output),
        gc_bounds=(0, 100),
        length_bounds=(0, 100),
        quality_threshold=20,
        log_file=str(tmp_path / "write.log"),
    )
    assert output.exists()
    assert written == 3


def test_filter_fastq_missing_input_logs_error_and_raises(tmp_path: Path):
    missing = tmp_path / "missing.fastq"
    output = tmp_path / "out.fastq"
    log_file = tmp_path / "error.log"

    with pytest.raises(FileNotFoundError):
        filter_fastq(
            str(missing),
            str(output),
            log_file=str(log_file),
        )

    assert log_file.exists()
    log_text = log_file.read_text(encoding="utf-8")
    assert "ERROR" in log_text
    assert "Input FASTQ file not found" in log_text