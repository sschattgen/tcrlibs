"""Smoke tests for the counting core."""

import csv
from pathlib import Path

from tcrlibs.counting import (
    get_reverse_complement,
    load_library,
    process_reads,
    resolve_trim_length,
    write_output,
)

DATA_DIR = Path(__file__).parent / "data"
LIB_CSV = DATA_DIR / "lib.csv"
SAMPLE_FASTQ = DATA_DIR / "sample.fastq"


def test_reverse_complement():
    assert get_reverse_complement("ACGT") == "ACGT"
    assert get_reverse_complement("AAAA") == "TTTT"
    assert get_reverse_complement("ACGTN") == "NACGT"
    assert get_reverse_complement("TTTTGGGGCC") == "GGCCCCAAAA"


def test_load_library_with_header():
    library = load_library(str(LIB_CSV))
    assert library.sequence_to_id == {
        "ACGTACGTAC": "seq1",
        "TTTTGGGGCC": "seq2",
    }
    # Uniform 10bp library -> target length auto-detected.
    assert library.target_length == 10


def test_load_library_without_header(tmp_path):
    # No header row; first row is data. Second column is long and lacks "seq".
    csv_path = tmp_path / "nohdr.csv"
    csv_path.write_text("id1,ACGTACGTAC\nid2,TTTTGGGGCC\n")
    library = load_library(str(csv_path))
    assert library.sequence_to_id == {
        "ACGTACGTAC": "id1",
        "TTTTGGGGCC": "id2",
    }
    assert library.target_length == 10


def test_load_library_varying_lengths(tmp_path):
    csv_path = tmp_path / "varied.csv"
    csv_path.write_text("ID,Sequence\na,ACGT\nb,ACGTACGT\n")
    library = load_library(str(csv_path))
    # Mixed lengths -> no single auto-detectable target length.
    assert library.target_length is None


def test_resolve_trim_length():
    library = load_library(str(LIB_CSV))
    # Explicit value wins.
    assert resolve_trim_length(42, library) == 42
    # Otherwise falls back to the library's uniform length.
    assert resolve_trim_length(None, library) == 10


def test_process_reads_counts_and_qc():
    library = load_library(str(LIB_CSV))
    trim_length = resolve_trim_length(None, library)
    result = process_reads(
        str(SAMPLE_FASTQ),
        library,
        trim_start=0,
        trim_length=trim_length,
    )

    assert result.total_reads == 4
    # read3 is all low-quality bases -> dropped by QC.
    assert result.skipped_qc == 1
    assert result.passed_reads == 3
    # read1 -> seq1, read2 -> seq2, read4 no match.
    assert result.matched_reads == 2
    assert result.counts == {"seq1": 1, "seq2": 1}


def test_write_output_sorted_desc(tmp_path):
    out = tmp_path / "out.csv"
    write_output(str(out), {"seq1": 5, "seq2": 10, "seq3": 0})

    with open(out) as f:
        rows = list(csv.reader(f))

    assert rows[0] == ["Library_ID", "Count"]
    # Sorted by count descending.
    assert rows[1] == ["seq2", "10"]
    assert rows[2] == ["seq1", "5"]
    assert rows[3] == ["seq3", "0"]
