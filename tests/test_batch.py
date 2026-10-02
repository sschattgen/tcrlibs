"""Smoke tests for the batch runner and the `batch` subcommand."""

import csv
import gzip
from pathlib import Path

from tcrlibs.batch import discover_samples, run_batch
from tcrlibs.cli import main

DATA_DIR = Path(__file__).parent / "data"
LIB_CSV = DATA_DIR / "lib.csv"

# Two tiny samples. SampleA matches seq1 once; SampleB matches seq2 twice.
SAMPLE_A_READS = [
    ("read1", "ACGTACGTAC", "IIIIIIIIII"),  # seq1
    ("read2", "AAAAAAAAAA", "IIIIIIIIII"),  # no match
]
SAMPLE_B_READS = [
    ("read1", "TTTTGGGGCC", "IIIIIIIIII"),  # seq2
    ("read2", "TTTTGGGGCC", "IIIIIIIIII"),  # seq2
]


def _write_fastq_gz(path: Path, reads):
    with gzip.open(path, "wt") as f:
        for read_id, seq, qual in reads:
            f.write(f"@{read_id}\n{seq}\n+\n{qual}\n")


def _make_fastq_dir(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    _write_fastq_gz(fastq_dir / "SampleA_S1_R1_001.fastq.gz", SAMPLE_A_READS)
    _write_fastq_gz(fastq_dir / "SampleB_S2_R1_001.fastq.gz", SAMPLE_B_READS)
    return fastq_dir


def _read_counts(path: Path):
    with open(path) as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["Library_ID", "Count"]
    return {row[0]: int(row[1]) for row in rows[1:]}


def test_discover_samples(tmp_path):
    fastq_dir = _make_fastq_dir(tmp_path)
    samples = discover_samples(str(fastq_dir))
    names = [name for _, name in samples]
    # Sample name is everything before "_R1_001" (so the "_S1" lane tag stays).
    assert names == ["SampleA_S1", "SampleB_S2"]


def test_run_batch(tmp_path):
    fastq_dir = _make_fastq_dir(tmp_path)
    out_dir = tmp_path / "outs"

    summary = run_batch(
        input_dir=str(fastq_dir),
        output_dir=str(out_dir),
        library_path=str(LIB_CSV),
        trim_start=0,
    )

    assert summary.processed == 2
    assert summary.failed == 0

    counts_a = _read_counts(out_dir / "SampleA_S1_counts.csv")
    counts_b = _read_counts(out_dir / "SampleB_S2_counts.csv")
    assert counts_a == {"seq1": 1, "seq2": 0}
    assert counts_b == {"seq1": 0, "seq2": 2}

    # Trimmed CSVs written per sample.
    assert (out_dir / "SampleA_S1_trimmed.csv").exists()
    assert (out_dir / "SampleB_S2_trimmed.csv").exists()


def test_batch_command_end_to_end(tmp_path):
    fastq_dir = _make_fastq_dir(tmp_path)
    out_dir = tmp_path / "outs"

    rc = main(
        [
            "batch",
            "-i",
            str(fastq_dir),
            "-o",
            str(out_dir),
            "-l",
            str(LIB_CSV),
            "--trim-start",
            "0",
        ]
    )
    assert rc == 0
    assert (out_dir / "SampleA_S1_counts.csv").exists()
    assert (out_dir / "SampleB_S2_counts.csv").exists()

def test_discover_samples_custom_pattern(tmp_path):
    fastq_dir = _make_fastq_dir(tmp_path)
    # Write a differently-named file that only matches a custom pattern.
    _write_fastq_gz(fastq_dir / "SampleC.R1.fastq.gz", SAMPLE_A_READS)

    samples = discover_samples(str(fastq_dir), pattern="*.R1.fastq.gz")
    names = [name for _, name in samples]
    assert names == ["SampleC"]

def test_run_batch_custom_pattern(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    _write_fastq_gz(fastq_dir / "SampleC.R1.fastq.gz", SAMPLE_B_READS)
    out_dir = tmp_path / "outs"

    summary = run_batch(
        input_dir=str(fastq_dir),
        output_dir=str(out_dir),
        library_path=str(LIB_CSV),
        pattern="*.R1.fastq.gz",
        trim_start=0,
    )
    assert summary.processed == 1
    assert (out_dir / "SampleC_counts.csv").exists()


def test_batch_command_custom_pattern(tmp_path):
    """CLI `--pattern` is forwarded to discovery and sample naming.

    **Validates: Requirements 6.1, 6.2**
    """
    import shutil

    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    # gzip the shared sample.fastq fixture under a non-default name.
    with open(DATA_DIR / "sample.fastq", "rb") as src, gzip.open(
        fastq_dir / "SampleC.R1.fastq.gz", "wb"
    ) as dst:
        shutil.copyfileobj(src, dst)
    out_dir = tmp_path / "outs"

    rc = main(
        [
            "batch",
            "-i",
            str(fastq_dir),
            "-o",
            str(out_dir),
            "-l",
            str(LIB_CSV),
            "--pattern",
            "*.R1.fastq.gz",
            "--trim-start",
            "0",
        ]
    )
    assert rc == 0
    # Suffix ".R1.fastq.gz" is stripped, leaving "SampleC".
    counts = _read_counts(out_dir / "SampleC_counts.csv")
    # read3 is dropped by QC, read4 matches nothing.
    assert counts == {"seq1": 1, "seq2": 1}
