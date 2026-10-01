"""End-to-end smoke test for the `count` subcommand."""

import csv
from pathlib import Path

from tcrlibs.cli import main

DATA_DIR = Path(__file__).parent / "data"
LIB_CSV = DATA_DIR / "lib.csv"
SAMPLE_FASTQ = DATA_DIR / "sample.fastq"


def test_count_command_end_to_end(tmp_path):
    out = tmp_path / "out.csv"
    trimmed = tmp_path / "trimmed.csv"

    rc = main(
        [
            "count",
            "-r",
            str(SAMPLE_FASTQ),
            "-l",
            str(LIB_CSV),
            "-o",
            str(out),
            "--output-trimmed",
            str(trimmed),
        ]
    )
    assert rc == 0
    assert out.exists()

    with open(out) as f:
        rows = list(csv.reader(f))

    assert rows[0] == ["Library_ID", "Count"]
    counts = {row[0]: int(row[1]) for row in rows[1:]}
    assert counts == {"seq1": 1, "seq2": 1}

    # Trimmed output records every passed read (3; read3 dropped by QC).
    with open(trimmed) as f:
        trimmed_rows = list(csv.reader(f))
    assert trimmed_rows[0] == ["Read_ID", "Trimmed_Sequence", "Matched_Library_ID"]
    assert len(trimmed_rows) == 1 + 3
