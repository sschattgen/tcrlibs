"""Smoke tests for prep-library target trimming."""

import csv
from pathlib import Path

import pytest

from tcrlibs.cli import main
from tcrlibs.prep import trim_targets

DATA_DIR = Path(__file__).parent / "data"
FULL_LIB = DATA_DIR / "full_lib.csv"
SEED = "GAGGACCTGAACAAGGTGTTTCCTCCAGAGGTGGCCGTGTTC"  # 42 bp


def _read_rows(path):
    with open(path) as f:
        return list(csv.reader(f))


def test_trim_targets_keeps_bases_before_seed(tmp_path):
    out = tmp_path / "clean.csv"
    # read_length 50, seed 42 bp -> keep 8 bases immediately before the seed.
    n = trim_targets(str(FULL_LIB), str(out), seed=SEED, read_length=50)
    assert n == 2

    rows = _read_rows(out)
    assert rows[0] == ["ID", "Sequence"]
    result = {row[0]: row[1] for row in rows[1:]}
    assert result == {"clone_1": "ACGTACGT", "clone_2": "TTGGCCAA"}
    # Every target is keep_length (8) bases.
    assert all(len(seq) == 8 for seq in result.values())


def test_trim_targets_case_insensitive_seed(tmp_path):
    out = tmp_path / "clean.csv"
    n = trim_targets(str(FULL_LIB), str(out), seed=SEED.lower(), read_length=50)
    assert n == 2
    rows = _read_rows(out)
    result = {row[0]: row[1] for row in rows[1:]}
    assert result == {"clone_1": "ACGTACGT", "clone_2": "TTGGCCAA"}


def test_trim_targets_rejects_short_read_length(tmp_path):
    out = tmp_path / "clean.csv"
    with pytest.raises(ValueError):
        trim_targets(str(FULL_LIB), str(out), seed=SEED, read_length=10)


def test_prep_command_end_to_end(tmp_path):
    out = tmp_path / "clean.csv"
    rc = main(
        [
            "prep-library",
            "--input",
            str(FULL_LIB),
            "--output",
            str(out),
            "--seed",
            SEED,
            "--read-length",
            "50",
        ]
    )
    assert rc == 0
    assert out.exists()
    rows = _read_rows(out)
    assert rows[0] == ["ID", "Sequence"]
