"""Sanity tests for tcrlibs.diagnose (fuller coverage in a later task)."""

from pathlib import Path

import pytest

from tcrlibs.counting import get_reverse_complement
from tcrlibs.diagnose import diagnose_library, diagnose_targets, format_report

DATA = Path(__file__).parent / "data"
SEED = "GAGGACCTGAACAAGGTGTTTCCTCCAGAGGTGGCCGTGTTC"


def test_forward_seed_on_full_lib():
    d = diagnose_library(str(DATA / "full_lib.csv"), SEED, read_length=len(SEED) + 8)
    assert d.total_rows == 2
    assert d.forward_hits == 2 and d.reverse_complement_hits == 0
    assert d.rows_one == 2 and d.rows_zero == 0
    assert d.orientation == "forward" and d.suggested_seed is None
    assert d.length_distribution == {8: 2}
    assert d.duplicates == {} and d.short_ids == []
    assert d.usable
    assert "Orientation: forward" in format_report(d)


def test_reverse_complement_seed_suggested(tmp_path):
    rc = get_reverse_complement(SEED)
    d = diagnose_library(str(DATA / "full_lib.csv"), rc.lower(), read_length=150)
    assert d.forward_hits == 0 and d.reverse_complement_hits == 2
    assert d.orientation == "reverse_complement"
    assert d.suggested_seed == SEED
    assert d.missing_ids == ["clone_1", "clone_2"]
    assert not d.usable
    assert list(tmp_path.iterdir()) == []  # nothing written


def test_multi_short_and_duplicates(tmp_path):
    csv = tmp_path / "lib.csv"
    csv.write_text("ID,Sequence\na,AAAAGGTTGG\nb,AAAAGG\nc,CGG\nd,TTTT\n")
    d = diagnose_library(str(csv), "gg", read_length=6)
    assert (d.rows_zero, d.rows_one, d.rows_multi) == (1, 2, 1)
    assert d.multi_ids == ["a"] and d.missing_ids == ["d"]
    assert d.short_ids == ["c"]
    assert d.duplicates == {"AAAA": ["a", "b"]}
    assert d.orientation == "forward"


def test_targets_report(tmp_path):
    csv = tmp_path / "t.csv"
    csv.write_text("ID,Sequence\n1,acgt \n2,ACGT\n3,AC\n")
    t = diagnose_targets(str(csv))
    assert t.total_rows == 3 and t.unique_sequences == 2
    assert t.length_distribution == {2: 1, 4: 2}
    assert t.duplicates == {"ACGT": ["1", "2"]}


@pytest.mark.parametrize(
    "seed,read_length", [("", 150), ("ACGT", 3)]
)
def test_invalid_args(seed, read_length):
    with pytest.raises(ValueError):
        diagnose_library(str(DATA / "full_lib.csv"), seed, read_length)


def test_missing_columns(tmp_path):
    csv = tmp_path / "bad.csv"
    csv.write_text("Name,Seq\nx,ACGT\n")
    with pytest.raises(ValueError):
        diagnose_library(str(csv), "AC", 10)
    with pytest.raises(ValueError):
        diagnose_targets(str(csv))


def test_absent_seed_not_found():
    d = diagnose_library(str(DATA / "full_lib.csv"), "NNNNNNNN", read_length=20)
    assert d.forward_hits == 0 and d.reverse_complement_hits == 0
    assert d.orientation == "not_found" and d.suggested_seed is None
    assert d.rows_zero == 2 and d.missing_ids == ["clone_1", "clone_2"]
    assert d.length_distribution == {} and d.duplicates == {}
    assert not d.usable
    report = format_report(d)
    assert "Orientation: not_found" in report and "Usable: no" in report


def test_tie_between_forward_and_rc_reports_forward(tmp_path):
    # ACGT is its own reverse complement, so forward and RC hits are equal.
    csv = tmp_path / "lib.csv"
    csv.write_text("ID,Sequence\nx,TTACGTAA\n")
    d = diagnose_library(str(csv), "acgt", read_length=6)
    assert d.forward_hits == d.reverse_complement_hits == 1
    assert d.orientation == "forward" and d.suggested_seed is None


def test_targets_csv_via_diagnose_library_and_no_writes(tmp_path):
    full = tmp_path / "full.csv"
    full.write_text("ID,Sequence\na,AAAAGGTTGG\nb,AAAAGG\nc,CGG\nd,TTTT\n")
    targets = tmp_path / "targets.csv"
    targets.write_text("ID,Sequence\nt1,AAAA\nt2, aaaa\nt3,AAAA \nt4,CC\n")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}

    d = diagnose_library(str(full), "GG", read_length=6, targets_csv=str(targets))

    after = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    assert after == before  # no files created or modified

    t = d.targets_report
    assert t is not None
    assert t.total_rows == 4 and t.unique_sequences == 2
    assert t.duplicates == {"AAAA": ["t1", "t2", "t3"]}
    assert t.length_distribution == {2: 1, 4: 3}

    # Full-library numbers alongside the targets report.
    assert d.length_distribution == {1: 1, 4: 2}
    assert d.keep_length == 4

    report = format_report(d)
    for snippet in [
        "Total rows: 4",
        "Rows containing seed: 3",
        "Keep length: 4bp",
        "0: 1",
        ">1: 1",
        "IDs: d",
        "IDs (first occurrence used): a",
        "Shorter than 4bp: 1",
        "IDs: c",
        "AAAA shared by 2 IDs: a, b",
        "Targets library",
        "Rows: 4",
        "Unique sequences: 2",
        "AAAA shared by 3 IDs: t1, t2, t3",
        "Usable: no",
    ]:
        assert snippet in report, snippet


def test_format_report_suggests_rc_seed():
    rc = get_reverse_complement(SEED)
    d = diagnose_library(str(DATA / "full_lib.csv"), rc, read_length=150)
    report = format_report(d)
    assert "Orientation: reverse_complement" in report
    assert f"Suggested seed: --seed {SEED}" in report
    assert "Rows containing reverse complement: 2" in report


def test_format_report_truncates_ids(tmp_path):
    csv = tmp_path / "lib.csv"
    rows = "\n".join(f"r{i},TTTT" for i in range(12))
    csv.write_text(f"ID,Sequence\n{rows}\n")
    d = diagnose_library(str(csv), "GG", read_length=6)
    assert len(d.missing_ids) == 12
    assert "(+2 more)" in format_report(d)


# --- CLI: tcrlibs diagnose-library ---------------------------------------

from tcrlibs.cli import main  # noqa: E402


def test_cli_usable_library_exits_zero(capsys):
    rc = main([
        "diagnose-library", "--input", str(DATA / "full_lib.csv"),
        "--seed", SEED, "--read-length", str(len(SEED) + 8),
    ])
    out = capsys.readouterr()
    assert rc == 0
    assert "Orientation: forward" in out.out and "Usable: yes" in out.out
    assert out.err == ""


def test_cli_reversed_seed_exits_one_with_suggestion(capsys):
    rc = main([
        "diagnose-library", "--input", str(DATA / "full_lib.csv"),
        "--seed", get_reverse_complement(SEED),
    ])
    out = capsys.readouterr()
    assert rc == 1
    assert "Orientation: reverse_complement" in out.out
    assert f"Suggested seed: --seed {SEED}" in out.out


@pytest.mark.parametrize("read_length", ["150", "2"])
def test_cli_errors_go_to_stderr(tmp_path, capsys, read_length):
    csv = tmp_path / "bad.csv"
    if read_length == "150":
        csv.write_text("Name,Seq\nx,ACGT\n")  # missing columns
    else:
        csv.write_text("ID,Sequence\nx,ACGT\n")  # read length < seed length
    rc = main(["diagnose-library", "--input", str(csv), "--seed", "ACG", "--read-length", read_length])
    out = capsys.readouterr()
    assert rc == 1
    assert out.err.startswith("Error: ")
    assert out.out == ""


def test_cli_targets_report_in_output(tmp_path, capsys):
    full = tmp_path / "full.csv"
    full.write_text("ID,Sequence\na,AAAAGG\nb,CCCCGG\n")
    targets = tmp_path / "targets.csv"
    targets.write_text("ID,Sequence\nt1,AAAA\nt2,aaaa\n")
    rc = main([
        "diagnose-library", "--input", str(full), "--seed", "GG",
        "--read-length", "6", "--targets", str(targets),
    ])
    out = capsys.readouterr()
    assert rc == 0
    assert "Targets library" in out.out
    assert "AAAA shared by 2 IDs: t1, t2" in out.out
