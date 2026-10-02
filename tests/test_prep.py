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


# --- Pure helpers (compute_targets / find_duplicates / _format_ids) ---

from tcrlibs.prep import (  # noqa: E402
    SeedNotFoundError,
    _format_ids,
    compute_targets,
    find_duplicates,
)


def test_compute_targets_first_occurrence_and_counts():
    recs = compute_targets(
        ["a", "b", "c", "d"],
        ["AAACCCSEEDTT", "GGSEEDXXSEED", "NOPE", "TTTTTTTTSEED"],
        "SEED",
        4,
    )
    a, b, c, d = recs
    assert (a.occurrences, a.prefix_length, a.target, a.is_short) == (1, 6, "ACCC", False)
    assert (b.occurrences, b.prefix_length, b.target, b.is_short) == (2, 2, "GG", True)
    assert (c.occurrences, c.prefix_length, c.target, c.is_short) == (0, 0, None, False)
    assert d.target == "TTTT"


def test_compute_targets_zero_keep_length_gives_empty_target():
    (rec,) = compute_targets(["a"], ["ACGTSEED"], "SEED", 0)
    assert rec.target == ""
    assert not rec.is_short


def test_find_duplicates_keeps_input_order():
    dups = find_duplicates(["a", "b", "c", "d", "e"], ["X", "Y", "X", "Z", "X"])
    assert dups == {"X": ["a", "c", "e"]}


def test_format_ids_truncates():
    assert _format_ids(["a", "b"]) == "a, b"
    ids = [f"id{i}" for i in range(12)]
    out = _format_ids(ids)
    assert out.startswith("id0, id1") and out.endswith("id9, ... (+2 more)")


def test_seed_not_found_error_is_value_error():
    err = SeedNotFoundError("msg", 3, ["a", "b", "c"], "ACGT")
    assert isinstance(err, ValueError)
    assert (err.missing_count, err.missing_ids, err.suggested_seed) == (3, ["a", "b", "c"], "ACGT")


# --- trim_targets validation / warnings (minimal; fuller coverage in task 1.5) ---


def _write_lib(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Sequence"])
        w.writerows(rows)
    return str(path)


def test_trim_targets_valid_input_is_silent(tmp_path, capsys):
    trim_targets(str(FULL_LIB), str(tmp_path / "o.csv"), seed=SEED, read_length=50)
    assert capsys.readouterr().err == ""


def test_trim_targets_missing_seed_raises_with_rc_hint(tmp_path):
    from tcrlibs.counting import get_reverse_complement

    rc = get_reverse_complement(SEED)
    lib = _write_lib(tmp_path / "lib.csv", [["a", "ACGTACGT" + SEED], ["b", "TTTT" + rc]])
    out = tmp_path / "o.csv"
    with pytest.raises(SeedNotFoundError) as exc:
        trim_targets(lib, str(out), seed=SEED, read_length=50)
    assert exc.value.missing_count == 1
    assert exc.value.suggested_seed == rc
    assert f"Try --seed {rc}" in str(exc.value)
    assert "--skip-missing-seed" in str(exc.value)
    assert not out.exists()


def test_trim_targets_skip_missing_seed_drops_and_warns(tmp_path, capsys):
    lib = _write_lib(
        tmp_path / "lib.csv",
        [["a", "ACGTACGT" + SEED], ["b", "NOSEEDHERE"], ["c", "GG" + SEED + SEED]],
    )
    out = tmp_path / "o.csv"
    n = trim_targets(lib, str(out), seed=SEED, read_length=50, skip_missing_seed=True)
    assert n == 2
    assert _read_rows(out) == [["ID", "Sequence"], ["a", "ACGTACGT"], ["c", "GG"]]
    err = capsys.readouterr().err
    assert "Skipped 1 rows lacking seed: b" in err
    assert "first occurrence used: c" in err
    assert "shorter than 8bp: c" in err


def test_trim_targets_warns_on_duplicates_and_empty_seed(tmp_path, capsys):
    lib = _write_lib(tmp_path / "lib.csv", [["a", "ACGTACGT" + SEED], ["b", "ACGTACGT" + SEED]])
    trim_targets(lib, str(tmp_path / "o.csv"), seed=SEED, read_length=50)
    assert "ACGTACGT shared by a, b" in capsys.readouterr().err
    with pytest.raises(ValueError):
        trim_targets(lib, str(tmp_path / "o2.csv"), seed="", read_length=50)


# --- trim_targets: fuller coverage (task 1.5) ---


def test_trim_targets_rc_seed_input_suggests_forward_seed(tmp_path):
    """Feeding the reverse-complemented seed suggests the forward seed."""
    from tcrlibs.counting import get_reverse_complement

    rc = get_reverse_complement(SEED)
    out = tmp_path / "o.csv"
    with pytest.raises(SeedNotFoundError) as exc:
        trim_targets(str(FULL_LIB), str(out), seed=rc, read_length=50)
    err = exc.value
    assert err.missing_count == 2
    assert err.missing_ids == ["clone_1", "clone_2"]
    assert err.suggested_seed == SEED
    msg = str(err)
    assert "Seed not found in 2 of 2 rows" in msg
    assert "matches 2 of these rows" in msg
    assert f"Try --seed {SEED}" in msg
    assert not out.exists()


def test_trim_targets_missing_seed_without_rc_match_has_no_suggestion(tmp_path):
    lib = _write_lib(tmp_path / "lib.csv", [["a", "ACGTACGT" + SEED], ["b", "ACGTNOSEED"]])
    out = tmp_path / "o.csv"
    with pytest.raises(SeedNotFoundError) as exc:
        trim_targets(lib, str(out), seed=SEED, read_length=50)
    assert exc.value.suggested_seed is None
    assert exc.value.missing_ids == ["b"]
    assert "Try --seed" not in str(exc.value)
    assert "reverse complement" not in str(exc.value)
    assert not out.exists()


def test_trim_targets_many_missing_ids_are_truncated(tmp_path):
    rows = [[f"m{i}", "ACGTACGT"] for i in range(12)]
    lib = _write_lib(tmp_path / "lib.csv", rows)
    with pytest.raises(SeedNotFoundError) as exc:
        trim_targets(lib, str(tmp_path / "o.csv"), seed=SEED, read_length=50)
    err = exc.value
    assert err.missing_count == 12
    assert err.missing_ids == [f"m{i}" for i in range(10)]
    msg = str(err)
    assert "Seed not found in 12 of 12 rows" in msg
    assert "m9, ... (+2 more)" in msg
    assert "m10" not in msg and "m11" not in msg


def test_trim_targets_skip_missing_seed_rc_hint_and_truncation(tmp_path, capsys):
    from tcrlibs.counting import get_reverse_complement

    rc = get_reverse_complement(SEED)
    rows = [["keep", "ACGTACGT" + SEED]] + [[f"m{i}", "TT" + rc] for i in range(12)]
    lib = _write_lib(tmp_path / "lib.csv", rows)
    out = tmp_path / "o.csv"
    n = trim_targets(lib, str(out), seed=SEED, read_length=50, skip_missing_seed=True)
    assert n == 1
    assert _read_rows(out) == [["ID", "Sequence"], ["keep", "ACGTACGT"]]
    err = capsys.readouterr().err
    assert "Skipped 12 rows lacking seed" in err
    assert "m9, ... (+2 more)" in err
    assert "matches 12 of these rows" in err
    assert f"Try --seed {rc}" in err


def test_trim_targets_skip_missing_seed_preserves_input_order(tmp_path):
    lib = _write_lib(
        tmp_path / "lib.csv",
        [["z", "AAAAAAAA" + SEED], ["x", "NOPE"], ["y", "CCCCCCCC" + SEED]],
    )
    out = tmp_path / "o.csv"
    n = trim_targets(lib, str(out), seed=SEED, read_length=50, skip_missing_seed=True)
    assert n == 2
    assert _read_rows(out)[1:] == [["z", "AAAAAAAA"], ["y", "CCCCCCCC"]]


def test_trim_targets_multi_short_duplicate_warnings_list_right_ids(tmp_path, capsys):
    lib = _write_lib(
        tmp_path / "lib.csv",
        [
            ["ok", "GGGGACGTACGT" + SEED],  # target ACGTACGT, clean
            ["multi", "TTTTTTTT" + SEED + "AC" + SEED],  # 2 occurrences
            ["short", "CAT" + SEED],  # 3bp < 8bp
            ["dup1", "CCCCCCCC" + SEED],
            ["dup2", "AACCCCCCCC" + SEED],  # same 8bp target as dup1
            ["dupA", "GATTACAA" + SEED],
            ["dupB", "GATTACAA" + SEED],
        ],
    )
    out = tmp_path / "o.csv"
    n = trim_targets(lib, str(out), seed=SEED, read_length=50)
    assert n == 7
    written = {r[0]: r[1] for r in _read_rows(out)[1:]}
    assert written["multi"] == "TTTTTTTT"  # first occurrence used
    assert written["short"] == "CAT"

    err = capsys.readouterr().err
    assert "seed occurs more than once in 1 rows; first occurrence used: multi" in err
    assert "1 targets shorter than 8bp: short" in err
    assert "duplicate target CCCCCCCC shared by dup1, dup2" in err
    assert "duplicate target GATTACAA shared by dupA, dupB" in err
    assert "ok" not in err


def test_trim_targets_short_warning_truncates_ids(tmp_path, capsys):
    rows = [[f"s{i:02d}", "ACGT"[: (i % 4) + 1] + "G" * (i // 4) + SEED] for i in range(12)]
    lib = _write_lib(tmp_path / "lib.csv", rows)
    trim_targets(lib, str(tmp_path / "o.csv"), seed=SEED, read_length=50)
    err = capsys.readouterr().err
    assert "12 targets shorter than 8bp" in err
    assert "s09, ... (+2 more)" in err


def test_trim_targets_empty_seed_raises_value_error(tmp_path):
    out = tmp_path / "o.csv"
    with pytest.raises(ValueError):
        trim_targets(str(FULL_LIB), str(out), seed="", read_length=50)
    assert not out.exists()


# --- CLI: prep-library missing-seed handling (task 2.2) ---


def test_prep_command_reversed_seed_errors_with_suggestion(tmp_path, capsys):
    """Requirement 2.5: missing seed -> rc 1, 'Error:' + suggested seed on stderr."""
    from tcrlibs.counting import get_reverse_complement

    rc_seed = get_reverse_complement(SEED)
    out = tmp_path / "o.csv"
    rc = main(
        ["prep-library", "--input", str(FULL_LIB), "--output", str(out),
         "--seed", rc_seed, "--read-length", "50"]
    )
    assert rc == 1
    err = capsys.readouterr().err
    assert "Error:" in err
    assert f"Try --seed {SEED}" in err
    assert not out.exists()


def test_prep_command_skip_missing_seed_writes_seeded_rows(tmp_path, capsys):
    """Requirement 2.7: --skip-missing-seed on a mixed library -> rc 0, seeded rows only."""
    lib = _write_lib(
        tmp_path / "lib.csv",
        [["a", "ACGTACGT" + SEED], ["b", "NOSEEDHERE"], ["c", "TTGGCCAA" + SEED]],
    )
    out = tmp_path / "o.csv"
    rc = main(
        ["prep-library", "--input", lib, "--output", str(out),
         "--seed", SEED, "--read-length", "50", "--skip-missing-seed"]
    )
    assert rc == 0
    assert _read_rows(out) == [["ID", "Sequence"], ["a", "ACGTACGT"], ["c", "TTGGCCAA"]]
    assert "Skipped 1 rows lacking seed: b" in capsys.readouterr().err
