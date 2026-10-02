"""Read-only diagnostics for a full library (and optionally a trimmed one).

``diagnose_library`` reports how the seed matches each row in both
orientations, how often it occurs, the resulting target lengths, and
duplicate targets. It shares the trimming math with ``prep.trim_targets``
(via ``prep.compute_targets``) so the two cannot disagree, and it never
writes any files.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from .counting import get_reverse_complement
from .prep import _format_ids, compute_targets, find_duplicates


@dataclass
class TargetsReport:
    """Summary of an already-trimmed targets CSV."""

    total_rows: int
    length_distribution: dict  # length -> row count
    duplicates: dict  # sequence -> IDs (size >= 2)
    unique_sequences: int  # what load_library would keep


@dataclass
class LibraryDiagnosis:
    """Everything ``diagnose-library`` reports about a full library."""

    total_rows: int
    seed: str
    reverse_complement_seed: str
    keep_length: int
    forward_hits: int
    reverse_complement_hits: int
    rows_zero: int
    rows_one: int
    rows_multi: int
    missing_ids: list = field(default_factory=list)
    multi_ids: list = field(default_factory=list)
    short_ids: list = field(default_factory=list)
    length_distribution: dict = field(default_factory=dict)
    duplicates: dict = field(default_factory=dict)
    suggested_seed: Optional[str] = None
    orientation: str = "not_found"  # "forward" | "reverse_complement" | "not_found"
    targets_report: Optional[TargetsReport] = None

    @property
    def usable(self) -> bool:
        """True when the library has rows and every row contains the seed."""
        return self.total_rows > 0 and self.rows_zero == 0


def _read_id_sequence_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "Sequence" not in df.columns or "ID" not in df.columns:
        raise ValueError(f"Input CSV must contain 'ID' and 'Sequence' columns: {path}")
    return df


def diagnose_targets(targets_csv: str) -> TargetsReport:
    """Report length distribution, duplicates and unique count of a targets CSV.

    Sequences are stripped and upper-cased to mirror ``load_library``.
    """
    df = _read_id_sequence_csv(targets_csv)
    ids = [str(i) for i in df["ID"]]
    seqs = [str(s).strip().upper() for s in df["Sequence"].fillna("")]
    return TargetsReport(
        total_rows=len(ids),
        length_distribution=dict(sorted(Counter(len(s) for s in seqs).items())),
        duplicates=find_duplicates(ids, seqs),
        unique_sequences=len(set(seqs)),
    )


def diagnose_library(
    input_csv: str,
    seed: str,
    read_length: int = 150,
    targets_csv: Optional[str] = None,
) -> LibraryDiagnosis:
    """Diagnose how ``seed`` matches a full library. Writes no files."""
    seed_seq = seed.upper()
    if not seed_seq:
        raise ValueError("Seed must be non-empty.")
    keep_length = read_length - len(seed_seq)
    if keep_length < 0:
        raise ValueError(
            f"read_length ({read_length}) is shorter than the seed ({len(seed_seq)}bp)."
        )

    df = _read_id_sequence_csv(input_csv)
    ids = [str(i) for i in df["ID"]]
    seqs = [str(s).upper() for s in df["Sequence"].fillna("")]
    rc_seed = get_reverse_complement(seed_seq)

    records = compute_targets(ids, seqs, seed_seq, keep_length)
    seeded = [r for r in records if r.occurrences > 0]

    forward_hits = len(seeded)
    rc_hits = sum(1 for s in seqs if rc_seed in s)

    if rc_hits > forward_hits:
        orientation, suggested = "reverse_complement", rc_seed
    elif forward_hits > 0:
        orientation, suggested = "forward", None
    else:
        orientation, suggested = "not_found", None

    return LibraryDiagnosis(
        total_rows=len(records),
        seed=seed_seq,
        reverse_complement_seed=rc_seed,
        keep_length=keep_length,
        forward_hits=forward_hits,
        reverse_complement_hits=rc_hits,
        rows_zero=sum(1 for r in records if r.occurrences == 0),
        rows_one=sum(1 for r in records if r.occurrences == 1),
        rows_multi=sum(1 for r in records if r.occurrences > 1),
        missing_ids=[r.id for r in records if r.occurrences == 0],
        multi_ids=[r.id for r in records if r.occurrences > 1],
        short_ids=[r.id for r in seeded if r.is_short],
        length_distribution=dict(
            sorted(Counter(len(r.target) for r in seeded).items())
        ),
        duplicates=find_duplicates([r.id for r in seeded], [r.target for r in seeded]),
        suggested_seed=suggested,
        orientation=orientation,
        targets_report=diagnose_targets(targets_csv) if targets_csv else None,
    )


def _format_distribution(dist: dict) -> str:
    if not dist:
        return "(none)"
    return ", ".join(f"{length}bp: {n}" for length, n in dist.items())


def _format_duplicates(dups: dict, max_ids: int) -> list:
    lines = [f"  Duplicate groups: {len(dups)}"]
    for seq, members in dups.items():
        lines.append(f"    {seq} shared by {len(members)} IDs: {_format_ids(members, max_ids)}")
    return lines


def format_report(d: LibraryDiagnosis, max_ids: int = 10) -> str:
    """Render a ``LibraryDiagnosis`` as a human-readable report."""
    lines = [
        "Library diagnosis",
        f"  Total rows: {d.total_rows}",
        f"  Seed: {d.seed}",
        f"  Reverse complement: {d.reverse_complement_seed}",
        f"  Keep length: {d.keep_length}bp",
        f"  Rows containing seed: {d.forward_hits}",
        f"  Rows containing reverse complement: {d.reverse_complement_hits}",
        f"  Orientation: {d.orientation}",
    ]
    if d.suggested_seed:
        lines.append(f"  Suggested seed: --seed {d.suggested_seed}")
    lines.append("Seed occurrences per row")
    lines.append(f"  0: {d.rows_zero}")
    if d.missing_ids:
        lines.append(f"    IDs: {_format_ids(d.missing_ids, max_ids)}")
    lines.append(f"  1: {d.rows_one}")
    lines.append(f"  >1: {d.rows_multi}")
    if d.multi_ids:
        lines.append(f"    IDs (first occurrence used): {_format_ids(d.multi_ids, max_ids)}")
    lines.append("Targets (rows containing seed)")
    lines.append(f"  Length distribution: {_format_distribution(d.length_distribution)}")
    lines.append(f"  Shorter than {d.keep_length}bp: {len(d.short_ids)}")
    if d.short_ids:
        lines.append(f"    IDs: {_format_ids(d.short_ids, max_ids)}")
    lines.extend(_format_duplicates(d.duplicates, max_ids))

    t = d.targets_report
    if t is not None:
        lines.append("Targets library")
        lines.append(f"  Rows: {t.total_rows}")
        lines.append(f"  Unique sequences: {t.unique_sequences}")
        lines.append(f"  Length distribution: {_format_distribution(t.length_distribution)}")
        lines.extend(_format_duplicates(t.duplicates, max_ids))

    lines.append(f"Usable: {'yes' if d.usable else 'no'}")
    return "\n".join(lines)
