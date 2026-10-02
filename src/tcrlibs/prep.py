"""Prepare a trimmed target library from full-length sequences.

Ported from the standalone ``generate_lib_targets.py`` notebook script. For
each sequence, we locate a constant "seed" region and keep the ``keep_length``
bases immediately preceding it, where ``keep_length = read_length - len(seed)``.
This yields targets matched to the sequencer's read length.

The input CSV must contain ``ID`` and ``Sequence`` columns (extra columns such
as ``Subpool`` are ignored). The output CSV has ``ID`` and ``Sequence``.
"""

import sys
from dataclasses import dataclass
from typing import Optional, Sequence

import pandas as pd

from .counting import get_reverse_complement


class SeedNotFoundError(ValueError):
    """Raised when library rows lack the seed.

    Subclasses ``ValueError`` so existing handlers keep working.
    """

    def __init__(
        self,
        message: str,
        missing_count: int,
        missing_ids: list,
        suggested_seed: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.missing_count = missing_count
        self.missing_ids = list(missing_ids)
        self.suggested_seed = suggested_seed


@dataclass(frozen=True)
class TargetRecord:
    """Result of trimming one library sequence against the seed."""

    id: str
    occurrences: int  # non-overlapping count of the seed in the sequence
    prefix_length: int  # bases before the FIRST seed occurrence (0 if absent)
    target: Optional[str]  # None if the seed is absent
    keep_length: int  # requested target length

    @property
    def is_short(self) -> bool:
        """True when a target exists but is shorter than ``keep_length``."""
        return self.target is not None and len(self.target) < self.keep_length


def compute_targets(
    ids: Sequence[str],
    sequences: Sequence[str],
    seed: str,
    keep_length: int,
) -> list:
    """Compute a ``TargetRecord`` per sequence, preserving input order.

    The target is the last ``keep_length`` bases before the first seed
    occurrence (or the whole prefix when it is shorter). ``keep_length == 0``
    yields an empty target.
    """
    if len(ids) != len(sequences):
        raise ValueError("ids and sequences must have the same length.")
    if not seed:
        raise ValueError("Seed must be non-empty.")
    if keep_length < 0:
        raise ValueError("keep_length must be >= 0.")

    records = []
    for id_, seq in zip(ids, sequences):
        occurrences = seq.count(seed)
        if occurrences == 0:
            records.append(TargetRecord(str(id_), 0, 0, None, keep_length))
            continue
        prefix = seq.split(seed, 1)[0]
        target = prefix[max(0, len(prefix) - keep_length):]
        records.append(
            TargetRecord(str(id_), occurrences, len(prefix), target, keep_length)
        )
    return records


def find_duplicates(ids: Sequence[str], sequences: Sequence[str]) -> dict:
    """Return ``{sequence: [ids in input order]}`` for sequences seen 2+ times."""
    groups: dict = {}
    for id_, seq in zip(ids, sequences):
        groups.setdefault(seq, []).append(str(id_))
    return {seq: members for seq, members in groups.items() if len(members) >= 2}


def _format_ids(ids: Sequence[str], limit: int = 10) -> str:
    """Join up to ``limit`` IDs, appending ``... (+N more)`` when truncated."""
    ids = [str(i) for i in ids]
    shown = ", ".join(ids[:limit])
    extra = len(ids) - limit
    if extra > 0:
        return f"{shown}, ... (+{extra} more)"
    return shown


def _warn(message: str) -> None:
    """Write a ``Warning:``-prefixed message to stderr."""
    print(f"Warning: {message}", file=sys.stderr)


def trim_targets(
    input_csv: str,
    output_csv: str,
    seed: str,
    read_length: int = 150,
    skip_missing_seed: bool = False,
) -> int:
    """Trim full sequences to read-length targets and write ``ID,Sequence``.

    Args:
        input_csv: path to the full library CSV (needs ID and Sequence columns).
        output_csv: path to write the trimmed library.
        seed: the constant region to split on (case-insensitive; upper-cased).
        read_length: the sequencer read length; the kept region is
            ``read_length - len(seed)`` bases preceding the seed.
        skip_missing_seed: drop rows lacking the seed (with a warning) instead
            of raising ``SeedNotFoundError``.

    Returns:
        The number of target rows written.

    Raises:
        ValueError: empty seed, ``read_length`` shorter than the seed, or
            missing columns.
        SeedNotFoundError: rows lack the seed and ``skip_missing_seed`` is
            False. No output file is written in that case.
    """
    seed_seq = seed.upper()
    if not seed_seq:
        raise ValueError("Seed must be non-empty.")
    keep_length = read_length - len(seed_seq)
    if keep_length < 0:
        raise ValueError(
            f"read_length ({read_length}) is shorter than the seed ({len(seed_seq)}bp)."
        )

    df = pd.read_csv(input_csv)
    if "Sequence" not in df.columns or "ID" not in df.columns:
        raise ValueError("Input CSV must contain 'ID' and 'Sequence' columns.")

    ids = [str(i) for i in df["ID"]]
    sequences = [str(s).upper() for s in df["Sequence"]]
    records = compute_targets(ids, sequences, seed_seq, keep_length)

    missing_idx = [i for i, r in enumerate(records) if r.occurrences == 0]
    if missing_idx:
        missing_ids = [records[i].id for i in missing_idx]
        rc_seed = get_reverse_complement(seed_seq)
        rc_hits = sum(1 for i in missing_idx if rc_seed in sequences[i])
        rc_hint = ""
        if rc_hits > 0:
            rc_hint = (
                f"The reverse complement of the seed matches {rc_hits} of these "
                "rows; the seed may be in the wrong orientation.\n"
                f"Try --seed {rc_seed}"
            )
        if not skip_missing_seed:
            message = (
                f"Seed not found in {len(missing_idx)} of {len(records)} rows "
                f"(e.g. {_format_ids(missing_ids)})."
            )
            if rc_hint:
                message += "\n" + rc_hint
            message += "\n(or pass --skip-missing-seed to drop rows without the seed)"
            raise SeedNotFoundError(
                message,
                missing_count=len(missing_idx),
                missing_ids=missing_ids[:10],
                suggested_seed=rc_seed if rc_hits > 0 else None,
            )
        warning = (
            f"Skipped {len(missing_idx)} rows lacking seed: {_format_ids(missing_ids)}"
        )
        if rc_hint:
            warning += "\n" + rc_hint
        _warn(warning)

    kept = [r for r in records if r.occurrences > 0]

    multi_ids = [r.id for r in kept if r.occurrences > 1]
    if multi_ids:
        _warn(
            f"seed occurs more than once in {len(multi_ids)} rows; "
            f"first occurrence used: {_format_ids(multi_ids)}"
        )

    short_ids = [r.id for r in kept if r.is_short]
    if short_ids:
        _warn(
            f"{len(short_ids)} targets shorter than {keep_length}bp: "
            f"{_format_ids(short_ids)}"
        )

    duplicates = find_duplicates([r.id for r in kept], [r.target for r in kept])
    for target, dup_ids in duplicates.items():
        _warn(f"duplicate target {target} shared by {_format_ids(dup_ids)}")

    # Written only after validation so a failed run leaves no partial CSV.
    kept_mask = [r.occurrences > 0 for r in records]
    out = df.loc[kept_mask, ["ID"]].copy()
    out["Sequence"] = [r.target for r in kept]
    out.to_csv(output_csv, index=False)

    return len(out)
