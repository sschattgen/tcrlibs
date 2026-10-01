"""Prepare a trimmed target library from full-length sequences.

Ported from the standalone ``generate_lib_targets.py`` notebook script. For
each sequence, we locate a constant "seed" region and keep the ``keep_length``
bases immediately preceding it, where ``keep_length = read_length - len(seed)``.
This yields targets matched to the sequencer's read length.

The input CSV must contain ``ID`` and ``Sequence`` columns (extra columns such
as ``Subpool`` are ignored). The output CSV has ``ID`` and ``Sequence``.
"""

import pandas as pd


def trim_targets(
    input_csv: str,
    output_csv: str,
    seed: str,
    read_length: int = 150,
) -> int:
    """Trim full sequences to read-length targets and write ``ID,Sequence``.

    Args:
        input_csv: path to the full library CSV (needs ID and Sequence columns).
        output_csv: path to write the trimmed library.
        seed: the constant region to split on (case-insensitive; upper-cased).
        read_length: the sequencer read length; the kept region is
            ``read_length - len(seed)`` bases preceding the seed.

    Returns:
        The number of target rows written.
    """
    seed_seq = seed.upper()
    keep_length = read_length - len(seed_seq)
    if keep_length < 0:
        raise ValueError(
            f"read_length ({read_length}) is shorter than the seed ({len(seed_seq)}bp)."
        )

    df = pd.read_csv(input_csv)
    if "Sequence" not in df.columns or "ID" not in df.columns:
        raise ValueError("Input CSV must contain 'ID' and 'Sequence' columns.")

    # For each sequence, take the portion before the seed and keep the last
    # keep_length bases of it. If the seed is absent, split returns the whole
    # string as the first element, matching the original script's behavior.
    split_seqs = df["Sequence"].str.upper().str.split(seed_seq)
    target_seqs = [parts[0][-keep_length:] for parts in split_seqs]

    out = df[["ID"]].copy()
    out["Sequence"] = target_seqs
    out.to_csv(output_csv, index=False)

    return len(out)
