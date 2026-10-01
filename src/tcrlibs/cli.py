"""Unified command-line interface for tcrlibs."""

import argparse
import sys

from . import __version__
from . import batch as batch_mod
from . import counting
from . import prep as prep_mod


def _add_count_parser(subparsers) -> None:
    p = subparsers.add_parser(
        "count",
        help="Count reads from one FASTQ against a library.",
        description="Count sequencing reads mapping to a reference library.",
    )
    p.add_argument("-r", "--reads", required=True, help="Path to input FASTQ or FASTQ.GZ file")
    p.add_argument("-l", "--library", required=True, help="Path to reference library CSV (Format: ID, Sequence)")
    p.add_argument("-o", "--output", required=True, help="Path to output CSV file")
    p.add_argument("--output-trimmed", help="Optional path to save all trimmed reads and their match status (CSV)")
    p.add_argument("--trim-start", type=int, default=0, help="Number of bases to trim from start of read (default: 0)")
    p.add_argument(
        "--trim-length",
        type=int,
        default=None,
        help="Length of sequence to keep after trimming (default: match library length)",
    )
    p.add_argument("--revcomp", action="store_true", help="Search for reverse complement of library sequences")
    p.add_argument(
        "--min-base-q",
        type=int,
        default=28,
        help='Minimum Phred quality score for a base to be considered "good" (default: 28)',
    )
    p.add_argument(
        "--max-bad-freq",
        type=float,
        default=10.0,
        help='Maximum percentage of "bad" bases allowed in a read before dropping it (default: 10%%)',
    )
    p.set_defaults(func=_run_count)


def _run_count(args) -> int:
    print(f"Loading library from {args.library}...")
    library = counting.load_library(args.library)
    print(f"Loaded {len(library.sequence_to_id)} unique sequences.")

    cut_len = counting.resolve_trim_length(args.trim_length, library)
    if args.trim_length is None and library.target_length is not None:
        print(f"Auto-detected trim length from library: {cut_len}bp")
    elif args.trim_length is None and library.target_length is None:
        print(
            "Warning: library contains sequences of varying lengths; "
            "no trim length set. Reads will not be trimmed to a fixed length."
        )

    print(f"Processing reads from {args.reads}...")
    print(f"QC Settings: dropping reads with >={args.max_bad_freq}% bases under Q{args.min_base_q}")
    if args.output_trimmed:
        print(f"Saving trimmed reads to {args.output_trimmed}...")

    result = counting.process_reads(
        args.reads,
        library,
        trim_start=args.trim_start,
        trim_length=cut_len,
        revcomp=args.revcomp,
        trimmed_output_path=args.output_trimmed,
        min_base_q=args.min_base_q,
        max_bad_freq=args.max_bad_freq,
    )

    print(f"\nFinished processing {result.total_reads} reads in {result.elapsed_seconds:.2f} seconds.")
    if result.total_reads > 0:
        print(f"Skipped (Low QC): {result.skipped_qc} ({(result.skipped_qc / result.total_reads) * 100:.2f}%)")
    print(f"Reads Processed:  {result.passed_reads}")
    if result.passed_reads > 0:
        pct = (result.matched_reads / result.passed_reads) * 100
        print(f"Total Matched:    {result.matched_reads} ({pct:.2f}% of passed reads)")

    print(f"Writing results to {args.output}...")
    counting.write_output(args.output, result.counts)
    print("Done.")
    return 0


def _add_batch_parser(subparsers) -> None:
    p = subparsers.add_parser(
        "batch",
        help="Count reads for every FASTQ in a directory against one library.",
        description="Batch process Illumina reads (one count run per *R1_001.fastq.gz sample).",
    )
    p.add_argument("-i", "--input-dir", required=True, help="Directory containing FASTQ files")
    p.add_argument("-o", "--output-dir", required=True, help="Directory to save output CSVs")
    p.add_argument("-l", "--library", required=True, help="Path to reference library CSV")
    p.add_argument("--trim-start", type=int, default=0, help="Bases to trim from start (default: 0)")
    p.add_argument("--trim-length", type=int, default=None, help="Length to keep (default: match library length)")
    p.add_argument("--revcomp", action="store_true", help="Search for reverse complement")
    p.add_argument("--min-base-q", type=int, default=28, help="Min base quality for QC (default: 28)")
    p.add_argument("--max-bad-freq", type=float, default=10.0, help="Max bad base percentage for QC (default: 10%%)")
    p.set_defaults(func=_run_batch)


def _run_batch(args) -> int:
    summary = batch_mod.run_batch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        library_path=args.library,
        trim_start=args.trim_start,
        trim_length=args.trim_length,
        revcomp=args.revcomp,
        min_base_q=args.min_base_q,
        max_bad_freq=args.max_bad_freq,
    )
    # Non-zero exit if every sample failed (and there was at least one).
    if summary.processed == 0 and summary.failed > 0:
        return 1
    return 0


def _add_prep_parser(subparsers) -> None:
    p = subparsers.add_parser(
        "prep-library",
        help="Trim full sequences to read-length targets.",
        description="Trim full library sequences to read-length targets using a constant seed region.",
    )
    p.add_argument("--input", required=True, help="Path to full library CSV (needs ID and Sequence columns)")
    p.add_argument("--output", required=True, help="Path to write the trimmed library CSV (ID, Sequence)")
    p.add_argument("--seed", required=True, help="Constant seed region to split on (case-insensitive)")
    p.add_argument("--read-length", type=int, default=150, help="Sequencer read length (default: 150)")
    p.set_defaults(func=_run_prep)


def _run_prep(args) -> int:
    n = prep_mod.trim_targets(
        input_csv=args.input,
        output_csv=args.output,
        seed=args.seed,
        read_length=args.read_length,
    )
    print(f"Wrote {n} trimmed targets to {args.output}.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tcrlibs",
        description="Tools for counting TCR sequencing reads against reference libraries.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    _add_count_parser(subparsers)
    _add_batch_parser(subparsers)
    _add_prep_parser(subparsers)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    func = getattr(args, "func", None)
    if func is None:
        parser.print_help()
        return 0

    return func(args)


if __name__ == "__main__":
    sys.exit(main())
