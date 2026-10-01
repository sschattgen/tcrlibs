"""Batch processing: count every FASTQ in a directory against one library.

Ported from the standalone ``batch_runner.py`` script. Instead of shelling out
to ``count_reads.py`` once per sample, this calls the counting functions
directly, loading the reference library a single time and reusing it across all
samples.
"""

import glob
import os
from dataclasses import dataclass
from typing import Optional

from . import counting

DEFAULT_PATTERN = "*R1_001.fastq.gz"

@dataclass
class BatchSummary:
    """Outcome of a batch run."""

    processed: int = 0
    failed: int = 0

def _derive_sample_name(filename: str, pattern: str) -> str:
    """Derive a sample name from a filename given the glob pattern.

    The literal text that follows the first ``*`` in the pattern (e.g.
    ``R1_001.fastq.gz`` for the default) marks the end of the sample name.
    Everything up to that marker, with a trailing ``_`` stripped, is the sample
    name. If the marker is not found, known FASTQ extensions are dropped as a
    fallback.
    """
    star = pattern.find("*")
    suffix = pattern[star + 1 :] if star != -1 else ""

    if suffix and suffix in filename:
        name = filename.split(suffix)[0]
    else:
        # No usable suffix marker; fall back to stripping known extensions.
        name = filename
        for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"):
            if name.endswith(ext):
                name = name[: -len(ext)]
                break

    return name.rstrip("_")

def discover_samples(input_dir: str, pattern: str = DEFAULT_PATTERN):
    """Return sorted (filepath, sample_name) pairs for files matching ``pattern``.

    ``pattern`` is a glob applied within ``input_dir`` (e.g. the default
    ``*R1_001.fastq.gz``). The sample name is derived from the portion of each
    filename preceding the pattern's literal suffix.
    """
    search_pattern = os.path.join(input_dir, pattern)
    files = sorted(glob.glob(search_pattern))

    samples = []
    for filepath in files:
        filename = os.path.basename(filepath)
        sample_name = _derive_sample_name(filename, pattern)
        samples.append((filepath, sample_name))
    return samples

def run_batch(
    input_dir: str,
    output_dir: str,
    library_path: str,
    pattern: str = DEFAULT_PATTERN,
    trim_start: int = 0,
    trim_length: Optional[int] = None,
    revcomp: bool = False,
    min_base_q: int = 28,
    max_bad_freq: float = 10.0,
) -> BatchSummary:
    """Process every file matching ``pattern`` in ``input_dir``.

    Writes ``{sample}_counts.csv`` and ``{sample}_trimmed.csv`` into
    ``output_dir`` for each sample. A failure on one sample is logged and the
    batch continues, matching the original behavior.
    """
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    samples = discover_samples(input_dir, pattern)
    if not samples:
        print(f"No files matching '{pattern}' found in {input_dir}")
        return BatchSummary()

    print(f"Found {len(samples)} samples to process.\n")

    # Load the library once and reuse it across all samples.
    print(f"Loading library from {library_path}...")
    library = counting.load_library(library_path)
    print(f"Loaded {len(library.sequence_to_id)} unique sequences.")

    cut_len = counting.resolve_trim_length(trim_length, library)
    if trim_length is None and library.target_length is not None:
        print(f"Auto-detected trim length from library: {cut_len}bp")

    summary = BatchSummary()

    for filepath, sample_name in samples:
        output_csv = os.path.join(output_dir, f"{sample_name}_counts.csv")
        trimmed_csv = os.path.join(output_dir, f"{sample_name}_trimmed.csv")

        print(f"--> Processing {sample_name}...")

        try:
            result = counting.process_reads(
                filepath,
                library,
                trim_start=trim_start,
                trim_length=cut_len,
                revcomp=revcomp,
                trimmed_output_path=trimmed_csv,
                min_base_q=min_base_q,
                max_bad_freq=max_bad_freq,
            )
            counting.write_output(output_csv, result.counts)

            passed = result.passed_reads
            matched_pct = (result.matched_reads / passed * 100) if passed > 0 else 0.0
            print(
                f"    {result.total_reads} reads, {result.skipped_qc} low-QC, "
                f"{result.matched_reads} matched ({matched_pct:.2f}% of passed)."
            )
            summary.processed += 1
        except Exception as e:  # noqa: BLE001 - keep batch going on any single-sample error
            print(f"Error processing {os.path.basename(filepath)}: {e}")
            summary.failed += 1
            # One bad file should not stop the whole batch.
            continue

    print("\nBatch processing complete.")
    return summary

