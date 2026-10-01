"""Core logic for counting sequencing reads against a reference library.

Ported from the standalone ``count_reads.py`` script. The functions here
return data structures rather than relying on ``sys.exit``/print-only side
effects, so they can be reused both by the ``count`` CLI subcommand and by the
``batch`` runner.
"""

import csv
import gzip
import time
from dataclasses import dataclass, field
from typing import Optional


def get_reverse_complement(seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


@dataclass
class Library:
    """A loaded reference library.

    Attributes:
        sequence_to_id: maps a (upper-cased) sequence to its library ID.
        target_length: the single sequence length if the library is uniform,
            otherwise ``None`` (exact matching then needs an explicit trim).
    """

    sequence_to_id: dict = field(default_factory=dict)
    target_length: Optional[int] = None

    def new_counts(self) -> dict:
        """Return a fresh {library_id: 0} counts dict for a run."""
        return {seq_id: 0 for seq_id in self.sequence_to_id.values()}


@dataclass
class CountResult:
    """Outcome of processing a FASTQ file."""

    counts: dict
    total_reads: int = 0
    skipped_qc: int = 0
    matched_reads: int = 0
    elapsed_seconds: float = 0.0

    @property
    def passed_reads(self) -> int:
        return self.total_reads - self.skipped_qc


def load_library(filepath: str) -> Library:
    """Load the reference library CSV (format: ID, Sequence).

    Mirrors the original heuristic: the first row is treated as a header when
    its second column looks like a header (short, or contains "seq").
    """
    sequence_to_id: dict = {}
    sequence_lengths = set()

    with open(filepath, mode="r") as f:
        reader = csv.reader(f)
        header = next(reader, None)

        # Heuristic to check if first row is a header.
        if header and (len(header[1]) < 10 or "seq" in header[1].lower()):
            # It was likely a header; leave the reader positioned after it.
            pass
        else:
            # It was not a header; rewind so it is processed as data.
            f.seek(0)

        for row in reader:
            if len(row) < 2:
                continue
            seq_id = row[0].strip()
            sequence = row[1].strip().upper()

            sequence_to_id[sequence] = seq_id
            sequence_lengths.add(len(sequence))

    target_length: Optional[int] = None
    if len(sequence_lengths) == 1:
        target_length = next(iter(sequence_lengths))

    return Library(sequence_to_id=sequence_to_id, target_length=target_length)


def process_reads(
    reads_path: str,
    library: Library,
    trim_start: int = 0,
    trim_length: Optional[int] = None,
    revcomp: bool = False,
    trimmed_output_path: Optional[str] = None,
    min_base_q: int = 28,
    max_bad_freq: float = 10.0,
) -> CountResult:
    """Parse a FASTQ(.gz) file and count exact matches with QC filtering.

    Reads are dropped when the fraction of bases below Q(``min_base_q``) is
    greater than or equal to ``max_bad_freq`` percent. Quality is assumed to be
    Phred+33 encoded. Trimming is applied before the optional reverse
    complement, matching the original behavior.
    """
    counts = library.new_counts()
    library_map = library.sequence_to_id

    if reads_path.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"  # read text mode
    else:
        open_func = open
        mode = "r"

    total_reads = 0
    skipped_qc = 0
    matched_reads = 0

    # Phred+33: a base is "bad" when ord(char) < min_base_q + 33.
    qc_char_limit = min_base_q + 33
    max_bad_ratio = max_bad_freq / 100.0

    trimmed_file = None
    if trimmed_output_path:
        trimmed_file = open(trimmed_output_path, "w")
        trimmed_file.write("Read_ID,Trimmed_Sequence,Matched_Library_ID\n")

    start_time = time.time()

    try:
        with open_func(reads_path, mode) as f:
            while True:
                # FASTQ format: 4 lines per record.
                identifier = f.readline()
                if not identifier:
                    break  # End of file.

                sequence = f.readline().strip().upper()
                f.readline()  # "+" separator line.
                quality = f.readline().strip()

                total_reads += 1

                # --- QC step ---
                low_quality_count = 0
                for char in quality:
                    if ord(char) < qc_char_limit:
                        low_quality_count += 1

                if len(quality) > 0 and (low_quality_count / len(quality)) >= max_bad_ratio:
                    skipped_qc += 1
                    continue  # Skip this read entirely.
                # ----------------

                # Trim sequence.
                if trim_length:
                    query_seq = sequence[trim_start : trim_start + trim_length]
                else:
                    query_seq = sequence[trim_start:]

                # Apply reverse complement if requested (after trimming).
                if revcomp:
                    query_seq = get_reverse_complement(query_seq)

                matched_id = ""
                if query_seq in library_map:
                    seq_id = library_map[query_seq]
                    counts[seq_id] += 1
                    matched_reads += 1
                    matched_id = seq_id

                if trimmed_file:
                    # Strip newline and the leading '@' from the FASTQ header.
                    clean_id = identifier.strip()[1:]
                    trimmed_file.write(f"{clean_id},{query_seq},{matched_id}\n")
    finally:
        if trimmed_file:
            trimmed_file.close()

    elapsed = time.time() - start_time

    return CountResult(
        counts=counts,
        total_reads=total_reads,
        skipped_qc=skipped_qc,
        matched_reads=matched_reads,
        elapsed_seconds=elapsed,
    )


def write_output(output_path: str, counts: dict) -> None:
    """Write counts to a CSV sorted by count descending."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Library_ID", "Count"])

        sorted_results = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        for seq_id, count in sorted_results:
            writer.writerow([seq_id, count])


def resolve_trim_length(
    explicit_trim_length: Optional[int], library: Library
) -> Optional[int]:
    """Pick the trim length: explicit value, else the library's uniform length."""
    if explicit_trim_length is not None:
        return explicit_trim_length
    return library.target_length
