# TCRlibs

Tools for counting TCR sequencing reads against reference libraries.

TCRlibs packages a small bioinformatics workflow as a single command-line tool:

- `prep-library` — trim full-length library sequences down to read-length targets
- `count` — count the reads in one FASTQ that exactly match a reference library
- `batch` — run `count` across every FASTQ in a directory, reusing one loaded library

## Installation

Install as an editable package (recommended for internal/lab use):

```bash
pip install -e .
```

For development (adds `pytest`):

```bash
pip install -e ".[dev]"
```

This installs a single `tcrlibs` command.

## Usage

```bash
tcrlibs --help
tcrlibs <command> --help
```

### prep-library

Trim full-length sequences to read-length targets. For each sequence, the tool
locates a constant seed region and keeps the `read_length - len(seed)` bases
immediately preceding it. The input CSV needs `ID` and `Sequence` columns
(extra columns such as `Subpool` are ignored); the output has `ID,Sequence`.

```bash
tcrlibs prep-library \
  --input library_table.csv \
  --output library_table_clean.csv \
  --seed GAGGACCTGAACAAGGTGTTTCCTCCAGAGGTGGCCGTGTTC \
  --read-length 150
```

### count

Count reads from a single FASTQ (or FASTQ.GZ) against a reference library.
Reads are quality-filtered (Phred+33), trimmed, optionally reverse-complemented,
then matched exactly against the library. If every library sequence is the same
length, the trim length is auto-detected.

```bash
tcrlibs count \
  -r sample_R1_001.fastq.gz \
  -l library_table_clean.csv \
  -o sample_counts.csv \
  --trim-start 42 \
  --revcomp
```

Options:

| Flag | Default | Description |
| --- | --- | --- |
| `-r`, `--reads` | required | Input FASTQ or FASTQ.GZ |
| `-l`, `--library` | required | Reference library CSV (`ID,Sequence`) |
| `-o`, `--output` | required | Output counts CSV |
| `--output-trimmed` | off | Also save every trimmed read and its match status |
| `--trim-start` | 0 | Bases to trim from the start of each read |
| `--trim-length` | library length | Length to keep after trimming |
| `--revcomp` | off | Match the reverse complement |
| `--min-base-q` | 28 | Minimum Phred score for a "good" base |
| `--max-bad-freq` | 10.0 | Max % of low-quality bases before dropping a read |

### batch

Process every `*R1_001.fastq.gz` file in a directory. The library is loaded once
and reused across all samples. Each sample writes `{sample}_counts.csv` and
`{sample}_trimmed.csv` into the output directory. A failure on one sample is
logged and the batch continues.

```bash
tcrlibs batch \
  -i /path/to/fastq \
  -o ./outs \
  -l library_table_clean.csv \
  --trim-start 42 \
  --revcomp
```

See `examples/run.zsh` for a complete batch invocation.

## Development

Run the test suite:

```bash
pytest
```

## Project layout

```
src/tcrlibs/
  cli.py        # argparse CLI: count, batch, prep-library
  counting.py   # library loading + read counting core
  batch.py      # directory discovery + per-sample loop
  prep.py       # seed-based target trimming
tests/          # pytest smoke tests + fixtures
scripts/        # original standalone scripts, kept for reference
```
