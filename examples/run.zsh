#!/bin/zsh
# Example batch invocation using the tcrlibs CLI.
#
# This replaces the original scripts/run.zsh (which invoked batch_runner.py
# and count_reads.py directly). Install the package first with:
#
#   pip install -e .
#
# Then adjust the paths below for your experiment and run this script.

cd /Volumes/common/TCR_libraries/screens/TCRCodex_Hiroyasu_expt1

tcrlibs batch \
  -i /Volumes/common/TIRTL/iSeq_rawData/20251118_FS10003214_33_BWB90303-2330/fastq \
  -o ./outs \
  -l ./scripts/library_table_clean.csv \
  --trim-start 42 \
  --revcomp
