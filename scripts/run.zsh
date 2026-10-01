#!/bin/zsh
cd /Volumes/common/TCR_libraries/screens/TCRCodex_Hiroyasu_expt1

python ./scripts/batch_runner.py \
-i /Volumes/common/TIRTL/iSeq_rawData/20251118_FS10003214_33_BWB90303-2330/fastq \
-o ./outs \
-l ./scripts/library_table_clean.csv \
-s ./scripts/count_reads.py \
--trim-start 42 \
--revcomp





