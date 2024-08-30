#!/bin/bash

echo "making an EMP Paired End Sequences object out of raw fastq files for QIIME2 analysis!"
qiime tools import \
    --type EMPPairedEndSequences \
    --input-path comparisons/data/first_set_qiime/SEQ016/raw_seqs \
    --output-path comparisons/data/first_set_qiime/SEQ016/oldExp_s3_paired_end_seqs.qza