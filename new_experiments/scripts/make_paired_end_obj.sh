#!/bin/bash

echo "making an EMP Paired End Sequences object out of raw fastq files for QIIME2 analysis!"
qiime tools import \
    --type EMPPairedEndSequences \
    --input-path ./data/SEQ096/raw_seqs \
    --output-path ./data/SEQ096/newExp_cult-stool_paired_end_seqs.qza