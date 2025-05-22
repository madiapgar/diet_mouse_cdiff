#!/bin/bash

##echo "making an EMP Paired End Sequences object out of raw fastq files for QIIME2 analysis!"
##qiime tools import \
    ##--type EMPPairedEndSequences \ 
    ##--input-path raw_seqs \
    ##--output-path newExp_d15-d3_paired_end.qza


##echo "running DADA2 on demultiplexed sequences!"
##qiime dada2 denoise-paired \
##    --i-demultiplexed-seqs ../data/SEQ088/newExp_d15-d3_demux.qza \
##    --p-trim-left-f 13 \
##    --p-trim-left-r 13 \
##    --p-trunc-len-f 230 \
##    --p-trunc-len-r 160 \
##    --o-table ../data/SEQ088/newExp_d15-d3_table.qza \
##    --o-representative-sequences ../data/SEQ088/newExp_d15-d3_rep_seqs.qza \
##    --o-denoising-stats ../data/SEQ088/newExp_d15-d3_denoise_stats.qza
##
##
##echo "filtering BIOM table and rep-seqs by the metadata file"
##qiime feature-table filter-samples \
##    --i-table ../data/SEQ088/newExp_d15-d3_table.qza \
##    --m-metadata-file ../data/SEQ088/newExp_d15-d3_metadata.txt \
##    --o-filtered-table ../data/SEQ088/newExp_d15-d3_filt_table.qza

echo "combining biom tables and rep seqs for both batches of the new experiments"
qiime feature-table merge \
    --i-tables ../data/CDD02_qiime/newExp_d15-d3_filt_table.qza \
    --i-tables ../data/CDD03_qiime/new_stool_table.qza \
    --o-merged-table ../data/qiime/newExp_d15-d3_allBatch_table.qza

qiime feature-table merge-seqs \
    --i-data ../data/CDD02_qiime/newExp_d15-d3_rep_seqs.qza \
    --i-data ../data/CDD03_qiime/new_stool_rep_seqs.qza \
    --o-merged-data ../data/qiime/newExp_d15-d3_allBatch_seqs.qza