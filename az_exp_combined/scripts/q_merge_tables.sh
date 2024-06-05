#!/bin/bash

echo "merging BIOM tables"
qiime feature-table merge \
    --i-tables ../data/misc/euk_filt_mergedDietAim1table_051523-Copy1.qza \
    --i-tables ../../new_experiments/data/qiime/newExp_d15-d3_filt_table.qza \
    --o-merged-table ../data/qiime/newExp_comp_d15_table.qza


echo "merging rep seqs"
qiime feature-table merge-seqs \
    --i-data ../data/misc/euk_filt-mergedDietAim1rep-seqs_051523-Copy1.qza \
    --i-data ../../new_experiments/data/qiime/newExp_d15-d3_rep_seqs.qza \
    --o-merged-data ../data/qiime/newExp_comp_d15_rep_seqs.qza