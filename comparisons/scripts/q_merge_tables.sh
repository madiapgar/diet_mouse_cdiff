#!/bin/bash


echo "combining biom tables and rep seqs at baseline for all studies"
qiime feature-table merge \
    --i-tables ../data/qiime/newExp_comp_d15_table.qza \
    --i-tables ../data/qiime/oldExp_table.qza \
    --o-merged-table ../data/qiime/unfilt_allExp_comp_d15_table.qza

qiime feature-table merge-seqs \
    --i-data ../data/qiime/newExp_comp_d15_rep_seqs.qza \
    --i-data ../data/qiime/oldExp_rep_seqs.qza \
    --o-merged-data ../data/qiime/allExp_comp_d15_seqs.qza


echo "filtering combined biom table by the metadata"
qiime feature-table filter-samples \
    --i-table ../data/qiime/unfilt_allExp_comp_d15_table.qza \
    --m-metadata-file ../data/misc/oldNew_comp_d15_metadata.tsv \
    --o-filtered-table ../data/qiime/allExp_comp_d15_table.qza