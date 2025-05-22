#!/bin/bash


echo "combining biom tables and rep seqs at baseline for all studies"
qiime feature-table merge \
    --i-tables ../../new_experiments/data/qiime/newExp_d15-d3_allBatch_table.qza \
    --i-tables ../../az_exp_combined/data/misc/comp_table.qza \
    --i-tables ../data/qiime/oldExp_table.qza \
    --o-merged-table ../data/qiime/unfilt_allExp_combined_table.qza

qiime feature-table merge-seqs \
    --i-data ../../new_experiments/data/qiime/newExp_d15-d3_allBatch_seqs.qza \
    --i-data ../../az_exp_combined/data/misc/comp_rep_seqs.qza \
    --i-data ../data/qiime/oldExp_rep_seqs.qza \
    --o-merged-data ../data/qiime/allExp_combined_d15-d3_seqs.qza


echo "filtering combined biom table by the metadata"
qiime feature-table filter-samples \
    --i-table ../data/qiime/unfilt_allExp_combined_table.qza \
    --m-metadata-file ../data/misc/oldNew_comp_d15-d3_metadata.tsv \
    --o-filtered-table ../data/qiime/allExp_combined_d15-d3_table.qza