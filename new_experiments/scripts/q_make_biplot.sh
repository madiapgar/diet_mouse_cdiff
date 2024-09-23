#!/bin/bash

echo "making a biplot using pcoa results (since that's what we have for everything else)"
## make relative frequency table from otu table (should use rarefied table but I can't use rarefied results)
qiime feature-table relative-frequency \
    --i-table ../data/qiime/core_outputs/rarefied_table.qza \
    --o-relative-frequency-table ../data/qiime/core_outputs/relative_rarefied_table.qza

## make biplot for unweighted unifrac pcoa
qiime diversity pcoa-biplot \
    --i-pcoa ../data/qiime/core_outputs/unweighted_unifrac_pcoa_results.qza \
    --i-features ../data/qiime/core_outputs/relative_rarefied_table.qza \
    --o-biplot ../data/qiime/core_outputs/uu_biplot_matrix.qza

## make biplot for weighted unifrac pcoa
qiime diversity pcoa-biplot \
    --i-pcoa ../data/qiime/core_outputs/weighted_unifrac_pcoa_results.qza \
    --i-features ../data/qiime/core_outputs/relative_rarefied_table.qza \
    --o-biplot ../data/qiime/core_outputs/wu_biplot_matrix.qza
