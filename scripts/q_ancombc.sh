#!/bin/bash

echo "----running ANCOMBC----"

qiime composition ancombc \
    --i-table ./data/qiime/tax_filt_actual.qza \
    --m-metadata-file ./data/misc/updated_metadata.tsv \
    --p-formula 'diet' \
    --p-p-adj-method 'BH' \
    --o-differentials ./data/qiime/test_ancombc.qza