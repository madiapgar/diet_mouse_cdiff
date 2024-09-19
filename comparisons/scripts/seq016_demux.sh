#!/bin/bash

## I had to demultiplex SEQ016 differently because the barcodes were golay-corrected due to the lab using different primers
## at that point in time 
## that's also why I had to add (only) the p-rev-comp-mapping-barcodes flag to the demultiplex command

echo "demultiplexing raw paired-end sequences"

qiime demux emp-paired \
    --m-barcodes-file ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_barcodes.txt \
    --m-barcodes-column BarcodeSequence \
    --i-seqs ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_paired_end_seqs.qza \
    --o-per-sample-sequences ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_demux.qza \
    --o-error-correction-details ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_demux_details.qza \
    --p-golay-error-correction \
    --p-rev-comp-mapping-barcodes


echo "creating demux visualization"

qiime demux summarize \
    --i-data ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_demux.qza \
    --o-visualization ./data/first_set_qiime/SEQ016/s16_raw_seqs/oldExp_s3-016_demux.qzv
