## step 1
## demux and dada2 of raw 16S sequences
import os

rule demux:
    input:
        in1 = expand(os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run_barcodes}_barcodes.txt"),
                                  run_barcodes=BARCODES),
        in2 = expand(os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_paired_end_seqs.qza"),
                        run=RAW_SEQS)
    output:
        out1 = os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_demux.qza"),
        out2 = os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_demux_details.qza")
    conda:
        QIIME
    shell:
        """
        qiime demux emp-paired \
            --m-barcodes-file {input.in1} \
            --m-barcodes-column BarcodeSequence \
            --i-seqs {input.in2} \
            --o-per-sample-sequences {output.out1} \
            --o-error-correction-details {output.out2} \
            --p-no-golay-error-correction
        """


rule demux_vis:
    input:
        os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_demux.qza")
    output:
        os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_demux.qzv")
    conda:
        QIIME
    shell:
        """
        qiime demux summarize \
            --i-data {input} \
            --o-visualization {output}
        """


rule dada2:
    input:
        os.path.join(DATASET_DIR, RAW_SEQ_DIR, "{run}_demux.qza")
    output:
        table = os.path.join(DATASET_DIR, "data/qiime/{run}_table.qza"),
        seqs = os.path.join(DATASET_DIR, "data/qiime/{run}_rep_seqs.qza"),
        stats = os.path.join(DATASET_DIR, "data/qiime/{run}_denoise_stats.qza")
    conda:
        QIIME
    params:
        trim_left_for=config["dada2_trim_left_for"],
        trim_left_rev=config["dada2_trim_left_rev"],
        trunc_len_for=config["dada2_trunc_len_for"],
        trunc_len_rev=config["dada2_trunc_len_rev"]
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --p-trim-left-f {params.trim_left_for} \
            --p-trim-left-r {params.trim_left_rev} \
            --p-trunc-len-f {params.trunc_len_for} \
            --p-trunc-len-r {params.trunc_len_rev} \
            --o-table {output.table} \
            --o-representative-sequences {output.seqs} \
            --o-denoising-stats {output.stats}
        """


rule merge_run_tables:
    input:
        table_list = expand(os.path.join(DATASET_DIR, "data/qiime/{run}_table.qza"),
                            run=RAW_SEQS),
        seqs_list = expand(os.path.join(DATASET_DIR, "data/qiime/{run}_rep_seqs.qza"),
                           run=RAW_SEQS)
    output:
        merged_table = os.path.join(DATASET_DIR, "data/qiime/merged_table.qza"),
        merged_seqs = os.path.join(DATASET_DIR, "data/qiime/merged_rep_seqs.qza")
    conda:
        QIIME
    shell:
        """
        qiime feature-table merge \
            --i-tables {input.table_list} \
            --o-merged-table {output.merged_table}

        qiime feature-table merge-seqs \
            --i-data {input.seqs_list} \
            --o-merged-data {output.merged_seqs}
        """