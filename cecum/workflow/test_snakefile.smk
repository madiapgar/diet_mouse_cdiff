run_dir=["SEQ069/s69", "SEQ070/s70", "SEQ071/s71"]
cecal_dir=["SEQ069/cecal_s69", "SEQ070/cecal_s70", "SEQ071/cecal_s71"]

rule all:
    input:
        expand("data/cecal_qiime_upper/{run}_demux.qza",
                run=run_dir),
        expand("data/cecal_qiime_upper/{run}_demux_details.qza",
                run=run_dir),
        expand("data/cecal_qiime_upper/{run}_demux.qzv",
                run=run_dir),
        expand("data/cecal_qiime/{run}_table.qza",
                run=run_dir),
        expand("data/cecal_qiime/{run}_rep_seqs.qza",
                run=run_dir),
        expand("data/cecal_qiime/{run}_denoise_stats.qza",
                run=run_dir)


rule demux:
    input:
        in1 = expand("data/cecal_qiime_upper/{run_cecal}_barcodes.txt",
                        run_cecal=cecal_dir),
        in2 = expand("data/cecal_qiime_upper/{run}_paired_end_seqs.qza",
                        run=run_dir)
    output:
        out1 = "data/cecal_qiime_upper/{run}_demux.qza",
        out2 = "data/cecal_qiime_upper/{run}_demux_details.qza"
    conda:
        "qiime2-2023.5"
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
        "data/cecal_qiime_upper/{run}_demux.qza"
    output:
        "data/cecal_qiime_upper/{run}_demux.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime demux summarize \
            --i-data {input} \
            --o-visualization {output}
        """


rule dada2:
    input:
       "data/cecal_qiime_upper/{run}_demux.qza"
    output:
        table = "data/cecal_qiime/{run}_table.qza",
        seqs = "data/cecal_qiime/{run}_rep_seqs.qza",
        stats = "data/cecal_qiime/{run}_denoise_stats.qza"
    conda:
        "qiime2-2023.5"
    params:
        trim_left_for=13,
        trim_left_rev=13,
        trunc_len_for=230,
        trunc_len_rev=160
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


        