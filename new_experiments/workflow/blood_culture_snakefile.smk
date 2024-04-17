configfile: "workflow/config.yml"


rule all:
    input:
        config["path_to_seq_dir"] + config["sampleType_prefix"] + "_paired_end_seqs.qza",
        "data/qiime/" + config["sampleType_prefix"] + "_demux.qza",
        "data/qiime/" + config["sampleType_prefix"] + "_demux_details.qza",
        "data/qiime/" + config["sampleType_prefix"] + "_demux.qzv",
        "data/qiime/" + config["sampleType_prefix"] + "_table.qza",
        "data/qiime/" + config["sampleType_prefix"] + "_rep_seqs.qza",
        "data/qiime/" + config["sampleType_prefix"] + "_denoise_stats.qza",
        "data/qiime/tree.qza",
        "data/qiime/placements.qza",
        "data/qiime/filt_" + config["sampleType_prefix"] + "_table.qza",
        "data/qiime/rem_" + config["sampleType_prefix"] + "_table.qza",
        "data/qiime/taxonomy.qza",
        "data/qiime/tax_filt.qza",
        "data/qiime/tax_filt.qzv",
        "data/qiime/filt_otu_table.qza",
        "data/qiime/tax_barplot.qzv"


rule import_paired_end:
    input:
        seqs_directory = config["raw_seq_dir"]
    output:
        seqs = config["path_to_seq_dir"] + config["sampleType_prefix"] + "_paired_end_seqs.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools import \
            --type EMPPairedEndSequences \
            --input-path {input.seqs_directory} \
            --output-path {output.seqs}
        """


rule demultiplex:
    input:
        barcodes = config["path_to_seq_dir"] + config["sampleType_prefix"] + "_barcodes.txt",
        seqs = config["path_to_seq_dir"] + config["sampleType_prefix"] + "_paired_end_seqs.qza"
    output:
        demux = "data/qiime/" + config["sampleType_prefix"] + "_demux.qza",
        demux_details = "data/qiime/" + config["sampleType_prefix"] + "_demux_details.qza",
        demux_vis = "data/qiime/" + config["sampleType_prefix"] + "_demux.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime demux emp-paired \
            --m-barcodes-file {input.barcodes} \
            --m-barcodes-column BarcodeSequence \
            --i-seqs {input.seqs} \
            --o-per-sample-sequences {output.demux} \
            --o-error-correction-details {output.demux_details} \
            --p-no-golay-error-correction
        
        qiime demux summarize \
            --i-data {output.demux} \
            --o-visualization {output.demux_vis}
        """


rule dada2:
    input:
        demux = "data/qiime/" + config["sampleType_prefix"] + "_demux.qza"
    output:
        table = "data/qiime/" + config["sampleType_prefix"] + "_table.qza",
        rep_seqs = "data/qiime/" + config["sampleType_prefix"] + "_rep_seqs.qza",
        denoise_stats = "data/qiime/" + config["sampleType_prefix"] + "_denoise_stats.qza"
    conda:
        "qiime2-2023.5"
    params:
        trim_left_for=config["dada2_trim_left_for"],
        trim_left_rev=config["dada2_trim_left_rev"],
        trunc_len_for=config["dada2_trunc_len_for"],
        trunc_len_rev=config["dada2_trunc_len_rev"]
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input.demux} \
            --p-trim-left-f {params.trim_left_for} \
            --p-trim-left-r {params.trim_left_rev} \
            --p-trunc-len-f {params.trunc_len_for} \
            --p-trunc-len-r {params.trunc_len_rev} \
            --o-table {output.table} \
            --o-representative-sequences {output.rep_seqs} \
            --o-denoising-stats {output.denoise_stats}
        """


rule get_reference_databases:
    output:
        "databases/sepp-refs-silva-128.qza",
        "databases/silva-138-99-515-806-nb-classifier.qza"
    shell:
        """
        wget https://data.qiime2.org/2023.5/common/sepp-refs-silva-128.qza -P ./databases/
        wget https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza -P ./databases/
        """


rule sepp_phylo_tree:
    input:
        seqs = "data/qiime/" + config["sampleType_prefix"] + "_rep_seqs.qza",
        silva_ref = "databases/sepp-refs-silva-128.qza"
    output:
        tree = "data/qiime/tree.qza",
        placements = "data/qiime/placements.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.seqs} \
            --i-reference-database {input.silva_ref} \
            --o-tree {output.tree} \
            --o-placements {output.placements}
        """


rule sepp_ASV_filtering:
    input:
        table = "data/qiime/" + config["sampleType_prefix"] + "_table.qza",
        tree = "data/qiime/tree.qza"
    output:
        filt_table = "data/qiime/filt_" + config["sampleType_prefix"] + "_table.qza",
        rem_table = "data/qiime/rem_" + config["sampleType_prefix"] + "_table.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime fragment-insertion filter-features \
            --i-table {input.table} \
            --i-tree {input.tree} \
            --o-filtered-table {output.filt_table} \
            --o-removed-table {output.rem_table}
        """


rule taxonomic_classification:
    input:
        silva_class = "databases/silva-138-99-515-806-nb-classifier.qza",
        seqs = "data/qiime/" + config["sampleType_prefix"] + "_rep_seqs.qza"
    output:
        taxonomy = "data/qiime/taxonomy.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.silva_class} \
            --i-reads {input.seqs} \
            --o-classification {output.taxonomy}
        """


rule filter_taxonomy:
    input:
        filt_table = "data/qiime/filt_" + config["sampleType_prefix"] + "_table.qza",
        taxonomy = "data/qiime/taxonomy.qza"
    output:
        tax_filt = "data/qiime/tax_filt.qza",
        tax_filt_vis = "data/qiime/tax_filt.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.filt_table} \
            --i-taxonomy {input.taxonomy} \
            --p-include p_ \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table {output.tax_filt}
        
        qiime feature-table summarize \
            --i-table {output.tax_filt} \
            --o-visualization {output.tax_filt_vis}
        """


rule filter_biom_w_metadata:
    input:
        tax_filt = "data/qiime/tax_filt.qza",
        metadata = config["metadata"]
    output:
        otu_table = "data/qiime/filt_otu_table.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.otu_table}
        """


rule make_taxa_barplot:
    input:
        table = "data/qiime/filt_otu_table.qza",
        taxonomy = "data/qiime/taxonomy.qza",
        metadata = config["metadata"]
    output:
        tax_barplot = "data/qiime/tax_barplot.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.tax_barplot}
        """