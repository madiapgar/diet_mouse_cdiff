rule all:
    input:
        "data/bloodCulture_qiime/blood_culture_demux.qza",
        "data/bloodCulture_qiime/blood_culture_demux_details.qza",
        "data/bloodCulture_qiime/blood_culture_demux.qzv",
        "data/bloodCulture_qiime/blood_culture_table.qza",
        "data/bloodCulture_qiime/blood_culture_rep_seqs.qza",
        "data/bloodCulture_qiime/blood_culture_denoise_stats.qza",
        "databases/sepp-refs-silva-128.qza",
        "databases/silva-138-99-515-806-nb-classifier.qza",
        "data/bloodCulture_qiime/tree.qza",
        "data/bloodCulture_qiime/placements.qza",
        "data/bloodCulture_qiime/filt_bc_table.qza",
        "data/bloodCulture_qiime/rem_bc_table.qza",
        "data/bloodCulture_qiime/taxonomy.qza",
        "data/bloodCulture_qiime/bc_tax_filt.qza",
        "data/bloodCulture_qiime/bc_tax_filt.qzv",
        "data/bloodCulture_qiime/filt_otu_table.qza",
        "data/bloodCulture_qiime/tax_barplot.qzv"


rule demultiplex:
    input:
        barcodes = "data/bloodCulture_qiime/bloodCulture_barcodes.txt",
        seqs = "data/bloodCulture_qiime/bloodCulture_paired_end_seqs.qza"
    output:
        demux = "data/bloodCulture_qiime/blood_culture_demux.qza",
        demux_details = "data/bloodCulture_qiime/blood_culture_demux_details.qza",
        demux_vis = "data/bloodCulture_qiime/blood_culture_demux.qzv"
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
        demux = "data/bloodCulture_qiime/blood_culture_demux.qza"
    output:
        table = "data/bloodCulture_qiime/blood_culture_table.qza",
        rep_seqs = "data/bloodCulture_qiime/blood_culture_rep_seqs.qza",
        denoise_stats = "data/bloodCulture_qiime/blood_culture_denoise_stats.qza"
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
        seqs = "data/bloodCulture_qiime/blood_culture_rep_seqs.qza",
        silva_ref = "databases/sepp-refs-silva-128.qza"
    output:
        tree = "data/bloodCulture_qiime/tree.qza",
        placements = "data/bloodCulture_qiime/placements.qza"
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
        table = "data/bloodCulture_qiime/blood_culture_table.qza",
        tree = "data/bloodCulture_qiime/tree.qza"
    output:
        filt_table = "data/bloodCulture_qiime/filt_bc_table.qza",
        rem_table = "data/bloodCulture_qiime/rem_bc_table.qza"
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
        seqs = "data/bloodCulture_qiime/blood_culture_rep_seqs.qza"
    output:
        taxonomy = "data/bloodCulture_qiime/taxonomy.qza"
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
        filt_table = "data/bloodCulture_qiime/filt_bc_table.qza",
        taxonomy = "data/bloodCulture_qiime/taxonomy.qza"
    output:
        tax_filt = "data/bloodCulture_qiime/bc_tax_filt.qza",
        tax_filt_vis = "data/bloodCulture_qiime/bc_tax_filt.qzv"
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
        tax_filt = "data/bloodCulture_qiime/bc_tax_filt.qza",
        metadata = "data/misc/blood_culture_metadata.tsv"
    output:
        meta_filt_tax = "data/bloodCulture_qiime/filt_otu_table.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.meta_filt_tax}
        """


rule make_taxa_barplot:
    input:
        table = "data/bloodCulture_qiime/filt_otu_table.qza",
        taxonomy = "data/bloodCulture_qiime/taxonomy.qza",
        metadata = "data/misc/blood_culture_metadata.tsv"
    output:
        tax_barplot = "data/bloodCulture_qiime/tax_barplot.qzv"
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