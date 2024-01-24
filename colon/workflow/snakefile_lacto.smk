rule all:
    input:
        "data/lacto_qiime/tree.qza",
        "data/lacto_qiime/placements.qza",
        "data/lacto_qiime/filt_lacto_table.qza",
        "data/lacto_qiime/rem_lacto_table.qza",
        "data/lacto_qiime/taxonomy.qza",
        "data/lacto_qiime/taxonomy_filtered.qza",
        "data/lacto_qiime/taxonomy_filtered.qzv",
        "data/lacto_qiime/tax_filt_actual.qza"

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
        lacto_rep_seqs = "data/misc/lactoOnly_mergedDietAim1rep-seqs_051523-Copy1.qza",
        silva_ref = "databases/sepp-refs-silva-128.qza"
    output:
        tree = "data/lacto_qiime/tree.qza",
        placements = "data/lacto_qiime/placements.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.lacto_rep_seqs} \
            --i-reference-database {input.silva_ref} \
            --o-tree {output.tree} \
            --o-placements {output.placements}
        """


rule sepp_ASV_filtering:
    input:
        lacto_table = "data/misc/lacto_only_slva_mergedDietAim1_051523-Copy1.qza",
        tree = "data/lacto_qiime/tree.qza"
    output:
        filt_lacto_table = "data/lacto_qiime/filt_lacto_table.qza",
        rem_lacto_table = "data/lacto_qiime/rem_lacto_table.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime fragment-insertion filter-features \
            --i-table {input.lacto_table} \
            --i-tree {input.tree} \
            --o-filtered-table {output.filt_lacto_table} \
            --o-removed-table {output.rem_lacto_table}
        """


rule taxonomic_classification:
    input:
        silva_class = "databases/silva-138-99-515-806-nb-classifier.qza",
        lacto_rep_seqs = "data/misc/lactoOnly_mergedDietAim1rep-seqs_051523-Copy1.qza"
    output:
        taxonomy = "data/lacto_qiime/taxonomy.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.silva_class} \
            --i-reads {input.lacto_rep_seqs} \
            --o-classification {output.taxonomy}
        """

## rerun after total sum scaling
## hopefully the second qiime command works 
rule filter_taxonomy:
    input:
        filt_lacto_table = "data/lacto_qiime/filt_lacto_table.qza",
        taxonomy = "data/lacto_qiime/taxonomy.qza"
    output:
        tax_filt_table = "data/lacto_qiime/taxonomy_filtered.qza",
        tax_filt_vis = "data/lacto_qiime/taxonomy_filtered.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.filt_lacto_table} \
            --i-taxonomy {input.taxonomy} \
            --p-include p_ \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table {output.tax_filt_table}
        
        qiime feature-table summarize \
            --i-table {output.tax_filt_table} \
            --o-visualization {output.tax_filt_vis}
        """

rule pre_core_metrics_filter:
    input:
        tax_filt_table = "data/lacto_qiime/taxonomy_filtered.qza",
        metadata = "data/misc/merged_metadata1.tsv"
    output:
        meta_filt_tax = "data/lacto_qiime/tax_filt_actual.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt_table} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.meta_filt_tax}
        """