## step 2
## phylogenetic classification

## table and rep seqs inputs will either be the merged ones or ones provided 

rule get_reference_databases:
    output:
        DATASET_DIR + "databases/sepp-refs-silva-128.qza",
        DATASET_DIR + "databases/silva-138-99-515-806-nb-classifier.qza"
    shell:
        """
        wget https://data.qiime2.org/2023.5/common/sepp-refs-silva-128.qza -P ./databases/
        wget https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza -P ./databases/
        """


rule sepp_phylo_tree:
    input:
        seqs = DATASET_DIR + REP_SEQS,
        silva_ref = DATASET_DIR + "databases/sepp-refs-silva-128.qza"
    output:
        tree = DATASET_DIR + "data/qiime/tree.qza",
        placements = DATASET_DIR + "data/qiime/placements.qza"
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
        table = DATASET_DIR + BIOM,
        tree = DATASET_DIR + "data/qiime/tree.qza"
    output:
        filt_table = DATASET_DIR + "data/qiime/filt_table.qza",
        rem_table = DATASET_DIR + "data/qiime/rem_table.qza"
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
        silva_class = DATASET_DIR + "databases/silva-138-99-515-806-nb-classifier.qza",
        seqs = DATASET_DIR + REP_SEQS
    output:
        taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza"
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
        filt_table = DATASET_DIR + "data/qiime/filt_table.qza",
        taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza"
    output:
        tax_filt = DATASET_DIR + "data/qiime/taxonomy_filtered.qza",
        tax_filt_vis = DATASET_DIR + "data/qiime/taxonomy_filtered.qzv"
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

rule create_lacto_table:
    input:
        filt_table = DATASET_DIR + "data/qiime/filt_table.qza",
        taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza",
        seqs = DATASET_DIR + REP_SEQS
    output:
        lacto_table = DATASET_DIR + "data/qiime/lacto_cecal_table.qza",
        lacto_rep_seqs = DATASET_DIR + "data/qiime/lacto_rep_seqs.qza" 
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.filt_table} \
            --i-taxonomy {input.taxonomy} \
            --p-include Lactococcus \
            --o-filtered-table {output.lacto_table}
        
        qiime feature-table filter-seqs \
            --i-data {input.seqs} \
            --i-table {output.lacto_table} \
            --o-filtered-data {output.lacto_rep_seqs}
        """


rule convert_to_fasta:
    input:
        lacto_rep_seqs = DATASET_DIR + "data/qiime/lacto_rep_seqs.qza"
    output:
        output_path = DATASET_DIR + "data/qiime/dna-sequences.fasta",
        lacto_fasta = DATASET_DIR + "data/qiime/lactoOnly_rep_seqs.fasta"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path {input.lacto_rep_seqs} \
            --output-path {output.output_path}
        
        mv {output.output_path} {output.lacto_fasta}
        """