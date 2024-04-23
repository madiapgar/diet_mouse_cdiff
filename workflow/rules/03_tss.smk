## step 3
## total sum scaling and re-filtration of taxonomy

## uses same biom table and rep seqs inputs as step 2!!

rule total_sum_scaling:
    input:
        biom = DATASET_DIR + BIOM,
        asvs = DATASET_DIR + "data/qiime/lactoOnly_rep_seqs.fasta"
    output:
        table = DATASET_DIR + "data/qiime/total_sum_scaling.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}total_sum_scaling.R --biom {input.biom} \
                                                            --sequence {input.asvs} \
                                                            --output {output.table}
        """


rule tss_tsv2biom:
    input:
        DATASET_DIR + "data/qiime/total_sum_scaling.tsv" 
    output:
        DATASET_DIR + "data/qiime/total_sum_scaling.biom"
    conda:
        "qiime2-2023.5"
    shell:
        """
        biom convert \
            -i {input} \
            -o {output} \
            --table-type "Table" \
            --to-hdf5
        """

rule tss_biom2qza:
    input:
        DATASET_DIR + "data/qiime/total_sum_scaling.biom"
    output:
        DATASET_DIR + "data/qiime/total_sum_scaling.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools import \
            --input-path {input} \
            --type 'FeatureTable[Frequency]' \
            --output-path {output}
        """


rule rep_seqs2fasta:
    input:
        DATASET_DIR + REP_SEQS
    output:
        DATASET_DIR + "data/qiime/fasta_files/dna-sequences.fasta"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {output}
        """


rule sepp_ASV_filtering2:
    input:
        table = DATASET_DIR + "data/qiime/total_sum_scaling.qza",
        tree = DATASET_DIR + "data/qiime/tree.qza"
    output:
        filt_table = DATASET_DIR + "data/qiime/total_sum_filt_table.qza",
        rem_table = DATASET_DIR + "data/qiime/total_sum_rem_table.qza"
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


rule filter_taxonomy2:
    input:
        table = DATASET_DIR + "data/qiime/total_sum_scaling.qza",
        taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza"
    output:
        tax_filt = DATASET_DIR + "data/qiime/tax_filt.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-include p_ \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table {output.tax_filt}
        """


rule filter_taxonomy_vis2:
    input:
        DATASET_DIR + "data/qiime/tax_filt.qza"
    output:
        DATASET_DIR + "data/qiime/tax_filt.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table summarize \
            --i-table {input} \
            --o-visualization {output}
        """