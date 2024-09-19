## step 3
## total sum scaling and re-filtration of taxonomy

## uses same biom table and rep seqs inputs as step 2!!
import os

rule create_lacto_table:
    input:
        filt_table = os.path.join(DATASET_DIR, "data/qiime/filt_table.qza"),
        taxonomy = os.path.join(DATASET_DIR, "data/qiime/taxonomy.qza"),
        seqs = os.path.join(DATASET_DIR, REP_SEQS)
    output:
        lacto_table = os.path.join(DATASET_DIR, "data/qiime/lacto_cecal_table.qza"),
        lacto_rep_seqs = os.path.join(DATASET_DIR, "data/qiime/lacto_rep_seqs.qza")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/qiime/lacto_rep_seqs.qza")
    output:
        os.path.join(DATASET_DIR, "data/qiime/lactoOnly_rep_seqs.fasta")
    conda:
        QIIME
    params:
        location=DATASET_DIR
    shell:
        """
        qiime tools export \
            --input-path ./{params.location}data/qiime/lacto_rep_seqs.qza \
            --output-path ./{params.location}data/qiime/fasta_files
        
        mv ./{params.location}data/qiime/fasta_files/dna-sequences.fasta \
        ./{params.location}data/qiime/fasta_files/lactoOnly_rep_seqs.fasta

        mv ./{params.location}data/qiime/fasta_files/lactoOnly_rep_seqs.fasta \
        ./{params.location}data/qiime/lactoOnly_rep_seqs.fasta
        """


rule total_sum_scaling:
    input:
        biom = os.path.join(DATASET_DIR, BIOM),
        asvs = os.path.join(DATASET_DIR, "data/qiime/lactoOnly_rep_seqs.fasta")
    output:
        table = os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.tsv")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}total_sum_scaling.R --biom {input.biom} \
                                                            --sequence {input.asvs} \
                                                            --output {output.table}
        """


rule tss_tsv2biom:
    input:
        os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.tsv")
    output:
        os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.biom")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.biom")
    output:
        os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.qza")
    conda:
        QIIME
    shell:
        """
        qiime tools import \
            --input-path {input} \
            --type 'FeatureTable[Frequency]' \
            --output-path {output}
        """


rule sepp_ASV_filtering2:
    input:
        table = os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.qza"),
        tree = os.path.join(DATASET_DIR, "data/qiime/tree.qza")
    output:
        filt_table = os.path.join(DATASET_DIR, "data/qiime/total_sum_filt_table.qza"),
        rem_table = os.path.join(DATASET_DIR, "data/qiime/total_sum_rem_table.qza")
    conda:
        QIIME
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
        table = os.path.join(DATASET_DIR, "data/qiime/total_sum_scaling.qza"),
        taxonomy = os.path.join(DATASET_DIR, "data/qiime/taxonomy.qza")
    output:
        tax_filt = os.path.join(DATASET_DIR, "data/qiime/total_sum_tax_filt.qza")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/qiime/total_sum_tax_filt.qza")
    output:
        os.path.join(DATASET_DIR, "data/qiime/total_sum_tax_filt.qzv")
    conda:
        QIIME
    shell:
        """
        qiime feature-table summarize \
            --i-table {input} \
            --o-visualization {output}
        """


rule pre_core_metrics_filter:
    input:
        tax_filt = os.path.join(DATASET_DIR, "data/qiime/total_sum_tax_filt.qza"),
        metadata = os.path.join(DATASET_DIR, METADATA)
    output:
        otu_table = os.path.join(DATASET_DIR, "data/qiime/total_sum_otu_table.qza")
    conda:
        QIIME
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.otu_table}
        """