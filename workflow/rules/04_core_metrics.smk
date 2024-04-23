## step 4
## core metrics analysis


rule pre_core_metrics_filter:
    input:
        tax_filt = DATASET_DIR + "data/qiime/tax_filt.qza",
        metadata = DATASET_DIR + METADATA
    output:
        otu_table = DATASET_DIR + "data/qiime/otu_table.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.otu_table}
        """


rule core_metrics_analysis:
    input:
        tree = DATASET_DIR + "data/qiime/tree.qza",
        otu_table = DATASET_DIR + "data/qiime/otu_table.qza",
        metadata = DATASET_DIR + METADATA
    output:
        output_dir = directory(DATASET_DIR + "data/core_outputs")
    conda:
        "qiime2-2023.5"
    params:
        sampling_depth=CORE_SAMPLING_DEPTH
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.tree} \
            --i-table {input.meta_filt_tax} \
            --p-sampling-depth  {params.sampling_depth} \
            --m-metadata-file {input.metadata} \
            --output-dir {output.output_dir}
        """

rule unzip_uw_distance_matrix:
    input:
        DATASET_DIR + "data/core_outputs/unweighted_unifrac_distance_matrix.qza" 
    output:
        DATASET_DIR + "data/core_outputs/uw_dist_matrix.tsv"
    conda:
        "qiime2-2023.5"
    params:
        location=DATASET_DIR
    shell:
        """
        qiime tools export \
            --input-path ./{params.location}data/core_outputs/unweighted_unifrac_distance_matrix.qza \
            --output-path ./{params.location}data/core_outputs/uw_dist_matrix
        
        mv ./{params.location}data/core_outputs/uw_dist_matrix/distance-matrix.tsv \
        ./{params.location}data/core_outputs/uw_dist_matrix/uw_dist_matrix.tsv

        mv ./{params.location}data/core_outputs/uw_dist_matrix/uw_dist_matrix.tsv \
        ./{params.location}data/core_outputs/
        """


rule unzip_w_distance_matrix:
    input:
        DATASET_DIR + "data/core_outputs/weighted_unifrac_distance_matrix.qza"
    output:
        DATASET_DIR + "data/core_outputs/w_dist_matrix.tsv"
    conda:
        "qiime2-2023.5"
    params:
        location=DATASET_DIR
    shell:
        """
        qiime tools export \
            --input-path ./{params.location}data/core_outputs/weighted_unifrac_distance_matrix.qza \
            --output-path ./{params.location}data/core_outputs/w_dist_matrix
        
        mv ./{params.location}data/core_outputs/w_dist_matrix/distance-matrix.tsv \
        ./{params.location}data/core_outputs/w_dist_matrix/w_dist_matrix.tsv

        mv ./{params.location}data/core_outputs/w_dist_matrix/w_dist_matrix.tsv \
        ./{params.location}data/core_outputs/ 
        """


rule unzip_shannon:
    input:
        DATASET_DIR + "data/core_outputs/shannon_vector.qza"
    output:
        DATASET_DIR + "data/core_outputs/shannon_entropy.tsv"
    conda:
        "qiime2-2023.5"
    params:
        location=DATASET_DIR
    shell:
        """
        qiime tools export \
            --input-path ./{params.location}data/core_outputs/shannon_vector.qza \
            --output-path ./{params.location}data/core_outputs/shannon_entropy
        
        mv ./{params.location}data/core_outputs/shannon_entropy/alpha-diversity.tsv \
        ./{params.location}data/core_outputs/shannon_entropy/shannon_entropy.tsv

        mv ./{params.location}data/core_outputs/shannon_entropy/shannon_entropy.tsv \
        ./{params.location}data/core_outputs/
        """


rule unzip_faith_pd:
    input:
        DATASET_DIR + "data/core_outputs/faith_pd_vector.qza"
    output:
        DATASET_DIR + "data/core_outputs/faith_pd.tsv"
    conda:
        "qiime2-2023.5"
    params:
        location=DATASET_DIR
    shell:
        """
        qiime tools export \
            --input-path ./{params.location}data/core_outputs/faith_pd_vector.qza \
            --output-path ./{params.location}data/core_outputs/faith_pd
        
        mv ./{params.location}data/core_outputs/faith_pd/alpha-diversity.tsv \
        ./{params.location}data/core_outputs/faith_pd/faith_pd.tsv

        mv ./{params.location}data/core_outputs/faith_pd/faith_pd.tsv \
        ./{params.location}data/core_outputs/
        """