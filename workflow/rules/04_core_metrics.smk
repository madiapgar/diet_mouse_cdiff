## step 4
## core metrics analysis
import os


rule pre_core_metrics_filter:
    input:
        tax_filt = os.path.join(DATASET_DIR, "data/qiime/tax_filt.qza"),
        metadata = os.path.join(DATASET_DIR, METADATA)
    output:
        otu_table = os.path.join(DATASET_DIR, "data/qiime/otu_table.qza")
    conda:
        QIIME
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.otu_table}
        """


rule core_metrics_analysis:
    input:
        tree = os.path.join(DATASET_DIR, "data/qiime/tree.qza"),
        otu_table = os.path.join(DATASET_DIR, "data/qiime/otu_table.qza"),
        metadata = os.path.join(DATASET_DIR, METADATA)
    output:
        output_dir = directory(os.path.join(DATASET_DIR, "data/core_outputs"))
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/core_outputs/unweighted_unifrac_distance_matrix.qza")
    output:
        os.path.join(DATASET_DIR, "data/core_outputs/uw_dist_matrix.tsv")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/core_outputs/weighted_unifrac_distance_matrix.qza")
    output:
        os.path.join(DATASET_DIR, "data/core_outputs/w_dist_matrix.tsv")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/core_outputs/shannon_vector.qza")
    output:
        os.path.join(DATASET_DIR, "data/core_outputs/shannon_entropy.tsv")
    conda:
        QIIME
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
        os.path.join(DATASET_DIR, "data/core_outputs/faith_pd_vector.qza")
    output:
        os.path.join(DATASET_DIR, "data/core_outputs/faith_pd.tsv")
    conda:
        QIIME
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