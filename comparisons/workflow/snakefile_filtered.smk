rule all:
    input:
        "data/comp_qiime/tax_s1_filt.qza",
        "data/s1_filt_core",
        "data/s1_filt_core/uw_dist_matrix.tsv",
        "data/s1_filt_core/w_dist_matrix.tsv",
        "data/s1_filt_core/shannon_entropy.tsv",
        "data/s1_filt_core/faith_pd.tsv",
        "s1_filt_plots/faith_pd.pdf",
        "s1_filt_plots/shannon_entropy.pdf",
        "s1_filt_stats/faith_total_results.tsv",
        "s1_filt_stats/faith_diet_results.tsv",
        "s1_filt_stats/faith_dunn.tsv",
        "s1_filt_stats/shannon_total_results.tsv",
        "s1_filt_stats/shannon_diet_results.tsv",
        "s1_filt_stats/shannon_dunn.tsv",
        "s1_filt_plots/faith_stat_vis.pdf",
        "s1_filt_plots/shannon_stat_vis.pdf",
        "s1_filt_plots/unweighted_unifrac_pcoa.pdf",
        "s1_filt_plots/weighted_unifrac_pcoa.pdf",
        "s1_filt_stats/w_adonis_results.tsv",
        "s1_filt_stats/uw_adonis_results.tsv",
        "s1_filt_stats/w_adonis_by_day.tsv",
        "s1_filt_stats/uw_adonis_by_day.tsv",
        "s1_filt_stats/wu_homogeneity.tsv",
        "s1_filt_stats/wu_homog_dunn.tsv",
        "s1_filt_stats/uu_homogeneity.tsv",
        "s1_filt_stats/uu_homog_dunn.tsv",
        "s1_filt_plots/wu_homogeneity.pdf",
        "s1_filt_plots/wu_homog_stats.pdf",
        "s1_filt_plots/uu_homogeneity.pdf",
        "s1_filt_plots/uu_homog_stats.pdf",
        "s1_filt_stats/uu_resiliency.tsv",
        "s1_filt_stats/uu_resil_dunn.tsv",
        "s1_filt_stats/wu_resiliency.tsv",
        "s1_filt_stats/wu_resil_dunn.tsv",
        "s1_filt_plots/wu_resiliency.pdf",
        "s1_filt_plots/wu_resil_stats.pdf",
        "s1_filt_plots/uu_resiliency.pdf",
        "s1_filt_plots/uu_resil_stats.pdf",
        "s1_filt_plots/family_abun1.pdf",
        "s1_filt_plots/family_abun2.pdf",
        "s1_filt_stats/family_abun_lm.tsv",
        "s1_filt_stats/family_abun_dunn.tsv",
        "s1_filt_plots/famAbun_stat_vis.pdf"



rule pre_core_metrics_filter:
    input:
        tss_tax_filt = "data/comp_qiime/tss_tax_filt.qza",
        metadata = "data/misc/s1_filt_comp_metadata.tsv"
    output:
        meta_filt_tax = "data/comp_qiime/tax_s1_filt.qza"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime feature-table filter-samples \
            --i-table {input.tss_tax_filt} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.meta_filt_tax}
        """


rule core_metrics_analysis:
##  NEED TO CHECK THE SAMPLING DEPTH BEFORE YOU RUN THIS STEP
    input:
        tree = "data/comp_qiime/tree.qza",
        meta_filt_tax = "data/comp_qiime/tax_s1_filt.qza",
        metadata = "data/misc/s1_filt_comp_metadata.tsv"
    output:
        output_dir = directory("data/s1_filt_core")
    conda:
        "qiime2-2023.5"
    params:
        sampling_depth=99631
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.tree} \
            --i-table {input.meta_filt_tax} \
            --p-sampling-depth {params.sampling_depth} \
            --m-metadata-file {input.metadata} \
            --output-dir {output}
        """


rule unzip_uw_distance_matrix:
    input:
        "data/s1_filt_core/unweighted_unifrac_distance_matrix.qza" 
    output:
        "data/s1_filt_core/uw_dist_matrix.tsv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path ./data/s1_filt_core/unweighted_unifrac_distance_matrix.qza \
            --output-path ./data/s1_filt_core/uw_dist_matrix
        
        mv ./data/s1_filt_core/uw_dist_matrix/distance-matrix.tsv \
        ./data/s1_filt_core/uw_dist_matrix/uw_dist_matrix.tsv

        mv ./data/s1_filt_core/uw_dist_matrix/uw_dist_matrix.tsv \
        ./data/s1_filt_core/
        """


rule unzip_w_distance_matrix:
    input:
        "data/s1_filt_core/weighted_unifrac_distance_matrix.qza"
    output:
        "data/s1_filt_core/w_dist_matrix.tsv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path ./data/s1_filt_core/weighted_unifrac_distance_matrix.qza \
            --output-path ./data/s1_filt_core/w_dist_matrix
        
        mv ./data/s1_filt_core/w_dist_matrix/distance-matrix.tsv \
        ./data/s1_filt_core/w_dist_matrix/w_dist_matrix.tsv

        mv ./data/s1_filt_core/w_dist_matrix/w_dist_matrix.tsv \
        ./data/s1_filt_core/ 
        """


rule unzip_shannon:
    input:
        "data/s1_filt_core/shannon_vector.qza"
    output:
        "data/s1_filt_core/shannon_entropy.tsv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path ./data/s1_filt_core/shannon_vector.qza \
            --output-path ./data/s1_filt_core/shannon_entropy
        
        mv ./data/s1_filt_core/shannon_entropy/alpha-diversity.tsv \
        ./data/s1_filt_core/shannon_entropy/shannon_entropy.tsv

        mv ./data/s1_filt_core/shannon_entropy/shannon_entropy.tsv \
        ./data/s1_filt_core/
        """


rule unzip_faith_pd:
    input:
        "data/s1_filt_core/faith_pd_vector.qza"
    output:
        "data/s1_filt_core/faith_pd.tsv"
    conda:
        "qiime2-2023.5"
    shell:
        """
        qiime tools export \
            --input-path ./data/s1_filt_core/faith_pd_vector.qza \
            --output-path ./data/s1_filt_core/faith_pd
        
        mv ./data/s1_filt_core/faith_pd/alpha-diversity.tsv \
        ./data/s1_filt_core/faith_pd/faith_pd.tsv

        mv ./data/s1_filt_core/faith_pd/faith_pd.tsv \
        ./data/s1_filt_core/
        """


rule alpha_diversity_plots:
    input:
        metadata = "data/misc/s1_filt_comp_metadata.tsv",
        faith_pd = "data/s1_filt_core/faith_pd.tsv",
        shannon = "data/s1_filt_core/shannon_entropy.tsv"
    output:
        faith_plot = "s1_filt_plots/faith_pd.pdf",
        shannon_plot = "s1_filt_plots/shannon_entropy.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/alpha_div_plots.R --metadata {input.metadata} \
                                          --faith_pd {input.faith_pd} \
                                          --shannon {input.shannon} \
                                          --output_faith {output.faith_plot} \
                                          --output_shannon {output.shannon_plot}
        """


rule alpha_diversity_stats:
    input:
        metadata = "data/misc/s1_filt_comp_metadata.tsv",
        faith_pd = "data/s1_filt_core/faith_pd.tsv",
        shannon = "data/s1_filt_core/shannon_entropy.tsv"
    output:
        faith_lm = "s1_filt_stats/faith_total_results.tsv",
        faith_lm_sec = "s1_filt_stats/faith_diet_results.tsv",
        faith_dunn = "s1_filt_stats/faith_dunn.tsv",
        shannon_lm = "s1_filt_stats/shannon_total_results.tsv",
        shannon_lm_sec = "s1_filt_stats/shannon_diet_results.tsv",
        shannon_dunn = "s1_filt_stats/shannon_dunn.tsv",
        faith_plot = "s1_filt_plots/faith_stat_vis.pdf",
        shannon_plot = "s1_filt_plots/shannon_stat_vis.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/alpha_div_stats.R --metadata {input.metadata} \
                                          --faith_pd {input.faith_pd} \
                                          --shannon {input.shannon} \
                                          --faith_lm {output.faith_lm} \
                                          --faith_lm_sec {output.faith_lm_sec} \
                                          --faith_dunn {output.faith_dunn} \
                                          --shannon_lm {output.shannon_lm} \
                                          --shannon_lm_sec {output.shannon_lm_sec} \
                                          --shannon_dunn {output.shannon_dunn} \
                                          --faith_plot {output.faith_plot} \
                                          --shannon_plot {output.shannon_plot}
        """


rule beta_diversity_plots:
    input:
        metadata = "data/misc/s1_filt_comp_metadata.tsv",
        unweighted_uni = "data/s1_filt_core/unweighted_unifrac_pcoa_results.qza",
        weighted_uni = "data/s1_filt_core/weighted_unifrac_pcoa_results.qza",
        faith_pd = "data/s1_filt_core/faith_pd.tsv",
        shannon = "data/s1_filt_core/shannon_entropy.tsv"
    output:
        unweighted_plot = "s1_filt_plots/unweighted_unifrac_pcoa.pdf",
        weighted_plot = "s1_filt_plots/weighted_unifrac_pcoa.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/beta_div_plots.R --metadata {input.metadata} \
                                         --unweighted_unifrac {input.unweighted_uni} \
                                         --weighted_unifrac {input.weighted_uni} \
                                         --faith_pd {input.faith_pd} \
                                         --shannon {input.shannon} \
                                         --output_uu {output.unweighted_plot} \
                                         --output_wu {output.weighted_plot}
        """


rule beta_diversity_stats:
    input:
        metadata = "data/misc/s1_filt_comp_metadata.tsv",
        uw_dist = "data/s1_filt_core/uw_dist_matrix.tsv",
        w_dist = "data/s1_filt_core/w_dist_matrix.tsv"
    output:
        w_adonis = "s1_filt_stats/w_adonis_results.tsv",
        uw_adonis = "s1_filt_stats/uw_adonis_results.tsv",
        w_adonis_day = "s1_filt_stats/w_adonis_by_day.tsv",
        uw_adonis_day = "s1_filt_stats/uw_adonis_by_day.tsv"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/beta_div_stats.R --metadata {input.metadata} \
                                         --uu_dist {input.uw_dist} \
                                         --wu_dist {input.w_dist} \
                                         --uu_adonis {output.uw_adonis} \
                                         --wu_adonis {output.w_adonis} \
                                         --uu_adonis_day {output.uw_adonis_day} \
                                         --wu_adonis_day {output.w_adonis_day}
        """


rule homogeneity:
    input:
        metadata = "data/misc/s1_filt_comp_metadata.tsv",
        uu_dist = "data/s1_filt_core/uw_dist_matrix.tsv",
        wu_dist = "data/s1_filt_core/w_dist_matrix.tsv"
    output:
        wu_lm = "s1_filt_stats/wu_homogeneity.tsv",
        wu_dunn = "s1_filt_stats/wu_homog_dunn.tsv",
        uu_lm = "s1_filt_stats/uu_homogeneity.tsv",
        uu_dunn = "s1_filt_stats/uu_homog_dunn.tsv",
        wu_plot = "s1_filt_plots/wu_homogeneity.pdf",
        wu_stat_plot = "s1_filt_plots/wu_homog_stats.pdf",
        uu_plot = "s1_filt_plots/uu_homogeneity.pdf",
        uu_stat_plot = "s1_filt_plots/uu_homog_stats.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/homog_calc.R --metadata {input.metadata} \
                                     --uu_dist {input.uu_dist} \
                                     --wu_dist {input.wu_dist} \
                                     --weighted_lm {output.wu_lm} \
                                     --weighted_dunn {output.wu_dunn} \
                                     --unweighted_lm {output.uu_lm} \
                                     --unweighted_dunn {output.uu_dunn} \
                                     --weighted_plot {output.wu_plot} \
                                     --weighted_stat_plot {output.wu_stat_plot} \
                                     --unweighted_plot {output.uu_plot} \
                                     --unweighted_stat_plot {output.uu_stat_plot}
        """


rule resiliency:
    input:
       metadata = "data/misc/s1_filt_comp_metadata.tsv",
       uu_dist = "data/s1_filt_core/uw_dist_matrix.tsv",
       wu_dist = "data/s1_filt_core/w_dist_matrix.tsv"
    output:
        uu_lm = "s1_filt_stats/uu_resiliency.tsv",
        uu_dunn = "s1_filt_stats/uu_resil_dunn.tsv",
        wu_lm = "s1_filt_stats/wu_resiliency.tsv",
        wu_dunn = "s1_filt_stats/wu_resil_dunn.tsv",
        wu_plot = "s1_filt_plots/wu_resiliency.pdf",
        wu_stat_plot = "s1_filt_plots/wu_resil_stats.pdf",
        uu_plot = "s1_filt_plots/uu_resiliency.pdf",
        uu_stat_plot = "s1_filt_plots/uu_resil_stats.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/resil_calc.R --metadata {input.metadata} \
                                     --uu_dist {input.uu_dist} \
                                     --wu_dist {input.wu_dist} \
                                     --weighted_lm {output.wu_lm} \
                                     --weighted_dunn {output.wu_dunn} \
                                     --unweighted_lm {output.uu_lm} \
                                     --unweighted_dunn {output.uu_dunn} \
                                     --weighted_plot {output.wu_plot} \
                                     --weighted_stat_plot {output.wu_stat_plot} \
                                     --unweighted_plot {output.uu_plot} \
                                     --unweighted_stat_plot {output.uu_stat_plot}
        """


rule family_abundance_plots:
    input:
        otu_table = "data/comp_qiime/tax_s1_filt.qza",
        taxonomy = "data/comp_qiime/taxonomy.qza",
        metadata = "data/misc/s1_filt_comp_metadata.tsv"
    output:
        plot1 = "s1_filt_plots/family_abun1.pdf",
        plot2 = "s1_filt_plots/family_abun2.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/family_abun_plots.R --metadata {input.metadata} \
                                            --otu {input.otu_table} \
                                            --taxonomy {input.taxonomy} \
                                            --plot1 {output.plot1} \
                                            --plot2 {output.plot2}
        """


rule family_abundance_stats:
    input:
       otu_table = "data/comp_qiime/tax_s1_filt.qza",
       taxonomy = "data/comp_qiime/taxonomy.qza",
       metadata = "data/misc/s1_filt_comp_metadata.tsv"
    output:
        lm = "s1_filt_stats/family_abun_lm.tsv",
        dunn = "s1_filt_stats/family_abun_dunn.tsv",
        stat_plot = "s1_filt_plots/famAbun_stat_vis.pdf"
    conda:
        "r_env"
    shell:
        """
        Rscript scripts/family_abun_stats.R --metadata {input.metadata} \
                                            --otu {input.otu_table} \
                                            --taxonomy {input.taxonomy} \
                                            --stat_plot {output.stat_plot} \
                                            --linear_model {output.lm} \
                                            --dunn {output.dunn}
        """