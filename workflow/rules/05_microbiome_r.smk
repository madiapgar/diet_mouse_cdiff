## step 5
## plots and stats in R for longitudinal datasets (measurements taken at multiple time points)

import os
import pandas as pd

rule alpha_diversity_plots:
    input:
        metadata = os.path.join(DATASET_DIR, PROCESSED_META),
        faith_pd = os.path.join(DATASET_DIR, "data/qiime/core_outputs/faith_pd.tsv"),
        shannon = os.path.join(DATASET_DIR, "data/qiime/core_outputs/shannon_entropy.tsv")
    output:
        faith_plot = os.path.join(DATASET_DIR, "plots/faith_pd.pdf"),
        shannon_plot = os.path.join(DATASET_DIR, "plots/shannon_entropy.pdf")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}alpha_div_plots.R --metadata {input.metadata} \
                                                            --faith_pd {input.faith_pd} \
                                                            --shannon {input.shannon} \
                                                            --output_faith {output.faith_plot} \
                                                            --output_shannon {output.shannon_plot}
        """


rule alpha_diversity_stats:
    input:
        metadata = os.path.join(DATASET_DIR, PROCESSED_META),
        faith_pd = os.path.join(DATASET_DIR, "data/qiime/core_outputs/faith_pd.tsv"),
        shannon = os.path.join(DATASET_DIR, "data/qiime/core_outputs/shannon_entropy.tsv")
    output:
        faith_lm = os.path.join(DATASET_DIR, "stats/faith_lm.tsv"),
        faith_dunn = os.path.join(DATASET_DIR, "stats/faith_dunn.tsv"),
        shannon_lm = os.path.join(DATASET_DIR, "stats/shannon_lm.tsv"),
        shannon_dunn = os.path.join(DATASET_DIR, "stats/shannon_dunn.tsv"),
        faith_plot = os.path.join(DATASET_DIR, "plots/faith_stat_vis.pdf"),
        shannon_plot = os.path.join(DATASET_DIR, "plots/shannon_stat_vis.pdf")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}alpha_div_stats.R --metadata {input.metadata} \
                                                            --faith_pd {input.faith_pd} \
                                                            --shannon {input.shannon} \
                                                            --faith_lm {output.faith_lm} \
                                                            --faith_dunn {output.faith_dunn} \
                                                            --shannon_lm {output.shannon_lm} \
                                                            --shannon_dunn {output.shannon_dunn} \
                                                            --faith_plot {output.faith_plot} \
                                                            --shannon_plot {output.shannon_plot}
        """
        


rule beta_diversity_plots:
    input:
        metadata = os.path.join(DATASET_DIR, PROCESSED_META),
        unweighted_uni = os.path.join(DATASET_DIR, "data/qiime/core_outputs/unweighted_unifrac_pcoa_results.qza"),
        weighted_uni = os.path.join(DATASET_DIR, "data/qiime/core_outputs/weighted_unifrac_pcoa_results.qza"),
        faith_pd = os.path.join(DATASET_DIR, "data/qiime/core_outputs/faith_pd.tsv"),
        shannon = os.path.join(DATASET_DIR, "data/qiime/core_outputs/shannon_entropy.tsv")
    output:
        unweighted_plot = os.path.join(DATASET_DIR, "plots/unweighted_unifrac_pcoa.pdf"),
        weighted_plot = os.path.join(DATASET_DIR, "plots/weighted_unifrac_pcoa.pdf")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}beta_div_plots.R --metadata {input.metadata} \
                                                            --unweighted_unifrac {input.unweighted_uni} \
                                                            --weighted_unifrac {input.weighted_uni} \
                                                            --faith_pd {input.faith_pd} \
                                                            --shannon {input.shannon} \
                                                            --output_uu {output.unweighted_plot} \
                                                            --output_wu {output.weighted_plot}
            """


rule beta_diversity_stats:
    input:
        metadata = os.path.join(DATASET_DIR, PROCESSED_META),
        uw_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/uw_dist_matrix.tsv"),
        w_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/w_dist_matrix.tsv")
    output:
        w_adonis = os.path.join(DATASET_DIR, "stats/w_adonis_results.tsv"),
        uw_adonis = os.path.join(DATASET_DIR, "stats/uw_adonis_results.tsv")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}beta_div_stats.R --metadata {input.metadata} \
                                                         --uu_dist {input.uw_dist} \
                                                         --wu_dist {input.w_dist} \
                                                         --uu_adonis {output.uw_adonis} \
                                                         --wu_adonis {output.w_adonis}
        """


rule relative_abundance_plots:
    input:
        otu_table = os.path.join(DATASET_DIR, WHICH_OTU),
        taxonomy = os.path.join(DATASET_DIR, "data/qiime/taxonomy.qza"),
        metadata = os.path.join(DATASET_DIR, PROCESSED_META)
    output:
        plot1 = os.path.join(DATASET_DIR, "plots/rel_abun1.pdf"),
        plot2 = os.path.join(DATASET_DIR, "plots/rel_abun2.pdf")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}family_abun_plots.R --metadata {input.metadata} \
                                                            --otu {input.otu_table} \
                                                            --taxonomy {input.taxonomy} \
                                                            --plot1 {output.plot1} \
                                                            --plot2 {output.plot2}
        """


rule relative_abundance_stats:
    input:
       otu_table = os.path.join(DATASET_DIR, WHICH_OTU),
       taxonomy = os.path.join(DATASET_DIR, "data/qiime/taxonomy.qza"),
       metadata = os.path.join(DATASET_DIR, PROCESSED_META)
    output:
        lm = os.path.join(DATASET_DIR, "stats/rel_abun_lm.tsv"),
        dunn = os.path.join(DATASET_DIR, "stats/rel_abun_dunn.tsv"),
        stat_plot = os.path.join(DATASET_DIR, "plots/relAbun_stat_vis.pdf")
    conda:
        "r_env"
    params:
        script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
    shell:
        """
        Rscript {params.script_location}family_abun_stats.R --metadata {input.metadata} \
                                                            --otu {input.otu_table} \
                                                            --taxonomy {input.taxonomy} \
                                                            --stat_plot {output.stat_plot} \
                                                            --linear_model {output.lm} \
                                                            --dunn {output.dunn}
        """


## conditionally running resiliency and homogeneity rules depending on whether -8 timepoint
## is present in the metadata
proc_meta = pd.read_csv(os.path.join(DATASET_DIR, PROCESSED_META), sep = '\t')

if proc_meta.query('day_post_inf == -8').shape[0] > 0 == True:
    rule homogeneity:
        input:
            metadata = os.path.join(DATASET_DIR, PROCESSED_META),
            uu_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/uw_dist_matrix.tsv"),
            wu_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/w_dist_matrix.tsv")
        output:
            wu_lm = os.path.join(DATASET_DIR, "stats/wu_homogeneity.tsv"),
            wu_dunn = os.path.join(DATASET_DIR, "stats/wu_homog_dunn.tsv"),
            uu_lm = os.path.join(DATASET_DIR, "stats/uu_homogeneity.tsv"),
            uu_dunn = os.path.join(DATASET_DIR, "stats/uu_homog_dunn.tsv"),
            wu_plot = os.path.join(DATASET_DIR, "plots/wu_homogeneity.pdf"),
            wu_stat_plot = os.path.join(DATASET_DIR, "plots/wu_homog_stats.pdf"),
            uu_plot = os.path.join(DATASET_DIR, "plots/uu_homogeneity.pdf"),
            uu_stat_plot = os.path.join(DATASET_DIR, "plots/uu_homog_stats.pdf")
        conda:
            "r_env"
        params:
            script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
        shell:
            """
            Rscript {params.script_location}homog_calc.R --metadata {input.metadata} \
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
            metadata = os.path.join(DATASET_DIR, PROCESSED_META),
            uu_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/uw_dist_matrix.tsv"),
            wu_dist = os.path.join(DATASET_DIR, "data/qiime/core_outputs/w_dist_matrix.tsv")
        output:
            uu_lm = os.path.join(DATASET_DIR, "stats/uu_resiliency.tsv"),
            uu_dunn = os.path.join(DATASET_DIR, "stats/uu_resil_dunn.tsv"),
            wu_lm = os.path.join(DATASET_DIR, "stats/wu_resiliency.tsv"),
            wu_dunn = os.path.join(DATASET_DIR, "stats/wu_resil_dunn.tsv"),
            wu_plot = os.path.join(DATASET_DIR, "plots/wu_resiliency.pdf"),
            wu_stat_plot = os.path.join(DATASET_DIR, "plots/wu_resil_stats.pdf"),
            uu_plot = os.path.join(DATASET_DIR, "plots/uu_resiliency.pdf"),
            uu_stat_plot = os.path.join(DATASET_DIR, "plots/uu_resil_stats.pdf")
        conda:
            "r_env"
        params:
            script_location=os.path.join(DATASET_DIR, R_SCRIPT_DIR)
        shell:
            """
            Rscript {params.script_location}resil_calc.R --metadata {input.metadata} \
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
else:
    print("Need sufficient time points to run resiliency and homogeneity calculations")
