## step 6
## plots and stats in R for non-longitudinal datasets (measurements only taken at ONE time point)


## UPDATE BILE ACID AND HYPOXIA RULES!!

rule sequencing_depth_calculation:
    input:
        biom = DATASET_DIR + "data/qiime/filt_table.qza",
        asvs = DATASET_DIR + "data/qiime/lactoOnly_rep_seqs.fasta"
    output:
        table = DATASET_DIR + "data/misc/seq_depth.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}seq_depth.R --biom {input.biom} \
                                                    --sequence {input.asvs} \
                                                    --output {output.table}
        """

rule metadata_processing:
    input:
        metadata = DATASET_DIR + METADATA,
        sampleID_key = DATASET_DIR + SAMPLEID_KEY,
        seq_depth = DATASET_DIR + "data/misc/seq_depth.tsv",
        id_file = DATASET_DIR + MOUSEID_FACIL_KEY
    output:
        table = DATASET_DIR + PROCESSED_META
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}metadata_processing.R --metadata {input.metadata} \
                                                                --cecal_key {input.sampleID_key} \
                                                                --seq_depth {input.seq_depth} \
                                                                --mouse_id_facil {input.id_file} \
                                                                --output {output.table}
        """


rule bile_acid_preProcessing:
    input:
        bile_acid = DATASET_DIR + BILE_ACID,
        sampleID_key = DATASET_DIR + SAMPLEID_KEY,
        id_file = DATASET_DIR + MOUSEID_FACIL_KEY
    output:
        proc_bile_acid = DATASET_DIR + "data/misc/corrected_bile_acid.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}bile_acid_processing.R --bile_acid {input.bile_acid} \
                                                                --cecal_key {input.sampleID_key} \
                                                                --mouse_id_facil {input.id_file} \
                                                                --output {output.proc_bile_acid}
        """


rule toxin_metab_histo_bileAcid_processing:
    input:
        metadata = DATASET_DIR + PROCESSED_META,
        toxin = DATASET_DIR + TOXIN,
        histo = DATASET_DIR + HISTO,
        metab = DATASET_DIR + METAB,
        bile_acid = DATASET_DIR + "data/misc/corrected_bile_acid.tsv"
    output:
        proc_neat_toxin = DATASET_DIR + "data/misc/processed_neatToxin.tsv",
        proc_dil_toxin = DATASET_DIR + "data/misc/processed_dilutedToxin.tsv",
        proc_metab = DATASET_DIR + "data/misc/processed_metabolomics.tsv",
        proc_histo = DATASET_DIR + "data/misc/processed_histopathology.tsv",
        proc_bile_acid = DATASET_DIR + "data/misc/processed_bile_acid.tsv",
        proc_bile_ratio = DATASET_DIR + "data/misc/processed_ratio_bileAcid.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}toxinMetab_histoBile_file_prep.R --metadata {input.metadata} \
                                                                            --toxin {input.toxin} \
                                                                            --histopathology {input.histo} \
                                                                            --metabolomics {input.metab} \
                                                                            --bile_acid {input.bile_acid} \
                                                                            --neat_toxin_out {output.proc_neat_toxin} \
                                                                            --dil_toxin_out {output.proc_dil_toxin} \
                                                                            --metab_out {output.proc_metab} \
                                                                            --histo_out {output.proc_histo} \
                                                                            --bile_acid_out {output.proc_bile_acid} \
                                                                            --bile_ratio_out {output.proc_bile_ratio}
        """


rule alpha_diversity_plots:
    input:
        metadata = DATASET_DIR + PROCESSED_META,
        faith_pd = DATASET_DIR + "data/core_outputs/faith_pd.tsv",
        shannon = DATASET_DIR + "data/core_outputs/shannon_entropy.tsv"
    output:
        faith_plot = DATASET_DIR + "plots/faith_pd.pdf",
        shannon_plot = DATASET_DIR + "plots/shannon_entropy.pdf"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
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
        metadata = DATASET_DIR + PROCESSED_META,
        faith_pd = DATASET_DIR + "data/core_outputs/faith_pd.tsv",
        shannon = DATASET_DIR + "data/core_outputs/shannon_entropy.tsv"
    output:
        faith_lm_sec = DATASET_DIR + "stats/faith_diet_results.tsv",
        faith_dunn = DATASET_DIR + "stats/faith_dunn.tsv",
        shannon_lm_sec = DATASET_DIR + "stats/shannon_diet_results.tsv",
        shannon_dunn = DATASET_DIR + "stats/shannon_dunn.tsv",
        faith_plot = DATASET_DIR + "plots/faith_stat_vis.pdf",
        shannon_plot = DATASET_DIR + "plots/shannon_stat_vis.pdf"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}alpha_div_stats.R --metadata {input.metadata} \
                                                            --faith_pd {input.faith_pd} \
                                                            --shannon {input.shannon} \
                                                            --faith_lm_sec {output.faith_lm_sec} \
                                                            --faith_dunn {output.faith_dunn} \
                                                            --shannon_lm_sec {output.shannon_lm_sec} \
                                                            --shannon_dunn {output.shannon_dunn} \
                                                            --faith_plot {output.faith_plot} \
                                                            --shannon_plot {output.shannon_plot}
        """


rule beta_diversity_plots:
    input:
        metadata = DATASET_DIR + PROCESSED_META,
        unweighted_uni = DATASET_DIR + "data/core_outputs/unweighted_unifrac_pcoa_results.qza",
        weighted_uni = DATASET_DIR + "data/core_outputs/weighted_unifrac_pcoa_results.qza",
        faith_pd = DATASET_DIR + "data/core_outputs/faith_pd.tsv",
        shannon = DATASET_DIR + "data/core_outputs/shannon_entropy.tsv"
    output:
        unweighted_plot = DATASET_DIR + "plots/unweighted_unifrac_pcoa.pdf",
        weighted_plot = DATASET_DIR + "plots/weighted_unifrac_pcoa.pdf"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
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
        metadata = DATASET_DIR + PROCESSED_META,
        uw_dist = DATASET_DIR + "data/core_outputs/uw_dist_matrix.tsv",
        w_dist = DATASET_DIR + "data/core_outputs/w_dist_matrix.tsv"
    output:
        w_adonis = DATASET_DIR + "stats/w_adonis_results.tsv",
        uw_adonis = DATASET_DIR + "stats/uw_adonis_results.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}beta_div_stats.R --metadata {input.metadata} \
                                                            --uu_dist {input.uw_dist} \
                                                            --wu_dist {input.w_dist} \
                                                            --uu_adonis {output.uw_adonis} \
                                                            --wu_adonis {output.w_adonis}
        """


rule family_abundance_plots:
    input:
        otu_table = DATASET_DIR + "data/qiime/otu_table.qza",
        taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza",
        metadata = DATASET_DIR + PROCESSED_META
    output:
        plot1 = DATASET_DIR + "plots/family_abun1.pdf",
        plot2 = DATASET_DIR + "plots/family_abun2.pdf"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}family_abun_plots.R --metadata {input.metadata} \
                                                            --otu {input.otu_table} \
                                                            --taxonomy {input.taxonomy} \
                                                            --plot1 {output.plot1} \
                                                            --plot2 {output.plot2}
        """


rule family_abundance_stats:
    input:
       otu_table = DATASET_DIR + "data/qiime/otu_table.qza",
       taxonomy = DATASET_DIR + "data/qiime/taxonomy.qza",
       metadata = DATASET_DIR + PROCESSED_META
    output:
        lm = DATASET_DIR + "stats/family_abun_lm.tsv",
        dunn = DATASET_DIR + "stats/family_abun_dunn.tsv",
        stat_plot = DATASET_DIR + "plots/famAbun_stat_vis.pdf"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}family_abun_stats.R --metadata {input.metadata} \
                                                            --otu {input.otu_table} \
                                                            --taxonomy {input.taxonomy} \
                                                            --stat_plot {output.stat_plot} \
                                                            --linear_model {output.lm} \
                                                            --dunn {output.dunn}
        """


rule histopathology:
    input:
        histo = DATASET_DIR + "data/misc/processed_histopathology.tsv"
    output:
        plot = DATASET_DIR + "plots/histopathology.pdf",
        lm = DATASET_DIR + "stats/histopathology_lm.tsv",
        dunn = DATASET_DIR + "stats/histopathology_dunn.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}histopathology.R --histo {input.histo} \
                                                            --histo_plot {output.plot} \
                                                            --linear_model {output.lm} \
                                                            --dunn {output.dunn}
        """


rule toxin:
    input:
        neat_toxin = DATASET_DIR + "data/misc/processed_neatToxin.tsv",
        dil_toxin = DATASET_DIR + "data/misc/processed_dilutedToxin.tsv"
    output:
        neat_plot = DATASET_DIR + "plots/neat_toxin.pdf",
        diluted_plot = DATASET_DIR + "plots/dil_toxin.pdf",
        neat_kruskal = DATASET_DIR + "stats/neatToxin_kruskal_test.tsv",
        neat_dunn = DATASET_DIR + "stats/neatToxin_dunn_test.tsv",
        dil_kruskal = DATASET_DIR + "stats/dilToxin_kruskal_test.tsv",
        dil_dunn = DATASET_DIR + "stats/dilToxin_dunn_test.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}toxin.R --neat_toxin {input.neat_toxin} \
                                                --dil_toxin {input.dil_toxin} \
                                                --neat_plot {output.neat_plot} \
                                                --diluted_plot {output.diluted_plot} \
                                                --neat_kruskal {output.neat_kruskal} \
                                                --neat_dunn {output.neat_dunn} \
                                                --diluted_kruskal {output.dil_kruskal} \
                                                --diluted_dunn {output.dil_dunn}
        """


rule metabolomics:
    input:
        metab = DATASET_DIR + "data/misc/processed_metabolomics.tsv"
    output:
        metab_plot = DATASET_DIR + "plots/metabolomics.pdf",
        metab_lm = DATASET_DIR + "stats/metab_linear_model.tsv",
        metab_dunn = DATASET_DIR + "stats/metab_dunn_test.tsv",
        metab_kruskal = DATASET_DIR + "stats/metab_kruskal_test.tsv"
    conda:
        "r_env"
    params:
        script_location=DATASET_DIR + R_SCRIPT_DIR
    shell:
        """
        Rscript {params.script_location}metab.R --metab {input.metab} \
                                                --metab_plot {output.metab_plot} \
                                                --metab_lm {output.metab_lm} \
                                                --metab_dunn {output.metab_dunn} \
                                                --metab_kruskal {output.metab_kruskal}
        """