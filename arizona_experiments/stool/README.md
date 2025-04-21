# University of Arizona Mouse Experiment: Stool Samples

Holds all of the general alpha/beta diversity and taxonomic relative abundance analysis done before stool results needed to be combined with cecal results due to lack of sampling in the high fiber diets at day 3 (`az_exp_combined/`). Survival results for the mice in these experiments is also included in this directory. 

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 1c-d
-   Supplemental Figure 6a-b

| Figure                   | Associated Scripts            |
|--------------------------|-------------------------------|
| Figure 1c-d              | `main_src/mouse_survival.qmd` |
| Supplemental Figure 6a-b | `main_src/mouse_survival.qmd` |

## Directory Key:

File structure of `stool/` looks like so:

``` bash
stool
├── README.md
├── additional_src
│   ├── alpha_beta_div_plots.Rmd
│   ├── alpha_div_stats.Rmd
│   ├── beta_div_stats.Rmd
│   ├── family_abundance.Rmd
│   ├── histopathology.qmd
│   ├── lacto_food_abun.Rmd
│   ├── micro_abun_correlations.Rmd
│   ├── picrust_plots.Rmd
│   ├── picrust_stats.Rmd
│   ├── resilience_homog.qmd
│   └── tax_barplot.Rmd
├── data
│   ├── lacto_qiime
│   │   ├── filt_lacto_table.qza
│   │   ├── placements.qza
│   │   ├── rem_lacto_table.qza
│   │   ├── tax_filt_actual.qza
│   │   ├── taxonomy.qza
│   │   ├── taxonomy_filtered.qza
│   │   ├── taxonomy_filtered.qzv
│   │   └── tree.qza
│   ├── misc
│   │   ├── aim1a_survival.csv
│   │   ├── combined_metadata.tsv
│   │   ├── dna-sequences.fasta
│   │   ├── euk_filt-mergedDietAim1rep-seqs_051523-Copy1.qza
│   │   ├── euk_filt_mergedDietAim1table_051523-Copy1.qza
│   │   ├── histo_updated_metadata.tsv
│   │   ├── lactoOnlydna-sequences.fasta
│   │   ├── merged_metadata1.tsv
│   │   ├── processed_metadata.tsv
│   │   ├── tss_missing_samples.tsv
│   │   ├── tss_sample_asv.tsv
│   │   ├── tss_seq_depth.tsv
│   │   ├── tss_stat_metadata.tsv
│   │   ├── updated_stool_metadata.csv
│   │   └── updated_stool_metadata.tsv
│   ├── picrust
│   │   └── tss3_meta_contrib.tsv
│   └── qiime
│       ├── core_outputs
│       │   ├── bray_curtis_distance_matrix.qza
│       │   ├── bray_curtis_emperor.qzv
│       │   ├── bray_curtis_pcoa_results.qza
│       │   ├── evenness_vector.qza
│       │   ├── faith_pd.tsv
│       │   ├── faith_pd_vector.qza
│       │   ├── jaccard_distance_matrix.qza
│       │   ├── jaccard_emperor.qzv
│       │   ├── jaccard_pcoa_results.qza
│       │   ├── observed_features_vector.qza
│       │   ├── rarefied_table.qza
│       │   ├── shannon_entropy.tsv
│       │   ├── shannon_vector.qza
│       │   ├── unweighted_unifrac_distance_matrix.qza
│       │   ├── unweighted_unifrac_emperor.qzv
│       │   ├── unweighted_unifrac_pcoa_results.qza
│       │   ├── uw_dist_matrix.tsv
│       │   ├── w_dist_matrix.tsv
│       │   ├── weighted_unifrac_distance_matrix.qza
│       │   ├── weighted_unifrac_emperor.qzv
│       │   └── weighted_unifrac_pcoa_results.qza
│       ├── placements.qza
│       ├── tax_filt_actual.qza
│       ├── taxonomy.qza
│       ├── taxonomy_filtered.qza
│       ├── total_sum_scaling.biom
│       ├── total_sum_scaling.qza
│       ├── total_sum_scaling.tsv
│       └── tree.qza
├── main_src
│   ├── metadata_wrangling.Rmd
│   └── mouse_survival.qmd
├── plots
│   ├── alpha_plots.pdf
│   ├── baiH.pdf
│   ├── baiH_stat_vis.pdf
│   ├── baiI.pdf
│   ├── beta_plots.pdf
│   ├── buk_stat_vis.pdf
│   ├── but_stat_vis.pdf
│   ├── buty_stat_vis.pdf
│   ├── butyrate_kinase.pdf
│   ├── butyryl_coa_transferase.pdf
│   ├── cdiff_rel_abun.pdf
│   ├── diet_baseline_plot.pdf
│   ├── diet_comp_lfc.pdf
│   ├── faith_pd.pdf
│   ├── faith_stats.pdf
│   ├── family_abun1.pdf
│   ├── family_abun1_w_stats.pdf
│   ├── fat_survival_curve.pdf
│   ├── fiber_survival_curve.pdf
│   ├── histo_batch_comp.pdf
│   ├── histo_facility_comp.pdf
│   ├── lacto_food_contam.pdf
│   ├── liveCult_survival.pdf
│   ├── shannon_entropy.pdf
│   ├── survival_curve_all.pdf
│   ├── survival_curve_stats.pdf
│   ├── tax_barplot.pdf
│   ├── tss_unweighted_unifrac.pdf
│   ├── tss_weighted_unifrac.pdf
│   ├── uu_homogeneity.pdf
│   ├── uu_resil_homog_stat_vis.pdf
│   ├── uu_resiliency.pdf
│   ├── w_resil_homog_stat_vis.pdf
│   ├── wu_homog_stats.pdf
│   ├── wu_homogeneity.pdf
│   ├── wu_resil_stats.pdf
│   └── wu_resiliency.pdf
├── stats
│   ├── bile_enzyme_lm.tsv
│   ├── buty_enzyme_lm.tsv
│   ├── faith_diet_results.tsv
│   ├── faith_dunn.tsv
│   ├── faith_total_results.tsv
│   ├── family_abun_dunn.tsv
│   ├── family_abun_lm.tsv
│   ├── histo_batch_dunn.tsv
│   ├── histo_facil_dunn.tsv
│   ├── histo_facil_tTest.tsv
│   ├── poster_family_abun_lm.tsv
│   ├── shannon_diet_results.tsv
│   ├── shannon_dunn.tsv
│   ├── shannon_total_results.tsv
│   ├── survival_hazardRatio_results.tsv
│   ├── uu_homog_dunn.tsv
│   ├── uu_homogeneity.tsv
│   ├── uu_resil_dunn.tsv
│   ├── uu_resiliency.tsv
│   ├── uw_adonis_by_day.tsv
│   ├── uw_adonis_day_nc.tsv
│   ├── uw_adonis_day_ns.tsv
│   ├── uw_adonis_nc_nseq.tsv
│   ├── uw_adonis_results.tsv
│   ├── w_adonis_by_day.tsv
│   ├── w_adonis_day_nc.tsv
│   ├── w_adonis_day_ns.tsv
│   ├── w_adonis_nc_nseq.tsv
│   ├── w_adonis_results.tsv
│   ├── wu_homog_dunn.tsv
│   ├── wu_homogeneity.tsv
│   ├── wu_resil_dunn.tsv
│   └── wu_resiliency.tsv
└── workflow_src
    ├── alpha_div_plots.R
    ├── alpha_div_stats.R
    ├── beta_div_plots.R
    ├── beta_div_stats.R
    ├── bile_acid_plots.R
    ├── butyrate_bile_stats.R
    ├── butyrate_plots.R
    ├── core_metrics_prep.ipynb
    ├── family_abun_plots.R
    ├── family_abun_stats.R
    ├── histopathology.R
    ├── homog_calc.R
    ├── hypoxia_plots.R
    ├── ko_contrib_filter.R
    ├── lacto_qiime.sbatch
    ├── metab.R
    ├── metadata_processing.R
    ├── q_ancombc.sh
    ├── resil_calc.R
    ├── sampleids_and_asvs.R
    ├── seq_depth.R
    ├── seq_depth_calc.R
    ├── survival.R
    ├── total_sum_scaling.R
    ├── toxin.R
    └── toxin_fixing.R

12 directories, 163 files
```

**Additional information on directory contents:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *data:*
    -   *lacto_qiime:* mouse purified diet sequencing data qiime analysis (lactococcus investigation)
    -   *misc:* raw data/metadata files used in the analysis
    -   *picrust:* picrust2 meta contrib file on this dataset
    -   *qiime:* qiime analysis files
        -   *core_outputs:* qiime core metrics analysis output files
-   *main_src:* scripts for the results featured in the paper
-   *plots:* all plots generated from `additional_src` and `main_src` analysis
-   *stats:* all stats generated from `additional_src` and `main_src` analysis
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis