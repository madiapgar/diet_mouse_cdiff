# University of Arizona Mouse Experiment: Stool and Cecal Content Samples Combined

Holds all of the general 16S sequencing alpha/beta diversity and taxonomic relative abundance analysis done with the stool and cecal content samples combined (cecal content samples were substituted for day 3 stool samples due to undersampling of the high fiber diets) and batch 1 of the experiments filtered out since those mice did not get infected with *C. difficile*. Preliminary analysis was done with the day 3 stool and cecal content samples together and it was determined that the cecal content samples were not significantly different at day 3 and could be substituted to increase power.

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 4a-b
-   Figure 5a-b
-   Supplemental Figure 4a-g
-   Supplemental Figure 5

> [!IMPORTANT]
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                                                                                                                                                                                                                                           | Plot Names                                                                                                                                                                               |
|--------------------|-----------------------------------|-----------------|
| Figure 4a-b              | a: `main_src/alpha_div_plots.Rmd` <br/> b: `main_src/alpha_div_stats.Rmd`                                                                                                                                                                                                                    | a: `faith_plot` <br/> b: `faith_stat_vis`                                                                                                                                                |
| Figure 5a-b              | a: `main_src/figure5_redo.Rmd` <br/> b: `main_src/figure5_redo.Rmd`                                                                                                                                                                                                                          | a: `lachno_plot`/`rumino_plot` <br/> b: `path_abun_plot` <br/> together: `genusAbun_plots_together_withLabs`                                                                                                                                |
| Supplemental Figure 4a-g | a: `main_src/lactococcus_relAbun.Rmd` <br/> b: `main_src/lactococcus_relAbun.Rmd` <br/> c: `main_src/beta_div_plots.Rmd` <br/> d: `main_src/resilience_homog.qmd` <br/> e: `main_src/resilience_homog.qmd` <br/> f: `main_src/resilience_homog.qmd` <br/> g: `main_src/resilience_homog.qmd` | a: `lacto_contam_plot` <br/> b: `lacto_stat_plot` <br/> c: `unweighted_pcoa` <br/> d: `uu_resil_plot` <br/> e: `uu_homog_plot` <br/> f:  `uu_resil_stat_vis` <br/> g: `uu_homog_stat_vis` |
| Supplemental Figure 5    | c: `main_src/figure5_redo.Rmd`                                                                                                                                                                                                                                                               | c: `abun_stat_plot`                                                                                                                                                                      |

## Directory Key:

**File structure of `az_exp_combined/` looks like so:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *data:*
    -   *comp_qiime:* qiime analysis outputs for stool and cecal samples combined
    -   *misc:* raw data/metadata files used in the analysis
    -   *s1_filt_core:* qiime core metrics analysis outputs (first batch is filtered out since mice were not infected with C. diff)
-   *main_src:* scripts for the results featured in the paper
-   *plots:* all plots generated from `additional_src` and `main_src` analysis
-   *stats:* all stats generated from `additional_src` and `main_src` analysis
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis

**Directory Tree:**

``` bash
az_exp_combined
├── README.md
├── additional_src
│   ├── dist_comparisons.Rmd
│   ├── family_abundance.Rmd
│   ├── micro_abun_correlations.Rmd
│   ├── sample_matrix.Rmd
│   ├── tax_barplot.Rmd
│   └── timeline_fig.Rmd
├── data
│   ├── comp_qiime
│   │   ├── tax_s1_filt.qza
│   │   ├── taxonomy.qza
│   │   └── taxonomy_filtered.qza
│   ├── misc
│   │   ├── comp_metadata.tsv
│   │   ├── comp_rep_seqs.qza
│   │   ├── comp_table.qza
│   │   └── s1_filt_comp_metadata.tsv
│   └── s1_filt_core
│       ├── bray_curtis_distance_matrix.qza
│       ├── bray_curtis_emperor.qzv
│       ├── bray_curtis_pcoa_results.qza
│       ├── evenness_vector.qza
│       ├── faith_pd.tsv
│       ├── faith_pd_vector.qza
│       ├── jaccard_distance_matrix.qza
│       ├── jaccard_emperor.qzv
│       ├── jaccard_pcoa_results.qza
│       ├── observed_features_vector.qza
│       ├── rarefied_table.qza
│       ├── shannon_entropy.tsv
│       ├── shannon_vector.qza
│       ├── unweighted_unifrac_distance_matrix.qza
│       ├── unweighted_unifrac_emperor.qzv
│       ├── unweighted_unifrac_pcoa_results.qza
│       ├── uw_dist_matrix.tsv
│       ├── w_dist_matrix.tsv
│       ├── weighted_unifrac_distance_matrix.qza
│       ├── weighted_unifrac_emperor.qzv
│       └── weighted_unifrac_pcoa_results.qza
├── main_src
│   ├── alpha_div_plots.Rmd
│   ├── alpha_div_stats.Rmd
│   ├── beta_div_plots.Rmd
│   ├── beta_div_stats.Rmd
│   ├── figure5_redo.Rmd
│   ├── lactococcus_relAbun.Rmd
│   └── resilience_homog.qmd
├── plots
│   ├── cdiff_rel_abun.pdf
│   ├── d3_filt_uu_pcoa.pdf
│   ├── d3_filt_wu_pcoa.pdf
│   ├── diet_baseline_plot.pdf
│   ├── diet_composition_lfc.pdf
│   ├── facAnaerobe_family_abun.pdf
│   ├── facility_uu_pcoa.pdf
│   ├── facility_wu_pcoa.pdf
│   ├── faith_pd.pdf
│   ├── faith_stats.pdf
│   ├── famAbun_stat_vis.pdf
│   ├── famAbun_together.pdf
│   ├── family_abun1.pdf
│   ├── genus_abun1.pdf
│   ├── lachno_genera_by_day.pdf
│   ├── lacto_relAbun_diet.pdf
│   ├── lacto_relAbun_withStats.pdf
│   ├── lacto_relbun_stats.pdf
│   ├── long_family_abun1.pdf
│   ├── long_family_abun1_w_stats.pdf
│   ├── mini_famAbun_stat_vis.pdf
│   ├── mini_family_abun1.pdf
│   ├── new_fig5_stats.pdf
│   ├── nonFilt_tax_barplot.pdf
│   ├── obAnaerobe_family_abun.pdf
│   ├── potential_new_fig5.pdf
│   ├── rumino_genera_by_day.pdf
│   ├── sample_spec_table.pdf
│   ├── shannon_entropy.pdf
│   ├── shannon_stats.pdf
│   ├── tax_barplot.pdf
│   ├── timeline_plot.pdf
│   ├── top_lachno_rumino.pdf
│   ├── unweighted_unifrac_pcoa.pdf
│   ├── updated_timeline.pdf
│   ├── uu_homog_stats.pdf
│   ├── uu_homogeneity.pdf
│   ├── uu_resil_homog_stat_vis.pdf
│   ├── uu_resil_stats.pdf
│   ├── uu_resiliency.pdf
│   ├── w_resil_homog_stat_vis.pdf
│   ├── weighted_unifrac_pcoa.pdf
│   ├── wu_homog_stats.pdf
│   ├── wu_homogeneity.pdf
│   ├── wu_resil_stats.pdf
│   └── wu_resiliency.pdf
├── stats
│   ├── faith_diet_results.tsv
│   ├── faith_total_results.tsv
│   ├── shannon_diet_results.tsv
│   ├── shannon_total_results.tsv
│   ├── short_family_abun_lm.tsv
│   ├── uu_homog_dunn.tsv
│   ├── uu_homogeneity.tsv
│   ├── uu_resil_dunn.tsv
│   ├── uu_resiliency.tsv
│   ├── uw_adonis_by_day.tsv
│   ├── uw_adonis_results.tsv
│   ├── w_adonis_by_day.tsv
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
    ├── comparison.sbatch
    ├── day3_filt.sbatch
    ├── family_abun_plots.R
    ├── family_abun_stats.R
    ├── homog_calc.R
    ├── mouse_count.R
    ├── q_merge_tables.sh
    ├── resil_calc.R
    └── total_sum_scaling.R
```