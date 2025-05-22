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