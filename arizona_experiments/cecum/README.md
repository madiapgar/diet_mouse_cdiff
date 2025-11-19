# University of Arizona Mouse Experiment: Cecal Content Samples

Holds all of the general 16S sequencing alpha/beta diversity, taxonomic relative abundance, histopathology, hypoxia, and metabolomic (SFCA, toxin, bile acids) analysis done for the cecal content samples collected upon sacrifice at day3.

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 2a-c
-   Figure 3a-e
-   Figure 4c
-   Figure 5c
-   Supplemental Figure 1a-b
-   Supplemental Figure 2a-b

> [!IMPORTANT] 
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                                                                                                                           | Plot Names                                                                                                                                      |
|----------------------|---------------------------|-----------------------|
| Figure 2a-c              | a: `main_src/toxin.Rmd` <br/> b: `main_src/histopathology.Rmd`<br/> c: `main_src/toxin_metab_histo_comp.Rmd`                                                                  | a: `neat_plot` <br/> b: `all_day_plot` <br/> c: `noDiet_histoCecum_toxStats_plot`                                                               |
| Figure 3a-e              | a: `main_src/bile_acid.Rmd` <br/> b: `main_src/bile_acid.Rmd` <br/> c: `main_src/figure3_redo.Rmd` <br/> d: `main_src/metabolomics.Rmd` <br/> e: `main_src/acetate_propionate_comp_redo.Rmd` | a: `dca_sep_plot` <br/> b: `ratio_bile_plot` <br/> c: `panelC_fig3_redoWithlab` <br/> d: `metab_plot` <br/> e: `all_metab_plotsWithLab` |
| Figure 4c                | c: `main_src/microbiome_cecum_comp.Rmd`                                                                                                                                      | c: `faith_combStats_plot`                                                                                                                       |
| Figure 5c                | c: `main_src/microbiome_cecum_comp.Rmd`                                                                                                                                      | c: `microbeProd_relAbunStats_plot`                                                                                                              |
| Supplemental Figure 1a-b | a: `main_src/histoScore_breakdown_plots.Rmd` <br/> b: `main_src/histoScore_breakdown_plots.Rmd`                                                 | a: `cecum_plot` <br/> b: `colon_plot`                                                                                                           |
| Supplemental Figure 2a-b | a: `main_src/bile_acid_comp_redo.Rmd` <br/> b: `main_src/mediation.Rmd`                                     | a: `supp_fig2a` <br/> b: `hflf_buty_plot` |

## Directory Key:

**File structure of `cecum/` looks like so:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *data:*
    -   *cecal_qiime:* qiime analysis outputs for cecal samples
    -   *core_outputs:* qiime core metrics analysis outputs
    -   *misc:* raw data/metadata files used in the analysis
-   *main_src:* scripts for the results featured in the paper
-   *plots:* all plots generated from `additional_src` and `main_src` analysis
-   *stats:* all stats generated from `additional_src` and `main_src` analysis
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis