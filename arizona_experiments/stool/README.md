# University of Arizona Mouse Experiment: Stool Samples

Holds all of the general 16S sequencing alpha/beta diversity and taxonomic relative abundance analysis done before stool results needed to be combined with cecal results due to lack of sampling in the high fiber diets at day 3 (`az_exp_combined/`). Survival results for the mice in these experiments are also included in this directory.

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 1c-d
-   Supplemental Figure 6a-b

> [!IMPORTANT]
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts            | Plot Names                                                      |
|------------------------|------------------------|------------------------|
| Figure 1c-d              | `main_src/mouse_survival.qmd` | c: `diet_plot_final` <br/> d: `surv_stat_vis`                   |
| Supplemental Figure 6a-b | `main_src/mouse_survival.qmd` | a: `liveCult_diet_plot_final` <br/> b: `liveCult_surv_stat_vis` |

## Directory Key:

**File structure of `stool/` looks like so:**

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