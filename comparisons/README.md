# All Experiment Comparison Analysis: Stool Samples

Holds all of the comparisons between stool 16S sequencing outputs at days -15 and 3 for all experiments done with this mouse model in the Lozupone lab. Experiments included are:
- **AMC-P:** [Previously published experiments](https://www.nature.com/articles/s41522-022-00276-1) done at AMC
- **UA:** Experiments conducted at the University of Arizona
- **AMC-FU:** Follow-up experiments conducted at AMC

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 7a-b

> [!IMPORTANT] 
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                            | Plot Names                                                       |
|--------------------|----------------------------------|------------------|
| Figure 7a-b              | a: `main_src/relAbun_plots.Rmd` <br/> b: `main_src/extra_relAbun_stats.Rmd` | a: `miniVendor_genusAbun_plot` <br/> b: `mini_d3_genus_relAbun_stat_plot`                 |

## Directory Key:

**File structure of `comparisons/` looks like so:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *baseline:* plots and stats generated from comparing baseline (day -15) stool samples for all experiments
    -   *plots:* all plots generated from `additional_src` and `main_src` analysis
    -   *stats:* all stats generated from `additional_src` and `main_src` analysis 
-   *baseline_day3:* plots and stats generated from comparing baseline (day -15) and day 3 stool samples for all experiments 
    -   *plots:* all plots generated from `additional_src` and `main_src` analysis
    -   *stats:* all stats generated from `additional_src` and `main_src` analysis 
-   *data:*
    -   *baseline_day3_qiime:* qiime analysis outputs for baseline (day -15) and day 3 for all experiments 
        -   *core_outputs:* qiime core metrics analysis output files 
    -   *baseline_qiime:* qiime analysis outputs for only baseline (day -15) for all experiments
        -   *core_outputs:* qiime core metrics analysis output files  
    -   *misc:* raw data/metadata files used in the analysis
-   *main_src:* scripts for the results featured in the paper
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis