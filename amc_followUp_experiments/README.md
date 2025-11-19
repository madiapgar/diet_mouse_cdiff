# AMC Follow-Up Mouse Experiments: Stool and Culture Samples

Holds all of the stool 16S sequencing outputs at days -15 and 3 along with the spleen, liver, and blood culture 16S sequencing outputs. Includes survival analysis and mouse weight data over the timeline. All data was collected from follow-up mouse experiments done at AMC (Anschutz Medical Campus). 

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 6a-b
-   Figure 7a-b
-   Supplemental Figure 6b-c

> [!IMPORTANT] 
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                            | Plot Names                                                       |
|--------------------|----------------------------------|------------------|
| Figure 6a-b              | a: `main_src/tax_biplot.Rmd` <br/> b: `main_src/cfu_surv_weight_comp.Rmd` | a: `uu_biplot`, `uu_biplot2_3` <br/> b: `facet_status_cfu_plot` |
| Figure 7a-b              | a: `main_src/plasma_cytokine_pcoa.Rmd` <br/> b: `main_src/plasma_cytokine_conc_plot.Rmd` | a: `p_pcoa_1_2` <br/> b: `cytokine_plot` |
| Supplemental Figure 6b-c | b: `main_src/survival.Rmd` <br/> c: `main_src/mouse_weight.Rmd`               | a: `newExp_dietVendor_surv_plot` <br/> b: `newExp_avWeight_plot` |

## Directory Key:

**File structure of `amc_followUp_experiments/` looks like so:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *data:*
    -   *misc:* raw data/metadata files used in the analysis
    -   *new_culture_qiime:* qiime analysis outputs for spleen, liver, and blood culture sequences 
    -   *qiime:* qiime analysis outputs for stool samples (days -15 and 3)
        - *core_outputs:* qiime core metrics analysis output files
-   *main_src:* scripts for the results featured in the paper
-   *plots:* all plots generated from `additional_src` and `main_src` analysis
-   *stats:* all stats generated from `additional_src` and `main_src` analysis
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis
