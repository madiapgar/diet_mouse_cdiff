# AMC Follow-Up Mouse Experiments: Stool and Culture Samples

Holds all of the stool sequencing outputs at days -15 and 3 along with the spleen, liver, and blood culture sequencing outputs. Includes survival analysis and mouse weight data over the timeline. All data was collected from follow-up mouse experiments done at AMC (Anschutz Medical Campus). 

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 6a-b
-   Supplemental Figure 7b-c

> [!IMPORTANT] 
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                            | Plot Names                                                       |
|--------------------|----------------------------------|------------------|
| Figure 6a-b              | a: `main_src/combined_cfu_plot.Rmd` <br/> b: `main_src/combined_cfu_plot.Rmd` | a: `everything_plot` <br/> b: `everything_stats`                 |
| Supplemental Figure 7b-c | b: `main_src/survival.Rmd` <br/> c: `main_src/mouse_weight.Rmd`               | a: `newExp_dietVendor_surv_plot` <br/> b: `newExp_avWeight_plot` |

## Directory Key:

File structure of `amc_followUp_experiments/` looks like so:

``` bash
amc_followUp_experiments
├── README.md
├── additional_src
│   ├── abun_day_comp.Rmd
│   ├── cdd01_culture_plots.Rmd
│   ├── cdd02-3_culture_plots.Rmd
│   ├── cfu_counts.Rmd
│   ├── kidney_histology.Rmd
│   ├── tax_barplot.Rmd
│   ├── tax_biplot.Rmd
│   ├── tax_biplot.md
│   └── tax_posCult_surv.Rmd
├── data
│   ├── misc
│   │   ├── blood_colony_count.txt
│   │   ├── blood_culture_metadata.tsv
│   │   ├── combined_CDD_cfuCounts.txt
│   │   ├── culture_genusAbun_table.tsv
│   │   ├── histologySubtotalsUnblinded_wKidney.csv
│   │   ├── newExp_cfus.txt
│   │   ├── newExp_d15-d3_metadata.txt
│   │   ├── newExp_d15-d3_seq_depth.tsv
│   │   ├── newExp_mouse_weightData.csv
│   │   ├── pos_culture_status.tsv
│   │   ├── proc_blood_culture_meta.tsv
│   │   ├── proc_cdd02-3_culture_metadata.tsv
│   │   ├── proc_colony_count.tsv
│   │   ├── proc_combined_CDD_cfuCounts.tsv
│   │   ├── proc_newExp_d15-d3_metadata.tsv
│   │   ├── proc_tax_CDD_cfuCounts.tsv
│   │   ├── survival_data.tsv
│   │   └── survival_status.tsv
│   ├── new_culture_qiime
│   │   ├── merged_rep_seqs.qza
│   │   ├── merged_table.qza
│   │   ├── taxOnly_otu_table.qza
│   │   └── taxonomy.qza
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
│       │   ├── relative_rarefied_table.qza
│       │   ├── shannon_entropy.tsv
│       │   ├── shannon_vector.qza
│       │   ├── unweighted_unifrac_distance_matrix.qza
│       │   ├── unweighted_unifrac_emperor.qzv
│       │   ├── unweighted_unifrac_pcoa_results.qza
│       │   ├── uu_biplot_matrix.qza
│       │   ├── uw_dist_matrix.tsv
│       │   ├── w_dist_matrix.tsv
│       │   ├── weighted_unifrac_distance_matrix.qza
│       │   ├── weighted_unifrac_emperor.qzv
│       │   ├── weighted_unifrac_pcoa_results.qza
│       │   └── wu_biplot_matrix.qza
│       ├── newExp_d15-d3_allBatch_seqs.qza
│       ├── newExp_d15-d3_allBatch_table.qza
│       ├── taxonomy.qza
│       └── total_sum_otu_table.qza
├── main_src
│   ├── combined_cfu_plot.Rmd
│   ├── mouse_weight.Rmd
│   └── survival.Rmd
├── plots
│   ├── allExp_culture_futureTest.pdf
│   ├── allExp_culture_results.pdf
│   ├── allExp_culture_withMicrobes.pdf
│   ├── allExp_culture_withMicrobes_stats.pdf
│   ├── avWeight_change_plot.pdf
│   ├── bloodCulture_abun1.pdf
│   ├── bloodCulture_abun2.pdf
│   ├── bloodCulture_tax_barplot.pdf
│   ├── cdd01_culture_results.pdf
│   ├── cdd02-3_cult_cfu_plot.pdf
│   ├── cdd02-3_cult_relAbun_diet_plot.pdf
│   ├── cdd02-3_cult_relAbun_location_plot.pdf
│   ├── cdd02_culture_results.pdf
│   ├── cdd03_culture_results.pdf
│   ├── cfus_bloodOnly.pdf
│   ├── dietExp_survival.pdf
│   ├── dietVendor_survival.pdf
│   ├── experimentVendor_survival.pdf
│   ├── facAn_deltas_plot.pdf
│   ├── faith_pd.pdf
│   ├── faith_stat_vis.pdf
│   ├── famAbun_stat_vis.pdf
│   ├── family_abun1.pdf
│   ├── family_abun2.pdf
│   ├── family_uu_biplot.pdf
│   ├── family_wu_biplot.pdf
│   ├── fullTax_uu_biplot.pdf
│   ├── fullTax_wu_biplot.pdf
│   ├── kidneyHisto_byExp_plot.pdf
│   ├── kidney_histo_plot.pdf
│   ├── mini_bloodCulture_abun.pdf
│   ├── newExp_stool_tax_barplot.pdf
│   ├── obAn_deltas_plot.pdf
│   ├── shannon_entropy.pdf
│   ├── shannon_stat_vis.pdf
│   ├── unweighted_unifrac_pcoa.pdf
│   ├── uu_biplot1_2.pdf
│   ├── uu_biplot2_3.pdf
│   ├── uu_biplot3_4.pdf
│   ├── vendor_uu_pcoa.pdf
│   ├── vendor_wu_pcoa.pdf
│   ├── weighted_unifrac_pcoa.pdf
│   ├── wu_biplot1_2.pdf
│   ├── wu_biplot2_3.pdf
│   └── wu_biplot3_4.pdf
├── stats
│   ├── avWeight_change_dunn.tsv
│   ├── faith_diet_results.tsv
│   ├── faith_dunn.tsv
│   ├── faith_lm.tsv
│   ├── faith_total_results.tsv
│   ├── family_abun_dunn.tsv
│   ├── family_abun_lm.tsv
│   ├── shannon_diet_results.tsv
│   ├── shannon_dunn.tsv
│   ├── shannon_lm.tsv
│   ├── shannon_total_results.tsv
│   ├── surv_diet.tsv
│   ├── surv_dietExp.tsv
│   ├── surv_dietVendor.tsv
│   ├── surv_dietVendor_dunn.tsv
│   ├── surv_expVendor.tsv
│   ├── surv_hazardRatio_results.tsv
│   ├── uw_adonis_by_day.tsv
│   ├── uw_adonis_results.tsv
│   ├── w_adonis_by_day.tsv
│   └── w_adonis_results.tsv
└── workflow_src
    ├── alpha_div_plots.R
    ├── alpha_div_stats.R
    ├── beta_div_plots.R
    ├── beta_div_stats.R
    ├── bloodCulture_meta_proc.R
    ├── blood_culture.sbatch
    ├── cdd02-3_cultMeta_proc.R
    ├── family_abun_plots.R
    ├── family_abun_stats.R
    ├── homog_calc.R
    ├── linux_meta_proc.sh
    ├── make_paired_end_obj.sh
    ├── newExp_d15-d3.sbatch
    ├── newExp_meta_proc.R
    ├── q_make_biplot.sh
    ├── q_newExp_d15-d3.sh
    ├── resil_calc.R
    ├── seq_depth.R
    ├── survival_data_prep.R
    └── total_sum_scaling.R
```

**Additional information on directory contents:**

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