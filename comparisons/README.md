# All Experiment Comparison Analysis: Stool Samples

Holds all of the comparisons between stool 16S sequencing outputs at days -15 and 3 for all experiments done with this mouse model in the Lozupone lab. Experiments included are:
- **AMC-P:** [Previously published experiments](https://www.nature.com/articles/s41522-022-00276-1#ref-CR67) done at AMC
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

**Directory Tree:**

``` bash
comparisons
├── README.md
├── additional_src
│   ├── az_newExp_diet_deathDate.Rmd
│   └── tax_biplot.Rmd
├── baseline
│   ├── plots
│   │   ├── faith_pd.pdf
│   │   ├── faith_stat_vis.pdf
│   │   ├── famAbun_stat_vis.pdf
│   │   ├── genusAbun_stat_vis.pdf
│   │   ├── relAbun_by_family1.pdf
│   │   ├── relAbun_by_family2.pdf
│   │   ├── relAbun_by_genus1.pdf
│   │   ├── relAbun_by_genus2.pdf
│   │   ├── shannon_entropy.pdf
│   │   ├── shannon_stat_vis.pdf
│   │   ├── unweighted_unifrac_pcoa.pdf
│   │   ├── uu_biplot1_2.pdf
│   │   ├── vendor_uu_pcoa.pdf
│   │   ├── vendor_wu_pcoa.pdf
│   │   ├── weighted_unifrac_pcoa.pdf
│   │   └── wu_biplot1_2.pdf
│   └── stats
│       ├── faith_dunn.tsv
│       ├── faith_lm.tsv
│       ├── family_abun_dunn.tsv
│       ├── family_abun_lm.tsv
│       ├── genus_abun_dunn.tsv
│       ├── genus_abun_lm.tsv
│       ├── shannon_dunn.tsv
│       ├── shannon_lm.tsv
│       ├── uw_adonis_results.tsv
│       └── w_adonis_results.tsv
├── baseline_day3
│   ├── plots
│   │   ├── deathDate_stat_plot.pdf
│   │   ├── diet_famAbun_plot.pdf
│   │   ├── diet_genusAbun_plot.pdf
│   │   ├── famAbun_stat_vis.pdf
│   │   ├── genusAbun_stat_vis.pdf
│   │   ├── miniDiet_genusAbun_plot.pdf
│   │   ├── miniVendor_genusAbun_plot.pdf
│   │   ├── mini_genusAbun_stat_vis.pdf
│   │   ├── vendor_famAbun_plot.pdf
│   │   └── vendor_genusAbun_plot.pdf
│   └── stats
│       └── allSurv_deathDate_dunn.tsv
├── data
│   ├── baseline_day3_qiime
│   │   ├── allExp_combined_d15-d3_seqs.qza
│   │   ├── allExp_combined_d15-d3_table.qza
│   │   ├── core_outputs
│   │   │   ├── bray_curtis_distance_matrix.qza
│   │   │   ├── bray_curtis_emperor.qzv
│   │   │   ├── bray_curtis_pcoa_results.qza
│   │   │   ├── evenness_vector.qza
│   │   │   ├── faith_pd.tsv
│   │   │   ├── faith_pd_vector.qza
│   │   │   ├── jaccard_distance_matrix.qza
│   │   │   ├── jaccard_emperor.qzv
│   │   │   ├── jaccard_pcoa_results.qza
│   │   │   ├── observed_features_vector.qza
│   │   │   ├── rarefied_table.qza
│   │   │   ├── shannon_entropy.tsv
│   │   │   ├── shannon_vector.qza
│   │   │   ├── unweighted_unifrac_distance_matrix.qza
│   │   │   ├── unweighted_unifrac_emperor.qzv
│   │   │   ├── unweighted_unifrac_pcoa_results.qza
│   │   │   ├── uw_dist_matrix.tsv
│   │   │   ├── w_dist_matrix.tsv
│   │   │   ├── weighted_unifrac_distance_matrix.qza
│   │   │   ├── weighted_unifrac_emperor.qzv
│   │   │   └── weighted_unifrac_pcoa_results.qza
│   │   ├── taxonomy.qza
│   │   └── total_sum_otu_table.qza
│   ├── baseline_qiime
│   │   ├── allExp_comp_d15_seqs.qza
│   │   ├── allExp_comp_d15_table.qza
│   │   ├── core_outputs
│   │   │   ├── bray_curtis_distance_matrix.qza
│   │   │   ├── bray_curtis_emperor.qzv
│   │   │   ├── bray_curtis_pcoa_results.qza
│   │   │   ├── evenness_vector.qza
│   │   │   ├── faith_pd.tsv
│   │   │   ├── faith_pd_vector.qza
│   │   │   ├── jaccard_distance_matrix.qza
│   │   │   ├── jaccard_emperor.qzv
│   │   │   ├── jaccard_pcoa_results.qza
│   │   │   ├── observed_features_vector.qza
│   │   │   ├── rarefied_table.qza
│   │   │   ├── relative_rarefied_table.qza
│   │   │   ├── shannon_entropy.tsv
│   │   │   ├── shannon_vector.qza
│   │   │   ├── unweighted_unifrac_distance_matrix.qza
│   │   │   ├── unweighted_unifrac_emperor.qzv
│   │   │   ├── unweighted_unifrac_pcoa_results.qza
│   │   │   ├── uu_biplot_matrix.qza
│   │   │   ├── uw_dist_matrix.tsv
│   │   │   ├── w_dist_matrix.tsv
│   │   │   ├── weighted_unifrac_distance_matrix.qza
│   │   │   ├── weighted_unifrac_emperor.qzv
│   │   │   ├── weighted_unifrac_pcoa_results.qza
│   │   │   └── wu_biplot_matrix.qza
│   │   ├── taxonomy.qza
│   │   └── total_sum_otu_table.qza
│   └── misc
│       ├── mouse_CF_set_5_mapping.txt
│       ├── mouse_set_3_mapping.txt
│       ├── mouse_set_4_mapping.txt
│       ├── newExp_comp_d15-d3_metadata.tsv
│       ├── newExp_comp_d15_metadata.tsv
│       ├── oldNew_comp_d15-d3_metadata.tsv
│       ├── oldNew_comp_d15_metadata.tsv
│       └── sample_ids.csv
├── main_src
│   ├── extra_relAbun_stats.Rmd
│   └── relAbun_plots.Rmd
└── workflow_src
    ├── alpha_div_plots.R
    ├── alpha_div_stats.R
    ├── beta_div_plots.R
    ├── beta_div_stats.R
    ├── family_abun_plots.R
    ├── family_abun_stats.R
    ├── make_paired_end_obj.sh
    ├── oldNew_meta_comb.R
    ├── q_make_biplot.sh
    ├── q_merge_tables.sh
    ├── seq016_demux.sh
    ├── seq_depth.R
    └── total_sum_scaling.R
```