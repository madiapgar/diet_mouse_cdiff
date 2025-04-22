# University of Arizona Mouse Experiment: Cecal Content Samples

Holds all of the general alpha/beta diversity, taxonomic relative abundance, histopathology, hypoxia, and metabolomic (SFCA, toxin, bile acids) analysis done for the cecal content samples collected upon sacrifice at day3.

## Paper Relevance:

Paper figures/tables generated from the contents of this directory:

-   Figure 2a-c
-   Figure 3a-e
-   Figure 4c
-   Figure 5c
-   Supplemental Figure 1a-b
-   Supplemental Figure 2a-c
-   Supplemental Figure 3a

> [!IMPORTANT] 
> Plot names in the table also match the variable names of plots included in the [apppleplots](https://github.com/madiapgar/apppleplots) R package!

| Figure                   | Associated Scripts                                                                                                                                                           | Plot Names                                                                                                                                      |
|----------------------|---------------------------|-----------------------|
| Figure 2a-c              | a: `main_src/toxin.Rmd` <br/> b: `main_src/histopathology.Rmd`<br/> c: `main_src/toxin_meta_histo_comp.Rmd`                                                                  | a: `neat_plot` <br/> b: `all_day_plot` <br/> c: `noDiet_histoCecum_toxStats_plot`                                                               |
| Figure 3a-e              | a: `main_src/bile_acid.Rmd` <br/> b: `main_src/bile_acid.Rmd` <br/> c: `main_src/figure3_redo.Rmd` <br/> d: `main_src/metabolomics.Rmd` <br/> e: `main_src/figure3_redo.Rmd` | a: `dca_sep_plot` <br/> b: `ratio_bile_plot` <br/> c: `panelC_fig3_redoWithlab` <br/> d: `butyrateOnly_plot` <br/> e: `panelE_fig3_redoWithlab` |
| Figure 4c                | c: `main_src/microbiome_cecum_comp.Rmd`                                                                                                                                      | c: `faith_combStats_plot`                                                                                                                       |
| Figure 5c                | c: `main_src/microbiome_cecum_comp.Rmd`                                                                                                                                      | c: `microbeProd_relAbunStats_plot`                                                                                                              |
| Supplemental Figure 1a-b | a: `main_src/histo_hypoxia_comp.Rmd` \*\*CONSIDER CHANGING <br/> b: `main_src/histo_hypoxia_comp.Rmd` \*\*CONSIDER CHANGING                                                  | a: `cecum_plot` <br/> b: `colon_plot`                                                                                                           |
| Supplemental Figure 2a-c | a: `main_src/acetate_propionate_comp_redo.Rmd` <br/> b: `main_src/metabolomics.Rmd` <br/> c: `main_src/acetate_propionate_comp_redo.Rmd`                                     | a: `supp_fig2a` <br/> b: `supp_fig2bc` <br/> c: `supp_fig2bc`                                                                                   |
| Supplemental Figure 3a   | a: `main_src/hypoxia.Rmd`                                                                                                                                                    | a: `hypoxia_cecum_plot`                                                                                                                         |

## Directory Key:

File structure of `cecum/` looks like so:

``` bash
cecum
├── README.md
├── additional_src
│   ├── alpha_div_plots.Rmd
│   ├── alpha_div_stats.Rmd
│   ├── beta_div_plots.Rmd
│   ├── beta_div_stats.Rmd
│   ├── bile_acid_comparisons.Rmd
│   ├── family_abun_plots.Rmd
│   ├── family_abun_stats.Rmd
│   ├── histo_percentages.Rmd
│   ├── hypoxia_ancombc.Rmd
│   ├── hypoxia_comparisons.Rmd
│   ├── metab_toxin_comp.Rmd
│   ├── tax_barplot.Rmd
│   └── toxin_metab_histo_comp.Rmd
├── data
│   ├── cecal_qiime
│   │   ├── placements.qza
│   │   ├── tax_filt_actual.qza
│   │   ├── taxonomy.qza
│   │   ├── taxonomy_filtered.qza
│   │   ├── taxonomy_filtered.qzv
│   │   ├── total_sum_filt_table.qza
│   │   ├── total_sum_rem_table.qza
│   │   └── tree.qza
│   ├── cecal_qiime_upper
│   │   ├── SEQ069
│   │   │   ├── S69_barcodes.txt
│   │   │   ├── cecal_s69_barcodes.txt
│   │   │   └── s69_paired_end_seqs.qza
│   │   ├── SEQ070
│   │   │   ├── S70_barcodes.txt
│   │   │   ├── cecal_s70_barcodes.txt
│   │   │   ├── s70_demux.qza
│   │   │   ├── s70_demux.qzv
│   │   │   ├── s70_demux_details.qza
│   │   │   └── s70_paired_end_seqs.qza
│   │   ├── SEQ071
│   │   │   ├── S71_barcodes.txt
│   │   │   ├── cecal_s71_barcodes.txt
│   │   │   └── s71_paired_end_seqs.qza
│   │   ├── meta2.txt
│   │   └── meta3.txt
│   ├── core_outputs
│   │   ├── bray_curtis_distance_matrix.qza
│   │   ├── bray_curtis_emperor.qzv
│   │   ├── bray_curtis_pcoa_results.qza
│   │   ├── evenness_vector.qza
│   │   ├── faith_pd.tsv
│   │   ├── faith_pd_vector.qza
│   │   ├── jaccard_distance_matrix.qza
│   │   ├── jaccard_emperor.qzv
│   │   ├── jaccard_pcoa_results.qza
│   │   ├── observed_features_vector.qza
│   │   ├── rarefied_table.qza
│   │   ├── shannon_entropy.tsv
│   │   ├── shannon_vector.qza
│   │   ├── unweighted_unifrac_distance_matrix.qza
│   │   ├── unweighted_unifrac_emperor.qzv
│   │   ├── unweighted_unifrac_pcoa_results.qza
│   │   ├── uw_dist_matrix.tsv
│   │   ├── w_dist_matrix.tsv
│   │   ├── weighted_unifrac_distance_matrix.qza
│   │   ├── weighted_unifrac_emperor.qzv
│   │   └── weighted_unifrac_pcoa_results.qza
│   ├── misc
│   │   ├── bile_acid.txt
│   │   ├── cecal_key.txt
│   │   ├── cecal_metadata.tsv
│   │   ├── cecal_processed_metadata.tsv
│   │   ├── cecum_histo_percentages.txt
│   │   ├── corrected_bile_acid.tsv
│   │   ├── diet_mouseID_only.tsv
│   │   ├── elenas_diets.tsv
│   │   ├── filt_cecal_metadata.tsv
│   │   ├── filt_cecal_processed_metadata.tsv
│   │   ├── filt_updated_cecal_metadata.tsv
│   │   ├── histoMetabToxin_results_comb.tsv
│   │   ├── histo_categories.txt
│   │   ├── histo_data.csv
│   │   ├── long_cecal_products.tsv
│   │   ├── metabolomics.csv
│   │   ├── mouseID_facil.tsv
│   │   ├── pimid_fluor.csv
│   │   ├── processed_bile_acid.tsv
│   │   ├── processed_dilutedToxin.tsv
│   │   ├── processed_histopathology.tsv
│   │   ├── processed_metabolomics.tsv
│   │   ├── processed_neatToxin.tsv
│   │   ├── processed_ratio_bileAcid.tsv
│   │   ├── seq_depth.tsv
│   │   ├── toxin.csv
│   │   ├── toxin_final_data.tsv
│   │   └── updated_cecal_metadata.tsv
│   └── picrust
│       └── meta_contrib.tsv
├── main_src
│   ├── acetate_propionate_comp_redo.Rmd
│   ├── bile_acid.Rmd
│   ├── bile_acid_comp_redo.Rmd
│   ├── figure3_redo.Rmd
│   ├── histo_hypoxia_comp.Rmd
│   ├── histopathology.Rmd
│   ├── hypoxia.Rmd
│   ├── metabolomics.Rmd
│   ├── microbiome_cecum_comp.Rmd
│   └── toxin.Rmd
├── plots
│   ├── alphaDiv_microbeProd_plot.pdf
│   ├── bileAcid_histoCecum_comp.pdf
│   ├── bileAcid_histoRatio_comp.pdf
│   ├── bileAcid_histo_comp.pdf
│   ├── bileAcid_metabInhib_comp.pdf
│   ├── bileAcid_metabPromot_comp.pdf
│   ├── bileAcid_metab_comp.pdf
│   ├── bileAcid_toxinInhib_comp.pdf
│   ├── bileAcid_toxinPromot_comp.pdf
│   ├── bileAcid_toxin_comp.pdf
│   ├── bileRatio_histo_comp.pdf
│   ├── bileRatio_toxin_comp.pdf
│   ├── bileRatio_toxin_comp_wStats.pdf
│   ├── bile_acid.pdf
│   ├── bile_acid_ratio.pdf
│   ├── butyrateOnly_metab.pdf
│   ├── butyrate_toxHisto_plot.pdf
│   ├── cecal_histo.pdf
│   ├── cecum_histoPerc.pdf
│   ├── colon_histo.pdf
│   ├── dcaSep_overall_bileAcid.pdf
│   ├── dca_toxHisto_plot.pdf
│   ├── dil_histo_toxin_comp.pdf
│   ├── dil_histo_toxin_stats.pdf
│   ├── dil_toxin.pdf
│   ├── faith_butyrate_comp.pdf
│   ├── faith_dca_comp.pdf
│   ├── faith_pd.pdf
│   ├── faith_stat_vis.pdf
│   ├── famAbun_stat_vis.pdf
│   ├── family_abun1.pdf
│   ├── family_abun2.pdf
│   ├── histo_categories.pdf
│   ├── histo_categories_hypox.pdf
│   ├── histo_metab_comp.pdf
│   ├── histo_metab_stats.pdf
│   ├── histopathology.pdf
│   ├── hypoxia_cecum_asvLevel.pdf
│   ├── hypoxia_cecum_famLevel.pdf
│   ├── hypoxia_diet_cecumONLY.pdf
│   ├── hypoxia_diet_facetLocation.pdf
│   ├── hypoxia_fiber_facetLocation.pdf
│   ├── hypoxia_histo_comp_plots.pdf
│   ├── hypoxia_location_facetDiet.pdf
│   ├── hypoxia_metab_comp_plots.pdf
│   ├── hypoxia_toxin_comp_plots.pdf
│   ├── metab_cecumHisto_stats.pdf
│   ├── metab_colonHisto_stats.pdf
│   ├── metab_tox_comp.pdf
│   ├── metabolomics.pdf
│   ├── neat_cecumHisto_stats.pdf
│   ├── neat_colonHisto_stats.pdf
│   ├── neat_histo_toxin_comp.pdf
│   ├── neat_histo_toxin_stats.pdf
│   ├── neat_toxin.pdf
│   ├── noButyrate_metab.pdf
│   ├── overall_bile_acid.pdf
│   ├── relAbun_microbeProd_plot.pdf
│   ├── shannon_entropy.pdf
│   ├── shannon_stat_vis.pdf
│   ├── supp_acetProp_comp_plot.pdf
│   ├── supp_bileAcid_comp_plot.pdf
│   ├── tax_barplot.pdf
│   ├── toxin_acetateStats_plot.pdf
│   ├── toxin_butyrateStats_plot.pdf
│   ├── toxin_propionateStats_plot.pdf
│   ├── unweighted_unifrac_pcoa.pdf
│   └── weighted_unifrac_pcoa.pdf
├── stats
│   ├── bile_acid_dunn.tsv
│   ├── cecum_histoPerc_dunn.tsv
│   ├── dilToxin_dunn_test.tsv
│   ├── dilToxin_kruskal_test.tsv
│   ├── faith_diet_results.tsv
│   ├── faith_dunn.tsv
│   ├── family_abun_dunn.tsv
│   ├── family_abun_lm.tsv
│   ├── histo_categories_hypox.tsv
│   ├── histopathology_dunn.tsv
│   ├── histopathology_lm.tsv
│   ├── hypoxia_dietLocation_dunn.tsv
│   ├── hypoxia_diet_dunn.tsv
│   ├── hypoxia_fiberContent_dunn.tsv
│   ├── hypoxia_fiberLocation_dunn.tsv
│   ├── metab_dunn_test.tsv
│   ├── metab_kruskal_test.tsv
│   ├── metab_linear_model.tsv
│   ├── neatToxin_dunn_test.tsv
│   ├── neatToxin_kruskal_test.tsv
│   ├── overall_bile_dunn.tsv
│   ├── ratio_bile_acid.tsv
│   ├── shannon_diet_results.tsv
│   ├── shannon_dunn.tsv
│   ├── uw_adonis_results.tsv
│   └── w_adonis_results.tsv
└── workflow_src
    ├── alpha_div_plots.R
    ├── alpha_div_stats.R
    ├── beta_div_plots.R
    ├── beta_div_stats.R
    ├── bile_acid_processing.R
    ├── family_abun_plots.R
    ├── family_abun_stats.R
    ├── histopathology.R
    ├── hypoxia_plots.R
    ├── metab.R
    ├── metadata_processing.R
    ├── seq_depth.R
    ├── total_sum_scaling.R
    ├── toxin.R
    └── toxinMetab_histoBile_file_prep.R

15 directories, 205 files
```

**Additional information on directory contents:**

-   *additional_src:* all the rest of the analysis that was done for this dataset that wasn't included in the paper
-   *data:*
    -   *cecal_qiime:* qiime analysis outputs for cecal samples
    -   *core_outputs:* qiime core metrics analysis outputs
    -   *misc:* raw data/metadata files used in the analysis
    -   *picrust:* picrust2 meta contrib file on this dataset
-   *main_src:* scripts for the results featured in the paper
-   *plots:* all plots generated from `additional_src` and `main_src` analysis
-   *stats:* all stats generated from `additional_src` and `main_src` analysis
-   *workflow_src:* R scripts used in the workflow pipeline to generate basic plots/stats for downstream analysis