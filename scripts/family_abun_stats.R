## 7-14-23
## this script contains the stats for the wanted microbe family abundances 
## by days relative to infection and diet

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(vegan)
library(viridis)
library(rstatix)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-o",
                    "--otu",
                    dest = "otu_table_FP",
                    help = "Filepath to OTU Table file in .qza format.")
parser$add_argument("-t",
                    "--taxonomy",
                    dest = "tax_FP",
                    help = "Filepath to taxonomy file in .qza format.")
parser$add_argument("-s",
                    "--stat_plot",
                    dest = "stat_plot_FP",
                    help = "Filepath to family abundance statistical plot in .pdf format.")
parser$add_argument("-lm",
                    "--linear_model",
                    dest = "lm_FP",
                    help = "Filepath to family abundance linear model results in .tsv format.")
parser$add_argument("-d",
                    "--dunn",
                    dest = "dunn_FP",
                    help = "Filepath to family abundance Dunn's Post Hoc test results in .tsv format.")

args <- parser$parse_args()

## input file paths and others
# otu_table_FP <- './data/qiime/taxonomy_filtered.qza'
# tax_FP <- './data/qiime/taxonomy.qza'
# metadata_FP <- './data/misc/processed_metadata.tsv'
wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Tannerellaceae', 'Muribaculaceae', 'Bacteroidaceae', 
                   'Marinifilaceae', 'Ruminococcaceae')

## function 1
family_abun_file_prep <- function(metadata_fp,
                                  tax_fp,
                                  otu_table_fp,
                                  tax_level,
                                  wanted_tax){
  ## metadata
  metadata <- read_tsv(metadata_FP)
  ## taxonomy
  taxonomy <- read_qza(tax_FP)$data %>% 
    parse_taxonomy() %>% 
    rownames_to_column('asv')
  ## otu table 
  otu_table <- read_qza(otu_table_FP)$data
  otu_table %>% 
    as_tibble(rownames = 'asv') %>% 
    gather(-asv, key = sampleid, value = abun) %>% 
    group_by(sampleid) %>% 
    mutate(rel_abun = abun/sum(abun)) %>% 
    mutate(p_abun = rel_abun + 0.000001) -> otu_table
  ## joining all tables together 
  otu_table %>% 
    left_join(metadata, by = 'sampleid') %>% 
    left_join(taxonomy, by = 'asv') -> abun_table
  abun_table %>% 
    group_by(sampleid, day_post_inf, diet, mouse_id, 
             purified_diet, high_fat, high_fiber, 
             seq_depth, .data[[tax_level]]) %>% 
    summarise(rel_abund = sum(rel_abun)) %>% 
    filter(.data[[tax_level]] %in% wanted_tax) %>% 
    mutate(mouse_fact = as.factor(mouse_id),
           day_fact = as.factor(day_post_inf)) -> abun_filt
  ## creating a list for my outputs
  my_list <- list(Metadata = metadata,
                  Taxonomy = taxonomy,
                  OTUTable = otu_table,
                  AbundanceTable = abun_filt)
  return(my_list)
}

## 2 
## preps dunns post hoc results for statistical visualization
stat_plot_prep <- function(biom_table,
                           dunn_test){
  biom_table %>% 
    group_by(diet, Family, day_post_inf) %>% 
    summarise(mean_abund = mean(rel_abund)) -> mean_table
  dunn_test %>% 
    merge(mean_table, 
          by.x = c('group1',
                   'day_post_inf',
                   'Family'),
          by.y = c('diet',
                   'day_post_inf',
                   'Family')) %>% 
    rename('group1_mean' = 'mean_abund') %>% 
    merge(mean_table,
          by.x = c('group2',
                   'day_post_inf',
                   'Family'),
          by.y = c('diet',
                   'day_post_inf',
                   'Family')) %>% 
    rename('group2_mean' = 'mean_abund') %>% 
    mutate(diff_means = (group1_mean - group2_mean),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  return(new_dunn)
}

## 3 
## statistical visualization 
stat_plot <- function(new_dunn){
  new_dunn %>% 
    filter(day_post_inf != -15) %>%
    ggplot(aes(x = group1, y = group2)) +
    geom_tile(aes(fill = stat_diff_means), alpha = 0.8, color = 'black') +
    scale_fill_gradient2(low = 'blue', high = 'green', name = 'Group 1 -\nGroup 2') +
    geom_text(aes(label = p.adj.signif)) +
    scale_x_discrete(labels = c('Chow',
                                'HFt/\nHFb',
                                'HFt/\nLFb',
                                'LFt/\nHFb')) +
    scale_y_discrete(labels = c('LFt / LFb',
                                'LFt / HFb',
                                'HFt / LFb',
                                'HFt / HFb')) +
    facet_grid(Family~day_post_inf,
               scales = 'free_x') +
    theme_bw(base_size = 16) +
    theme(strip.text.y = element_text(angle = 0)) +
    xlab('Group 1') +
    ylab('Group 2') -> stat_vis
  return(stat_vis)
}

## prepping the file needed to run the linear model
abun_files <- family_abun_file_prep(args$metadata_FP,
                                    args$tax_FP,
                                    args$otu_table_FP,
                                    wanted_level,
                                    wanted_family)

## can still pull out the metadata, taxonomy, and otu table as well
abun_filt <- abun_files$AbundanceTable

## performing the linear model
## looking at each microbe family by days relative to infection and the different diets
## takes out the intercept term and only shows significant p values 
abun_filt %>%
  na.omit() %>% 
  group_by(Family, day_post_inf) %>% 
  do(tidy(lm(rel_abund ~ (purified_diet * seq_depth) + high_fat + high_fiber,
             data =.))) %>% 
  adjust_pvalue(method = 'BH') %>% 
  na.omit() %>% 
  filter(term != '(Intercept)',
         p.value <= 0.05) -> family_abun_lm

## performing Kruskal-Wallis and Dunn's Post Hoc test
abun_filt %>% 
  na.omit() %>% 
  group_by(Family, day_post_inf) %>% 
  do(tidy(kruskal.test(rel_abund ~ diet,
                       data = .))) -> kruskal_test

abun_filt %>% 
  na.omit() %>% 
  group_by(Family, day_post_inf) %>% 
  dunn_test(rel_abund ~ diet,
            p.adjust.method = 'BH',
            data = .) -> abun_dunn_test

## prepping for and putting together the statistical visualization based on dunns post hoc test
new_dunn_test <- stat_plot_prep(abun_filt,
                                abun_dunn_test)

abun_stat_vis <- stat_plot(new_dunn_test)

## saving my outputs as a .tsv
write_tsv(family_abun_lm,
          args$lm_FP)
write_tsv(abun_dunn_test,
          args$dunn_FP)

## saving statistical visualization
ggsave(args$stat_plot_FP,
       plot = abun_stat_vis, 
       width = 18, 
       height = 8)
