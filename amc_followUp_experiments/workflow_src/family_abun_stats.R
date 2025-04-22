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
# otu_table_FP <- './new_experiments/data/qiime/total_sum_otu_table.qza'
# tax_FP <- './new_experiments/data/qiime/taxonomy.qza'
# metadata_FP <- './new_experiments/data/misc/proc_newExp_d15-d3_metadata.tsv'
wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Morganellaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Tannerellaceae', 'Bacteroidaceae', 'Ruminococcaceae', 'Peptostreptococcaceae')

## function 1
family_abun_file_prep <- function(metadata_fp,
                                  tax_fp,
                                  otu_table_fp,
                                  tax_level,
                                  wanted_tax){
  ## metadata
  metadata <- read_tsv(metadata_fp)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'                                        
  ## taxonomy
  taxonomy <- read_qza(tax_fp)$data %>% 
    parse_taxonomy() %>% 
    rownames_to_column('asv')
  ## otu table 
  otu_table <- read_qza(otu_table_fp)$data
  otu_table %>% 
    as_tibble(rownames = 'asv') %>% 
    gather(-asv, key = sampleid, value = abun) %>% 
    group_by(sampleid) %>% 
    mutate(rel_abun = abun/sum(abun)) %>% 
    mutate(rel_abun = rel_abun + 0.000001) -> otu_table
  ## joining all tables together 
  otu_table %>% 
    left_join(metadata, by = 'sampleid') %>% 
    left_join(taxonomy, by = 'asv') -> abun_table
  abun_table %>% 
    group_by(sampleid, day_post_inf, diet, mouse_id, 
             purified_diet, high_fat, high_fiber, 
             .data[[tax_level]]) %>% 
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
## runs statistics on the large assembled data table for stat visualization
abun_stats <- function(filtered_table,
                       tax_level){
  ## linear modeling 
  filtered_table %>%
    na.omit() %>% 
    group_by(.data[[tax_level]], day_post_inf) %>% 
    do(glance(lm(rel_abund ~ purified_diet * high_fat * high_fiber,
                 data =.))) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(adj.p = p.adjust(p.value, 
                            method = "BH"),
           test_id = paste(.data[[tax_level]], day_post_inf, sep = "_")) -> lm_full
    #filter(adj.p <= 0.05) -> lm_full
  filtered_table %>%
    na.omit() %>% 
    group_by(.data[[tax_level]], day_post_inf) %>% 
    mutate(test_id = paste(.data[[tax_level]], day_post_inf, sep = "_")) %>% 
    filter(test_id %in% lm_full$test_id) %>% 
    do(tidy(lm(rel_abund ~ purified_diet * high_fat * high_fiber,
               data =.))) %>% 
    na.omit() %>% 
    filter(term != '(Intercept)') -> linear_model
  linear_model['signif'] <- symnum(linear_model$p.value,
                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                   symbols = c("****", "***", "**", "*", "ns"),
                                   abbr.colnames = FALSE,
                                   na = "")
  ## kruskal wallis test 
  filtered_table %>% 
    na.omit() %>% 
    group_by(.data[[tax_level]], day_post_inf) %>% 
    do(tidy(kruskal.test(rel_abund ~ diet,
                         data = .))) %>% 
    ungroup() %>% 
    arrange(p.value) %>% 
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(.data[[tax_level]], day_post_inf, sep = "_")) -> kruskal
    #filter(p.adj <= 0.05) -> kruskal
  ## dunn's post hoc test
  filtered_table %>% 
    na.omit() %>% 
    group_by(.data[[tax_level]], day_post_inf) %>%
    mutate(test_id = paste(.data[[tax_level]], day_post_inf, sep = "_")) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(rel_abund ~ diet,
              p.adjust.method = 'BH',
              data = .) -> dunn
  ## editing my dunn's post hoc test to include the difference in means between groups 
  filtered_table %>% 
    group_by(diet, .data[[tax_level]], day_post_inf) %>% 
    summarise(mean_rel_abund = mean(rel_abund)) -> mean_abun
  dunn %>% 
    merge(mean_abun, 
          by.x = c('group1',
                   'day_post_inf',
                   'Family'),
          by.y = c('diet',
                   'day_post_inf',
                   'Family')) %>% 
    rename('group1_rel_abun' = 'mean_rel_abund') %>% 
    merge(mean_abun,
          by.x = c('group2',
                   'day_post_inf',
                   'Family'),
          by.y = c('diet',
                   'day_post_inf',
                   'Family')) %>% 
    rename('group2_rel_abun' = 'mean_rel_abund') %>% 
    mutate(diff_means = (group1_rel_abun - group2_rel_abun),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  ## creating a list of my outputs
  my_list <- list(LinearModel = linear_model,
                  KruskalTest = kruskal,
                  DunnPostHoc = new_dunn)
  return(my_list)
}


## 3 
## statistical visualization 
stat_plot <- function(new_dunn,
                      tax_level){
  new_dunn %>% 
    ggplot(aes(x = group1, y = group2)) +
    geom_tile(aes(fill = stat_diff_means), alpha = 0.8, color = 'black') +
    scale_fill_gradient2(low = 'blue', high = 'green', name = 'Group 1 -\nGroup 2') +
    geom_text(aes(label = p.adj.signif)) +
    scale_x_discrete(labels = c('Chow',
                                'HFt/\nHFb',
                                'HFt/\nLFb',
                                'LFt/\nHFb')) +
    scale_y_discrete(labels = c('HFt / HFb',
                                'HFt / LFb',
                                'LFt / HFb',
                                'LFt / LFb')) +
    facet_grid(.data[[tax_level]]~day_post_inf,
               scales = 'free_x') +
    theme_bw(base_size = 20) +
    theme(strip.text.y = element_text(angle = 0)) +
    xlab('Group 1') +
    ylab('Group 2') -> stat_plot
  return(stat_plot)
}

## prepping the file needed to run the linear model
abun_files <- family_abun_file_prep(args$metadata_FP,
                                    args$tax_FP,
                                    args$otu_table_FP,
                                    wanted_level,
                                    wanted_family)

## can still pull out the metadata, taxonomy, and otu table as well
abun_filt <- abun_files$AbundanceTable

## performing statistical analysis
abun_stats <- abun_stats(abun_filt,
                         wanted_level)

kruskal_test <- abun_stats$KruskalTest
new_dunn_test <- abun_stats$DunnPostHoc
family_abun_lm <- abun_stats$LinearModel


## putting together the statistical visualization based on dunns post hoc test
abun_stat_vis <- stat_plot(new_dunn_test,
                           wanted_level)

## saving my outputs as a .tsv
write_tsv(family_abun_lm,
          args$lm_FP)
write_tsv(new_dunn_test,
          args$dunn_FP)

## saving statistical visualization
ggsave(args$stat_plot_FP,
       plot = abun_stat_vis, 
       width = 12, 
       height = 8)
