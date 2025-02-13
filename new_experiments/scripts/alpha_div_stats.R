## 6-26-23
## Qiime2 core diversity analysis statistical analysis for alpha diversity metrics

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(ape)
library(rstatix)
library(ggh4x)
library(vegan)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-f",
                    "--faith_pd",
                    dest = "faith_pd_FP",
                    help = "Filepath to Faith's PD file in .tsv format.")
parser$add_argument("-s",
                    "--shannon",
                    dest = "shannon_FP",
                    help = "Filepath to Shannon Entropy file in .tsv format.")
parser$add_argument("-flm",
                    "--faith_lm",
                    dest = "faith_lm_FP",
                    help = "Filepath to Faith's PD total linear model results in .tsv format.")
parser$add_argument("-fd",
                    "--faith_dunn",
                    dest = "faith_dunn_FP",
                    help = "Filepath to Faith's PD Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-slm",
                    "--shannon_lm",
                    dest = "shannon_lm_FP",
                    help = "Filepath to Shannon Entropy total linear model results in .tsv format.")
parser$add_argument("-sd",
                    "--shannon_dunn",
                    dest = "shannon_dunn_FP",
                    help = "Filepath to Shannon Entropy Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-fp",
                    "--faith_plot",
                    dest = "faith_plot_FP",
                    help = "Filepath to Faith's PD plot in .pdf format.")
parser$add_argument("-sp",
                    "--shannon_plot",
                    dest = "shannon_plot_FP",
                    help = "Filepath to Shannon Entropy plot in .pdf format.")

args <- parser$parse_args()

## input file paths
# metadata_FP <- './data/misc/processed_metadata.tsv'
# faith_pd_FP <- './data/qiime/core_outputs/faith_pd.tsv'
# shannon_FP <- './data/qiime/core_outputs/shannon_entropy.tsv'


## functions in order that they're used
## 1
## alpha diversity file prep 
## alpha diversity file prep 
alpha_div_prep <- function(file_path1,
                           file_path2,
                           metadata_fp){
  ## faith's pd 
  faith_pd <- read_tsv(file_path1)
  names(faith_pd)[names(faith_pd) == '#SampleID'] <- 'sampleid'
  ## metadata file for both
  stat_meta <- read_tsv(metadata_fp)
  names(stat_meta)[names(stat_meta) == '#SampleID'] <- 'sampleid'
  ## joining faith's pd and metadata file together into one table
  stat_meta %>% 
    filter(sampleid %in% faith_pd$sampleid) %>% 
    left_join(faith_pd, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> faith_biom
  ## shannon entropy
  shannon <- read_tsv(file_path2)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  ## joining shannon and metadata file together into one table 
  stat_meta %>% 
    filter(sampleid %in% shannon$sampleid) %>% 
    left_join(shannon, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> shannon_biom
  ## creating a list for outputs 
  my_list <- list(FaithPD = faith_biom,
                  Shannon = shannon_biom, 
                  Metadata = stat_meta)
  return(my_list)
}

## 2
## stats calculations
## faith's pd 
faith_div_stats <- function(biom_table){
  ## sectioned out by diet 
  biom_table %>% 
    group_by(day_post_inf) %>% 
    do(glance(lm(faith_pd ~ purified_diet * high_fat * high_fiber,
                 data = .))) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(adj.p = p.adjust(p.value, 
                            method = "BH"),
           test_id = paste(day_post_inf)) -> lm_full
    #filter(adj.p <= 0.05) -> lm_full
  biom_table %>% 
    group_by(day_post_inf) %>% 
    mutate(test_id = paste(day_post_inf)) %>% 
    filter(test_id %in% lm_full$test_id) %>% 
    do(tidy(lm(faith_pd ~ purified_diet * high_fat * high_fiber,
               data = .))) %>%
    filter(term != '(Intercept)') -> sectioned_lm
  sectioned_lm['signif'] <- symnum(sectioned_lm$p.value,
                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                                   symbols = c("****", "***", "**", "*", "+", "ns"),
                                   abbr.colnames = FALSE,
                                   na = "")
  ## not sectioned out by diet 
  ## haven't used these results much so decided not to do anything to this
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(faith_pd ~ diet,
               data = .))) -> not_sectioned_lm
  #not_sectioned_lm %>% 
    #filter(p.value <= 0.05) -> not_sectioned_lm
  ## kruskal wallis and dunns post hoc tests
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    do(tidy(kruskal.test(faith_pd ~ diet,
                         data = .))) %>% 
    ungroup() %>% 
    arrange(p.value) %>% 
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(day_post_inf)) -> kruskal
    #filter(p.adj <= 0.05) -> kruskal
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    mutate(test_id = paste(day_post_inf)) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(faith_pd ~ diet,
              p.adjust.method = 'BH',
              data = .) -> dunn
  ## creating a list 
  my_list <- list(DietSpecific = sectioned_lm,
                  OverallDiet = not_sectioned_lm,
                  KruskalTest = kruskal,
                  DunnPostHoc = dunn)
  return(my_list)
}

## 3
## shannon entropy 
shannon_div_stats <- function(biom_table){
  ## alpha_cat is what the alpha div column is called (faith_pd or shannon_entropy)
  ## sectioned out by diet 
  biom_table %>% 
    group_by(day_post_inf) %>% 
    do(glance(lm(shannon_entropy ~ purified_diet * high_fat * high_fiber,
                 data = .))) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(adj.p = p.adjust(p.value, 
                            method = "BH"),
           test_id = paste(day_post_inf)) -> lm_full
    #filter(adj.p <= 0.05) -> lm_full
  biom_table %>% 
    group_by(day_post_inf) %>% 
    mutate(test_id = paste(day_post_inf)) %>% 
    filter(test_id %in% lm_full$test_id) %>% 
    do(tidy(lm(shannon_entropy ~ purified_diet * high_fat * high_fiber,
               data = .))) %>%
    filter(term != '(Intercept)') -> sectioned_lm
  sectioned_lm['signif'] <- symnum(sectioned_lm$p.value,
                                   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                                   symbols = c("****", "***", "**", "*", "+", "ns"),
                                   abbr.colnames = FALSE,
                                   na = "")
  ## not sectioned out by diet 
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(shannon_entropy ~ diet,
               data = .))) -> not_sectioned_lm
  #not_sectioned_lm %>% 
    #filter(p.value <= 0.05) -> not_sectioned_lm
  ## kruskal wallis and dunns post hoc tests
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    do(tidy(kruskal.test(shannon_entropy ~ diet,
                         data = .))) %>% 
    ungroup() %>% 
    arrange(p.value) %>% 
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(day_post_inf)) -> kruskal
    #filter(p.adj <= 0.05) -> kruskal
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    mutate(test_id = paste(day_post_inf)) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(shannon_entropy ~ diet,
              p.adjust.method = 'BH',
              data = .) -> dunn
  ## creating a list 
  my_list <- list(DietSpecific = sectioned_lm,
                  OverallDiet = not_sectioned_lm,
                  KruskalTest = kruskal,
                  DunnPostHoc = dunn)
  return(my_list)
}

## 4 
## preps dunns post hoc results for statistical visualization
stat_plot_prep <- function(biom_table,
                           dunn_test,
                           value){
  biom_table %>% 
    group_by(diet, day_post_inf) %>% 
    summarise(means = mean(.data[[value]])) -> mean_table
  dunn_test %>% 
    merge(mean_table, 
          by.x = c('group1',
                   'day_post_inf'),
          by.y = c('diet',
                   'day_post_inf')) %>% 
    rename('group1_means' = 'means') %>% 
    merge(mean_table,
          by.x = c('group2',
                   'day_post_inf'),
          by.y = c('diet',
                   'day_post_inf')) %>% 
    rename('group2_means' = 'means') %>% 
    mutate(diff_means = (group1_means - group2_means),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  return(new_dunn)
}

## 5 
## statistical visualization 
stat_plot <- function(new_dunn){
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
    facet_grid(~day_post_inf,
               scales = 'free_x') +
    theme_bw(base_size = 16) +
    theme(strip.text.y = element_text(angle = 0)) +
    xlab('Group 1') +
    ylab('Group 2') -> stat_vis
  return(stat_vis)
}

## use of functions 
## alpha diversity analysis  
alpha_files <- alpha_div_prep(args$faith_pd_FP,
                              args$shannon_FP,
                              args$metadata_FP)

faith <- alpha_files$FaithPD
shannon <- alpha_files$Shannon
metadata <- alpha_files$Metadata

## faith's pd stats and visualization
faith_stats <- faith_div_stats(faith)
sectioned_faith_lm <- faith_stats$DietSpecific
faith_lm <- faith_stats$OverallDiet
faith_kruskal <- faith_stats$KruskalTest
faith_dunn <- faith_stats$DunnPostHoc

new_faith_dunn <- stat_plot_prep(faith,
                                 faith_dunn,
                                 'faith_pd')

faith_stat_vis <- stat_plot(new_faith_dunn)

## shannon entropy stats and visualization 
shannon_stats <- shannon_div_stats(shannon)
sectioned_shannon_lm <- shannon_stats$DietSpecific
shannon_lm <- shannon_stats$OverallDiet
shannon_kruskal <- shannon_stats$KruskalTest
shannon_dunn <- shannon_stats$DunnPostHoc

new_shannon_dunn <- stat_plot_prep(shannon,
                                   shannon_dunn,
                                   'shannon_entropy')

shannon_stat_vis <- stat_plot(new_shannon_dunn)

## writing out results as a .tsv file 
write_tsv(sectioned_faith_lm, args$faith_lm_FP)
write_tsv(new_faith_dunn, args$faith_dunn_FP)
write_tsv(sectioned_shannon_lm, args$shannon_lm_FP)
write_tsv(new_shannon_dunn, args$shannon_dunn_FP)

## saving my statistical visualizations
ggsave(args$faith_plot_FP,
       plot = faith_stat_vis, 
       width = 8, 
       height = 3)

ggsave(args$shannon_plot_FP,
       plot = shannon_stat_vis, 
       width = 8, 
       height = 3)
