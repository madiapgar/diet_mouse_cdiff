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
parser$add_argument("-flms",
                    "--faith_lm_sec",
                    dest = "faith_lm_sec_FP",
                    help = "Filepath to Faith's PD sectioned linear model results in .tsv format.")
parser$add_argument("-fd",
                    "--faith_dunn",
                    dest = "faith_dunn_FP",
                    help = "Filepath to Faith's PD Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-slm",
                    "--shannon_lm",
                    dest = "shannon_lm_FP",
                    help = "Filepath to Shannon Entropy total linear model results in .tsv format.")
parser$add_argument("-slms",
                    "--shannon_lm_sec",
                    dest = "shannon_lm_sec_FP",
                    help = "Filepath to Shannon Entropy sectioned linear model results in .tsv format.")
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
unwanted_samples <- c('Mock20220615A', 'Mock_1A', 'Mock_2A',
                      'Mock_3A', 'Mock_4A', 'Mock_5A', 'Mock_6A',
                      'Mock_7A', 'PCR Blank0',
                      'PCR Blank1', 'Mock_7', 'Mock_6',
                      'Mock_5', 'Mock_4', 'PCR blank')

## functions in order that they're used
## 1
## alpha diversity file prep 
alpha_div_prep <- function(file_path1,
                           file_path2,
                           sample_filter,
                           metadata_fp){
  ## faith's pd 
  alpha_faith <- read_tsv(file_path1)
  names(alpha_faith)[names(alpha_faith) == '#SampleID'] <- 'sampleid'
  alpha_faith %>% 
    filter(!(sampleid %in% sample_filter)) -> faith_pd
  ## metadata file for both
  stat_meta <- read_tsv(metadata_fp)
  stat_meta %>% 
    filter(!(sampleid %in% sample_filter)) -> stat_meta
  ## joining faith's pd and metadata file together into one table
  stat_meta %>% 
    filter(sampleid %in% faith_pd$sampleid) %>% 
    left_join(faith_pd, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> faith_biom
  ## shannon entropy
  alpha_shannon <- read_tsv(file_path2)
  names(alpha_shannon)[names(alpha_shannon) == '...1'] <- 'sampleid'
  alpha_shannon %>% 
    filter(!(sampleid %in% sample_filter)) -> shannon
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
  ## alpha_cat is what the alpha div column is called (faith_pd or shannon_entropy)
  ## sectioned out by diet 
  biom_table %>% 
    group_by(day_post_inf) %>% 
    do(tidy(lm(faith_pd ~ (purified_diet * seq_depth) + high_fat + high_fiber + study,
               data = .))) -> sectioned_lm
  sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> sectioned_lm
  ## not sectioned out by diet 
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(faith_pd ~ diet * seq_depth,
               data = .))) -> not_sectioned_lm
  not_sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> not_sectioned_lm
  ## kruskal wallis and dunns post hoc tests
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    do(tidy(kruskal.test(faith_pd ~ diet,
                         data = .))) -> kruskal
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
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
    do(tidy(lm(shannon_entropy ~ (purified_diet * seq_depth) + high_fat + high_fiber + study,
               data = .))) -> sectioned_lm
  sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> sectioned_lm
  ## not sectioned out by diet 
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(shannon_entropy ~ diet * seq_depth,
               data = .))) -> not_sectioned_lm
  not_sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> not_sectioned_lm
  ## kruskal wallis and dunns post hoc tests
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
    do(tidy(kruskal.test(shannon_entropy ~ diet,
                         data = .))) -> kruskal
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf) %>% 
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
                              unwanted_samples,
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
write_tsv(faith_lm, args$faith_lm_FP)
write_tsv(sectioned_faith_lm, args$faith_lm_sec_FP)
write_tsv(new_faith_dunn, args$faith_dunn_FP)
write_tsv(shannon_lm, args$shannon_lm_FP)
write_tsv(sectioned_shannon_lm, args$shannon_lm_sec_FP)
write_tsv(new_shannon_dunn, args$shannon_dunn_FP)

## saving my statistical visualizations
ggsave(args$faith_plot_FP,
       plot = faith_stat_vis, 
       width = 13, 
       height = 3)

ggsave(args$shannon_plot_FP,
       plot = shannon_stat_vis, 
       width = 13, 
       height = 3)
