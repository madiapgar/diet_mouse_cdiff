## 6-26-23
## resiliency calculations for unweighted and weighted UniFrac distance matrices
## outputs a plot and the statistical results

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(rstatix)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-uu",
                    "--uu_dist",
                    dest = "uu_dist_fp",
                    help = "Filepath to Unweighted UniFrac Distance Matrix in .tsv format.")
parser$add_argument("-wu",
                    "--wu_dist",
                    dest = "wu_dist_fp",
                    help = "Filepath to Weighted UniFrac Distance Matrix in .tsv format")
parser$add_argument("-wlm",
                    "--weighted_lm",
                    dest = "wu_lm_fp",
                    help = "Filepath to Weighted Distance Matrix linear model results in .tsv format.")
parser$add_argument("-wd",
                    "--weighted_dunn",
                    dest = "wu_dunn_fp",
                    help = "Filepath to Weighted UniFrac Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-ulm",
                    "--unweighted_lm",
                    dest = "uu_lm_fp",
                    help = "Filepath to Unweighted Distance Matrix linear model results in .tsv format.")
parser$add_argument("-ud",
                    "--unweighted_dunn",
                    dest = "uu_dunn_fp",
                    help = "Filepath to Unweighted UniFrac Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-wp",
                    "--weighted_plot",
                    dest = "wu_plot_fp",
                    help = "Filepath to Weighted UniFrac resiliency plot in .pdf format.")
parser$add_argument("-wsp",
                    "--weighted_stat_plot",
                    dest = "wu_stat_plot_fp",
                    help = "Filepath to Weighted UniFrac resiliency statistical plot in .pdf format.")
parser$add_argument("-up",
                    "--unweighted_plot",
                    dest = "uu_plot_fp",
                    help = "Filepath to Unweighted UniFrac resiliency plot in .pdf format.")
parser$add_argument("-usp",
                    "--unweighted_stat_plot",
                    dest = "uu_stat_plot_fp",
                    help = "Filepath to Unweighted UniFrac resiliency statistical plot in .pdf format.")

args <- parser$parse_args()

## input file paths
# metadata_FP <- './data/misc/processed_metadata.tsv'
# uu_dist_fp <- './data/qiime/core_outputs/uw_dist_matrix.tsv'
# wu_dist_fp <- './data/qiime/core_outputs/w_dist_matrix.tsv'

## lists to redo the diet names on the facet labels of the ggplot created below 
diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')

names(diet_labs) <- c('Chow', 'HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')

## functions in order of usage
## 1 
## this function combines your metadata file and your distance matrix 
## it also pulls out the wanted diet variable from the metadata file 
## requires the input of a distance matrix (dist_mat), metadata file (meta),
## and the list of unique diets (mydiet)
## used in subsequent resil/homog_diet_assembly functions
dist_filter <- function(dist_mat, meta, mydiet){
  meta %>% 
    filter(diet == mydiet) -> min_meta
  
  dist_mat <- data.frame(dist_mat)
  rownames(dist_mat) <- dist_mat[,1]
  dist_mat <- dist_mat[,-1]
  
  sample_idx <- rownames(dist_mat) %in% min_meta$sampleid 
  
  min_dist_mat <- dist_mat[sample_idx, sample_idx]
  min_dist_mat <- as_tibble(min_dist_mat, rownames = 'sampleid')
  colnames(min_dist_mat)[2:ncol(min_dist_mat)] <- min_dist_mat$sampleid
  return(min_dist_mat)
}

## 2
## this function takes the output from dist_filter and the metadata file 
## and pulls all unique distance matrix comparisons out for that particular diet 
## and places them in a new column 
## it also esta the baseline that all sample pairs will be compared to 
## requires the distance matrix output from the dist_filter (min_dist_mat),
## and the metadata file (meta)
resil_dist_comp <- function(min_dist_mat, meta){
  meta %>% 
    select(sampleid, diet, day_post_inf, high_fat, high_fiber, purified_diet,
           seq_depth, mouse_id) -> min_meta
  
  min_dist_mat %>% 
    rename(sampleid_x = sampleid) %>% 
    gather(-sampleid_x, key = sampleid_y, value = dist) %>% 
    merge(min_meta, by.x = 'sampleid_x', by.y = 'sampleid') %>% 
    merge(min_meta, by.x = 'sampleid_y', by.y = 'sampleid') -> long_dist_mat
  
  long_dist_mat %>% 
    mutate(key = paste0(pmin(sampleid_x, sampleid_y), 
                        pmax(sampleid_x, sampleid_y), 
                        sep = "_")) %>% 
    distinct(key, .keep_all = TRUE) %>% 
    filter(sampleid_x != sampleid_y,
           ## for homogenicity testing, do day_post_inf.x == day_post_inf.y 
           ## and don't filter out day -15 
           day_post_inf.x == -8,
           day_post_inf.y > -8) %>% 
    select(-diet.y) %>% 
    rename(diet = diet.x) -> long_dist_mat
  
  return(long_dist_mat)
}

## 3
## this is a for loop that does the above analysis for each diet variable and places them 
## in a list of dataframes that can be bound together to create one plot for comparison
## requires the input of a distance matrix (dist_mat), metadata file (meta),
## and the list of unique diets (mydiet)
resil_diet_assembly <- function(dist_mat, 
                                meta,
                                mydiet){
  resil <- list()
  for (i in 1:length(mydiet)){
    tmp_dist_mat <- dist_filter(dist_mat, meta, mydiet[i])
    tmp_resil <- resil_dist_comp(tmp_dist_mat, meta)
    resil[[i]] <- tmp_resil
  }
  resil <- bind_rows(resil)
  return(resil)
}

## 4
## all functions for calculating statistical analysis and creating a plot visualization
stats <- function(biom_table){
  ## linear modeling 
  biom_table %>% 
    filter(day_post_inf.y != -15) %>% 
    group_by(day_post_inf.y) %>% 
    do(glance(lm(dist ~ (purified_diet.y * seq_depth.y) + high_fat.y * high_fiber.y,
                 data = .))) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(adj.p = p.adjust(p.value, 
                            method = "BH"),
           test_id = paste(day_post_inf.y)) %>%
    filter(p.value <= 0.05) -> lm_full
  
  biom_table %>% 
    group_by(day_post_inf.y) %>% 
    mutate(test_id = paste(day_post_inf.y)) %>% 
    filter(test_id %in% lm_full$test_id) %>% 
    do(tidy(lm(dist ~ (purified_diet.y * seq_depth.y) + high_fat.y * high_fiber.y,
               data = .))) %>%
    filter(term != '(Intercept)') %>% 
    na.omit() -> lm_results
  
  lm_results['signif'] <- symnum(lm_results$p.value,
                                 cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                                 symbols = c("****", "***", "**", "*", "+", "ns"),
                                 abbr.colnames = FALSE,
                                 na = "")
  ## kruskal wallis and dunns post hoc test
  biom_table %>% 
    na.omit() %>% 
    filter(day_post_inf.y != -15) %>% 
    group_by(day_post_inf.y) %>% 
    do(tidy(kruskal.test(dist ~ diet,
                         data = .))) %>% 
    ungroup() %>% 
    arrange(p.value) %>% 
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(day_post_inf.y)) %>%
    filter(p.adj <= 0.05) -> kruskal
  
  biom_table %>% 
    na.omit() %>% 
    group_by(day_post_inf.y) %>% 
    mutate(test_id = paste(day_post_inf.y)) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(dist ~ diet,
              p.adjust.method = 'BH',
              data = .) -> dunn
  ## creating a named list of results
  my_list <- list(LinearModel = lm_results,
                  KruskalTest = kruskal,
                  DunnPostHoc = dunn)
  return(my_list)
}

## 5
## prepping the dunn test for the statistical visualization
stat_plot_prep <- function(biom_table,
                           dunn_test){
  biom_table %>% 
    group_by(diet, day_post_inf.y) %>% 
    summarise(mean_dist = mean(dist)) -> mean_dist
  dunn_test %>% 
    merge(mean_dist, 
          by.x = c('group1',
                   'day_post_inf.y'),
          by.y = c('diet',
                   'day_post_inf.y')) %>% 
    rename('group1_dist' = 'mean_dist') %>% 
    merge(mean_dist,
          by.x = c('group2',
                   'day_post_inf.y'),
          by.y = c('diet',
                   'day_post_inf.y')) %>% 
    rename('group2_dist' = 'mean_dist') %>% 
    mutate(diff_means = (group1_dist - group2_dist),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  return(new_dunn)
}

## 5
## actual statistical visualization 
stat_plot <- function(new_dunn){
  new_dunn %>% 
    filter(day_post_inf.y != -15) %>%
    ggplot(aes(x = group1, y = group2)) +
    geom_tile(aes(fill = stat_diff_means), alpha = 0.8, color = 'black') +
    scale_fill_gradient2(low = 'green', high = 'blue', name = 'Group 1 -\nGroup 2',
                         trans = 'reverse') +
    geom_text(aes(label = p.adj.signif)) +
    scale_x_discrete(labels = c('Chow',
                                'HFt/\nHFb',
                                'HFt/\nLFb',
                                'LFt/\nHFb')) +
    scale_y_discrete(labels = c('HFt / HFb',
                                'HFt / LFb',
                                'LFt / HFb',
                                'LFt / LFb')) +
    facet_grid(~day_post_inf.y) +
    theme_bw(base_size = 16) +
    theme(strip.text.y = element_text(angle = 0)) +
    xlab('Group 1') +
    ylab('Group 2') -> stat_vis
  return(stat_vis)
}

## file prep for resil. plot and stats
meta <- read_tsv(args$metadata_FP)
names(meta)[names(meta) == '#SampleID'] <- 'sampleid'

uu_dist <- read_tsv(args$uu_dist_fp)
wu_dist <- read_tsv(args$wu_dist_fp)

## pulling out unique diets 
diets <- unique(meta$diet)

## resil calculations 
## unweighted unifrac
uu_resil <- resil_diet_assembly(uu_dist,
                                meta,
                                diets)

## weighted unifrac
wu_resil <- resil_diet_assembly(wu_dist,
                                meta,
                                diets)

## unweighted UniFrac plot and stats 
## plot 
uu_resil %>% 
  ggplot(aes(x = day_post_inf.y, y = dist)) +
  geom_boxplot((aes(group = day_post_inf.y)), outlier.shape = NA) +
  geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
  geom_line(aes(group = mouse_id.y), alpha = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_smooth(se = FALSE, method = 'loess') + 
  scale_x_continuous(breaks = c(-3, 0, 3)) + 
  scale_y_reverse() +
  facet_wrap(~diet, 
             labeller = labeller(diet = diet_labs),
             nrow = 1) +
  theme_bw(base_size = 16) +
  xlab('Days Relative to Infection') +
  ylab('Unweighted UniFrac Distance\n(Pairwise Distances from Day -8)') +
  ggtitle("Microbiome Resilience Over Time") -> uu_resil_plot

## stats
uu_resil_stats <- stats(uu_resil)

uu_resil_lm <- uu_resil_stats$LinearModel
uu_resil_kruskal <- uu_resil_stats$KruskalTest
uu_resil_dunn <- uu_resil_stats$DunnPostHoc

## statistical visualization
stat_plot_prep(uu_resil,
               uu_resil_dunn) -> new_uu_resil_dunn

stat_plot(new_uu_resil_dunn) -> uu_resil_stat_vis


## weighted UniFrac plot and stats
## plot
wu_resil %>% 
  ggplot(aes(x = day_post_inf.y, y = dist)) +
  geom_boxplot((aes(group = day_post_inf.y)), outlier.shape = NA) +
  geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
  geom_line(aes(group = mouse_id.y), alpha = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_smooth(se = FALSE, method = 'loess') + 
  scale_x_continuous(breaks = c(-3, 0, 3)) + 
  scale_y_reverse() +
  facet_wrap(~diet, 
             labeller = labeller(diet = diet_labs),
             nrow = 1) +
  theme_bw(base_size = 16) +
  xlab('Days Relative to Infection') +
  ylab('Weighted UniFrac Distance\n(Pairwise Distances from Day -8)') +
  ggtitle("Microbiome Resilience Over Time") -> wu_resil_plot

## stats
wu_resil_stats <- stats(wu_resil)

wu_resil_lm <- wu_resil_stats$LinearModel
wu_resil_kruskal <- wu_resil_stats$KruskalTest
wu_resil_dunn <- wu_resil_stats$DunnPostHoc

## statistical visualization
stat_plot_prep(wu_resil,
               wu_resil_dunn) -> new_wu_resil_dunn

stat_plot(new_wu_resil_dunn) -> wu_resil_stat_vis

## saving my output plots and stats
## plots
ggsave(args$wu_plot_fp,
       plot = wu_resil_plot,
       width = 12,
       height = 4)
ggsave(args$wu_stat_plot_fp,
       plot = wu_resil_stat_vis,
       width = 14,
       height = 4)
ggsave(args$uu_plot_fp,
       plot = uu_resil_plot,
       width = 12,
       height = 4)
ggsave(args$uu_stat_plot_fp,
       plot = uu_resil_stat_vis,
       width = 14,
       height = 4)

## stats
write_tsv(uu_resil_lm,
          args$uu_lm_fp)
write_tsv(new_uu_resil_dunn,
          args$uu_dunn_fp)
write_tsv(wu_resil_lm,
          args$wu_lm_fp)
write_tsv(new_wu_resil_dunn,
          args$wu_dunn_fp)
