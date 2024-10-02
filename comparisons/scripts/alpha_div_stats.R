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
library(glue)

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
                    help = "Filepath to Faith's PD sectioned linear model results in .tsv format.")
parser$add_argument("-fd",
                    "--faith_dunn",
                    dest = "faith_dunn_FP",
                    help = "Filepath to Faith's PD Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-slm",
                    "--shannon_lm",
                    dest = "shannon_lm_FP",
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

## file paths
# metadata_FP <- './comparisons/data/misc/newExp_comp_d15_metadata.tsv'
# faith_pd_FP <- './comparisons/data/qiime/core_outputs/faith_pd.tsv'
# shannon_FP <- './comparisons/data/qiime/core_outputs/shannon_entropy.tsv'


## functions in order that they're used
## 1
## alpha diversity file prep 
## alpha diversity file prep 
alpha_div_prep <- function(faith_fp,
                           shannon_fp,
                           metadata_fp){
  ## faith's pd 
  faith_pd <- read_tsv(faith_fp)
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
  shannon <- read_tsv(shannon_fp)
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
## linear model calculations
linear_model <- function(input_table,
                         grouped_by,
                         adjust_method,
                         filter_adj_p_value = FALSE,
                         formula_left,
                         formula_right){
  ## linear modeling
  funky_formula <- paste(formula_left, formula_right, sep = "~")
  pre_lm <- input_table %>%
    na.omit() %>%
    group_by(input_table[grouped_by]) %>%
    do(glance(lm(as.formula(funky_formula),
                 data = .))) %>%
    ungroup() %>%
    na.omit() %>%
    mutate(p.adj = p.adjust(p.value,
                            method = adjust_method))
  
  if (length(grouped_by) > 1) {
    mini_pre_lm <- pre_lm %>%
      select(grouped_by)
    
    pre_lm <- cbind(pre_lm, test_id=do.call(paste, c(mini_pre_lm, sep = "_")))
    
    mini_input <- input_table %>%
      select(grouped_by)
    
    alt_input <- cbind(input_table, test_id=do.call(paste, c(mini_input, sep = "_")))
    
  } else {
    pre_lm <- pre_lm %>%
      mutate(test_id = paste(.data[[grouped_by[1]]]))
    
    alt_input <- input_table %>%
      mutate(test_id = paste(.data[[grouped_by[1]]]))
  }
  
  if (filter_adj_p_value == TRUE) {
    pre_lm <- pre_lm %>%
      filter(p.adj <= 0.05)
    
    linear_model_results <- alt_input %>%
      na.omit() %>%
      group_by(input_table[grouped_by]) %>%
      filter(test_id %in% pre_lm$test_id) %>%
      do(tidy(lm(as.formula(funky_formula),
                 data = .))) %>%
      na.omit() %>%
      filter(term != '(Intercept)')
  } else {
    pre_lm
    
    linear_model_results <- alt_input %>%
      na.omit() %>%
      group_by(input_table[grouped_by]) %>%
      filter(test_id %in% pre_lm$test_id) %>%
      do(tidy(lm(as.formula(funky_formula),
                 data = .))) %>%
      na.omit() %>%
      filter(term != '(Intercept)')
  }
  
  ## assigning p-value significance
  linear_model_results['signif'] <- symnum(linear_model_results$p.value,
                                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                           symbols = c("****", "***", "**", "*", "ns"),
                                           abbr.colnames = FALSE,
                                           na = "")
  return(linear_model_results)
}


## 3 
## kruskal wallis and dunns post hoc test calculations
kruskal_dunn_stats <- function(input_table,
                               grouped_by,
                               adjust_method,
                               filter_adj_p_value = FALSE,
                               formula_left,
                               formula_right){
  # kruskal test
  kruskal_results <- input_table %>%
    group_by(input_table[grouped_by]) %>%
    do(tidy(kruskal.test(.data[[formula_left]] ~ .data[[formula_right]],
                          data = .))) %>%
    ungroup() %>%
    arrange(p.value) %>%
    mutate(p.adj = p.adjust(p.value,
                            method = adjust_method))
  
  
  if (length(grouped_by) > 1) {
    mini_kruskal_results <- kruskal_results %>%
      select(grouped_by)
    
    kruskal_results <- cbind(kruskal_results, test_id=do.call(paste, c(mini_kruskal_results,
                                                                       sep = "_")))
    
    mini_input <- input_table %>%
      select(grouped_by)
    
    alt_input <- cbind(input_table, test_id=do.call(paste, c(mini_input, sep = "_")))
    
  } else {
    kruskal_results <- kruskal_results %>%
      mutate(test_id = paste(.data[[grouped_by[1]]]))
    
    alt_input <- input_table %>%
      mutate(test_id = paste(.data[[grouped_by[1]]]))
  }
  
  ## dunns post hoc test
  ## need to use reformulate/glue to have function variables work with the dunn test formula
  rightSide_name <- formula_right
  funky_formula <- reformulate(glue("{rightSide_name}"),
                               glue("{formula_left}"))
  
  if (filter_adj_p_value == TRUE) {
    kruskal_results <- kruskal_results %>%
      filter(p.adj <= 0.05)
    
    dunn_results <- alt_input %>%
      group_by(alt_input[grouped_by]) %>%
      filter(test_id %in% kruskal_results$test_id) %>%
      dunn_test(funky_formula,
                p.adjust.method = adjust_method,
                data = .) %>%
      add_xy_position(scales = 'free',
                      fun = 'max')
  } else {
    kruskal_results
    
    dunn_results <- alt_input %>%
      group_by(alt_input[grouped_by]) %>%
      dunn_test(funky_formula,
                p.adjust.method = adjust_method,
                data = .) %>%
      add_xy_position(scales = 'free',
                               fun = 'max')
  }
  
  ## list of outputs
  my_list <- list(KruskalTest = kruskal_results,
                  DunnTest = dunn_results)
  return(my_list)
}


## 4 
## preps dunns post hoc results for statistical visualization
stat_plot_prep <- function(filtered_table,
                           first_group,
                           mean_value,
                           dunn_test){
  filtered_table %>% 
    group_by(.data[[first_group]]) %>% 
    summarise(mean = mean(.data[[mean_value]])) -> mean_table
  
  dunn_test %>% 
    merge(mean_table, 
          by.x = c('group1'),
          by.y = c(first_group)) %>% 
    rename_with(~paste0('group1_', mean_value, recycle0 = TRUE), contains('mean')) %>% 
    merge(mean_table,
          by.x = c('group2'),
          by.y = c(first_group)) %>% 
    rename_with(~paste0('group2_', mean_value, recycle0 = TRUE), contains('mean')) -> int_dunn
  
  group1_col <- paste0('group1_', mean_value)
  group2_col <- paste0('group2_', mean_value)
  
  int_dunn %>% 
    mutate(diff_means = (.data[[group1_col]] - .data[[group2_col]]),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  
  return(new_dunn)
}


## 5 
## statistical visualization 
stat_plot <- function(new_dunn,
                      title){
  new_dunn %>% 
    ggplot(aes(x = group1, y = group2)) +
    geom_tile(aes(fill = stat_diff_means), alpha = 0.8, color = 'black') +
    scale_fill_gradient2(low = 'blue', high = 'green', name = 'Group 1 -\nGroup 2') +
    geom_text(aes(label = p.adj.signif)) +
    theme_bw(base_size = 20) +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_discrete(labels = c('New Anschutz (2024)',
                                'U of Arizona')) +
    scale_x_discrete(labels = c('Old Anschutz (2020)',
                                'New Anschutz (2024)')) +
    xlab('Group 1') +
    ylab('Group 2') +
    ggtitle(title) -> stat_vis
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
faith_lm <- linear_model(input_table = faith,
                         grouped_by = 'vendor',
                         adjust_method = 'BH',
                         filter_adj_p_value = FALSE,
                         formula_left = 'faith_pd',
                         formula_right = 'experiment_set')

faith_stats <- kruskal_dunn_stats(input_table = faith,
                                  grouped_by = 'vendor',
                                  adjust_method = 'BH',
                                  filter_adj_p_value = FALSE,
                                  formula_left = 'faith_pd',
                                  formula_right = 'experiment_set')

faith_kruskal <- faith_stats$KruskalTest
faith_dunn <- faith_stats$DunnTest

new_faith_dunn <- stat_plot_prep(filtered_table = faith,
                                 first_group = 'experiment_set',
                                 mean_value = 'faith_pd',
                                 dunn_test = faith_dunn)

faith_stat_vis <- stat_plot(new_faith_dunn,
                            "All Exp Faith's PD")

## shannon entropy stats and visualization 
shannon_lm <- linear_model(input_table = shannon,
                           grouped_by = 'vendor',
                           adjust_method = 'BH',
                           filter_adj_p_value = FALSE,
                           formula_left = 'shannon_entropy',
                           formula_right = 'experiment_set')

shannon_stats <- kruskal_dunn_stats(input_table = shannon,
                                    grouped_by = 'vendor',
                                    adjust_method = 'BH',
                                    filter_adj_p_value = FALSE,
                                    formula_left = 'shannon_entropy',
                                    formula_right = 'experiment_set')

shannon_kruskal <- shannon_stats$KruskalTest
shannon_dunn <- shannon_stats$DunnTest

new_shannon_dunn <- stat_plot_prep(filtered_table = shannon,
                                   first_group = 'experiment_set',
                                   mean_value = 'shannon_entropy',
                                   dunn_test = shannon_dunn)

shannon_stat_vis <- stat_plot(new_shannon_dunn,
                              "All Exp Shannon Entropy")

## writing out results as a .tsv file 
write_tsv(faith_lm, args$faith_lm_FP)
write_tsv(new_faith_dunn, args$faith_dunn_FP)
write_tsv(shannon_lm, args$shannon_lm_FP)
write_tsv(new_shannon_dunn, args$shannon_dunn_FP)

## saving my statistical visualizations
ggsave(args$faith_plot_FP,
       plot = faith_stat_vis, 
       width = 9, 
       height = 4)

ggsave(args$shannon_plot_FP,
       plot = shannon_stat_vis, 
       width = 9, 
       height = 4)
