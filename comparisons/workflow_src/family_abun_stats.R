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
library(glue)

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

wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Bacteroidaceae', 'Tannerellaceae', 'Morganellaceae')
wanted_genus <- c('Akkermansia', 'Enterococcus', 'Escherichia-Shigella', 'Proteus', 
                  'Bacteroides', 'Lactobacillus', 'Staphylococcus', 'Muribaculaceae')

group1_labs <- c('New Anschutz (2024)')
names(group1_labs) <- c('new_exp_anschutz')

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
    group_by(sampleid, diet, mouse_id,
             experiment_set, vendor, mouse_sex,
             .data[[tax_level]]) %>% 
    summarise(rel_abund = sum(rel_abun)) %>% 
    filter(.data[[tax_level]] %in% wanted_tax) %>% 
    mutate(mouse_fact = as.factor(mouse_id)) -> abun_filt
  ## creating a list for my outputs
  my_list <- list(Metadata = metadata,
                  Taxonomy = taxonomy,
                  OTUTable = otu_table,
                  AbundanceTable = abun_filt)
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
    mutate(p.adj = p.adjust(.data$p.value,
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
      filter(.data$p.adj <= 0.05)
    
    linear_model_results <- alt_input %>%
      na.omit() %>%
      group_by(input_table[grouped_by]) %>%
      filter(.data$test_id %in% pre_lm$test_id) %>%
      do(tidy(lm(as.formula(funky_formula),
                 data = .))) %>%
      na.omit() %>%
      filter(term != '(Intercept)')
  } else {
    pre_lm
    
    linear_model_results <- alt_input %>%
      na.omit() %>%
      group_by(input_table[grouped_by]) %>%
      filter(.data$test_id %in% pre_lm$test_id) %>%
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
    arrange(.data$p.value) %>%
    mutate(p.adj = p.adjust(.data$p.value,
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
      filter(.data$p.adj <= 0.05)
    
    dunn_results <- alt_input %>%
      group_by(alt_input[grouped_by]) %>%
      filter(.data$test_id %in% kruskal_results$test_id) %>%
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
                           second_group,
                           mean_value,
                           dunn_test){
  filtered_table %>% 
    group_by(.data[[first_group]], .data[[second_group]]) %>%
    summarise(mean = mean(.data[[mean_value]])) -> mean_table
  
  dunn_test %>% 
    merge(mean_table, 
          by.x = c('group1',
                   second_group),
          by.y = c(first_group,
                   second_group)) %>%
    rename_with(~paste0('group1_', mean_value, recycle0 = TRUE), contains('mean')) %>% 
    merge(mean_table,
          by.x = c('group2',
                   second_group),
          by.y = c(first_group,
                   second_group)) %>%
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
    geom_tile(aes(fill = stat_diff_means), alpha = 1, color = 'black') +
    scale_fill_gradient2(low = 'blue', high = 'green', name = 'Group 1 -\nGroup 2') +
    geom_text(aes(label = p.adj.signif)) +
    facet_wrap(~Family,
               nrow = 2) +
    theme_bw(base_size = 20) +
    ggtitle(label = title,
            subtitle = 'Microbe Genus') +
    theme(strip.text.y = element_text(angle = 0),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_discrete(labels = c('New Anschutz (2024)',
                                'U of Arizona')) +
    scale_x_discrete(labels = c('Old Anschutz (2020)',
                                'New Anschutz (2024)')) +
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

## performing statistical analysis
family_abun_lm <- linear_model(input_table = abun_filt,
                               grouped_by = wanted_level,
                               adjust_method = 'BH',
                               filter_adj_p_value = FALSE,
                               formula_left = 'rel_abund',
                               formula_right = 'experiment_set * vendor')

abun_stats <- kruskal_dunn_stats(input_table = abun_filt,
                                 grouped_by = wanted_level,
                                 adjust_method = 'BH',
                                 filter_adj_p_value = FALSE,
                                 formula_left = 'rel_abund',
                                 formula_right = 'experiment_set')

abun_kruskal_test <- abun_stats$KruskalTest
abun_dunn_test <- abun_stats$DunnTest

new_dunn_test <- stat_plot_prep(filtered_table = abun_filt,
                                first_group = 'experiment_set',
                                second_group = wanted_level,
                                mean_value = 'rel_abund',
                                dunn_test = abun_dunn_test)

abun_stat_vis <- stat_plot(new_dunn_test,
                           "All Exp Microbe Relative Abundance")


## saving my outputs as a .tsv
write_tsv(family_abun_lm,
          args$lm_FP)
write_tsv(new_dunn_test,
          args$dunn_FP)

## saving statistical visualization
ggsave(args$stat_plot_FP,
       plot = abun_stat_vis, 
       width = 14, 
       height = 8)
