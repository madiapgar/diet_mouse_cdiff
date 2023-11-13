## 7-14-23
## creating a script that runs linear modeling statistics on the picrust 
## taxon functional abundance outputs for butyrate and secondary bile acid enzymes

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
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-k",
                    "--ko",
                    dest = "ko_contrib_FP",
                    help = "Filepath to KO metagenome contrib file in .tsv format.")
parser$add_argument("-t",
                    "--taxonomy",
                    dest = "tax_FP",
                    help = "Filepath to taxonomy file in .qza format.")
parser$add_argument("-bt",
                    "--buty_lm",
                    dest = "butyrate_lm_FP",
                    help = "Filepath to butyrate enzyme linear model results in .tsv format.")
parser$add_argument("-bd",
                    "--buty_dunn",
                    dest = "butyrate_dunn_FP",
                    help = "Filepath to butyrate enzyme Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-ba",
                    "--bile_lm",
                    dest = "bile_lm_FP",
                    help = "Filepath to bile acid enzyme linear model results in .tsv format.")
parser$add_argument("-bid",
                    "--bile_dunn",
                    dest = "bile_dunn_FP",
                    help = "Filepath to bile acid enzyme Dunn's Post Hoc test results in .tsv format.")
parser$add_argument("-bs",
                    "--buty_stat_vis",
                    dest = "buty_stat_FP",
                    help = "Filepath to butyrate enzyme statistical visualization in .pdf format.")
parser$add_argument("-bis",
                    "--bile_stat_vis",
                    dest = "bile_stat_FP",
                    help = "Filepath to bile acid enzyme statistical visualization in .pdf format.")

args <- parser$parse_args()

## input file paths and KOs
# metadata_FP <- './cecum/data/misc/cecal_processed_metadata.tsv'
# ko_contrib_FP <- './cecum/data/picrust/meta_contrib.tsv'
# tax_FP <- './cecum/data/cecal_qiime/taxonomy.qza'
but_kos <- c('K00929','K01034')
bile_kos <- c('K15873', 'K15874')

buty_ko_labs <- c('Butyrate\nKinase',
                  'Butyryl-CoA\nTransferase')

names(buty_ko_labs) <- c('K00929',
                         'K01034')

bile_ko_labs <- c('baiH',
                  'baiI')

names(bile_ko_labs) <- c('K15873', 
                         'K15874')

## needed functions
## 1
stat_file_prep <- function(metadata_fp,
                           ko_contrib_fp,
                           tax_fp,
                           ko_list,
                           taxonomy_level){
  ## taxonomy file
  read_qza(file = tax_fp)$data %>% 
    parse_taxonomy() %>% 
    as_tibble(rownames = 'taxon') %>% 
    rename_all(tolower) -> taxonomy
  ## metadata file
  metadata <- read_tsv(metadata_fp)
  names(metadata)[names(metadata) == 'sampleid'] <- 'sample'
  ## ko meta contrib file
  kometa_contrib <- read_tsv(file = ko_contrib_fp)
  names(kometa_contrib)[names(kometa_contrib) == 'function'] <- 'ko'
  ## putting all components together into one giant table
  kometa_contrib %>% 
    left_join(taxonomy, by = 'taxon') %>% 
    left_join(metadata, by = 'sample') -> stat_biom
  ## messing with biom table format so that the zeroes are represented
  stat_biom %>% 
    select(sample, ko, taxon_function_abun, study, diet, facility, purified_diet, high_fat,
           high_fiber, mouse_id, seq_depth, any_of(taxonomy_level)) %>% 
    filter(ko %in% ko_list) %>% 
    group_by(ko, sample, diet, study, facility, purified_diet, high_fat, high_fiber, mouse_id, seq_depth,
             .data[[taxonomy_level]]) %>% 
    summarise(taxon_function_abun = sum(taxon_function_abun)) %>% 
    spread(diet, taxon_function_abun, fill = 0) %>% 
    gather(-ko, -sample, -facility, -study, -.data[[taxonomy_level]], -purified_diet, -high_fat, 
           -high_fiber, -mouse_id, -seq_depth,
           key = diet, value = taxon_function_abun) -> biom_long
  return(biom_long)
}

## 2
## function for running statisical analysis on picrust kos of interest 
stat_calc <- function(biom_table,
                           taxonomy_level){
  ## linear modeling
  biom_table %>% 
    na.omit() %>% 
    group_by(ko, .data[[taxonomy_level]]) %>% 
    do(glance(lm(taxon_function_abun ~ (purified_diet * seq_depth) + high_fat * high_fiber,
                 data = .))) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(adj.p = p.adjust(p.value, 
                            method = "BH"),
           test_id = paste(ko, .data[[taxonomy_level]], sep = "_")) -> lm_full
    # filter(adj.p <= 0.05) -> lm_full
  biom_table %>% 
    na.omit() %>% 
    group_by(ko, .data[[taxonomy_level]]) %>% 
    mutate(test_id = paste(ko, .data[[taxonomy_level]], sep = "_")) %>% 
    filter(test_id %in% lm_full$test_id) %>% 
    do(tidy(lm(taxon_function_abun ~ (purified_diet * seq_depth) + high_fat * high_fiber,
               data = .))) %>%
    na.omit() %>% 
    filter(term != '(Intercept)') -> lm
  lm['signif'] <- symnum(lm$p.value,
                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("****", "***", "**", "*", "+", "ns"),
                         abbr.colnames = FALSE,
                         na = "")
  ## kruskal wallis and dunns post hoc tests
  biom_table %>% 
    na.omit() %>% 
    group_by(ko, .data[[taxonomy_level]]) %>% 
    do(tidy(kruskal.test(taxon_function_abun ~ diet,
                         data = .))) %>% 
    ungroup() %>% 
    arrange(p.value) %>% 
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(ko, .data[[taxonomy_level]], sep = "_")) -> kruskal
  biom_table %>% 
    na.omit() %>% 
    group_by(ko, .data[[taxonomy_level]]) %>% 
    mutate(test_id = paste(ko, .data[[taxonomy_level]], sep = "_")) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(taxon_function_abun ~ diet,
              p.adjust.method = 'BH',
              data = .) -> dunn
  ## creating a list 
  my_list <- list(LinearModel = lm,
                  KruskalTest = kruskal,
                  DunnPostHoc = dunn)
  return(my_list)
}


## 3
## creating a function for this so I don't have to keep doing each one by hand 
stat_plot_prep <- function(biom_table,
                           dunn_test,
                           taxonomy_level,
                           value){
  biom_table %>% 
    group_by(diet, .data[[taxonomy_level]]) %>% 
    summarise(means = mean(.data[[value]])) -> mean_table
  dunn_test %>% 
    merge(mean_table, 
          by.x = c('group1',
                   taxonomy_level),
          by.y = c('diet',
                   taxonomy_level)) %>% 
    rename('group1_means' = 'means') %>% 
    merge(mean_table,
          by.x = c('group2',
                   taxonomy_level),
          by.y = c('diet',
                   taxonomy_level)) %>% 
    rename('group2_means' = 'means') %>% 
    mutate(diff_means = (group1_means - group2_means),
           stat_diff_means = if_else(p.adj > 0.05, 0, diff_means)) -> new_dunn
  return(new_dunn)
}

## 4
stat_plot <- function(new_dunn,
                      taxonomy_level,
                      labs,
                      title){
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
    facet_grid(ko~.data[[taxonomy_level]],
               scales = 'free_x',
               labeller = labeller(ko = labs)) +
    theme_bw(base_size = 16) +
    theme(strip.text.y = element_text(angle = 0)) +
    xlab('Group 1') +
    ylab('Group 2') +
    ggtitle(title) -> stat_vis
  return(stat_vis)
}

## file prep 
## butyrate 
but_long <- stat_file_prep(args$metadata_FP,
                           args$ko_contrib_FP,
                           args$tax_FP,
                           but_kos,
                           'class')

## bile acids
bile_long <- stat_file_prep(args$metadata_FP,
                            args$ko_contrib_FP,
                            args$tax_FP,
                            bile_kos,
                            'genus')

## butyrate statistical analysis
buty_stats <- stat_calc(but_long,
                        'class')

buty_lm <- buty_stats$LinearModel
buty_kruskal <- buty_stats$KruskalTest
buty_dunn <- buty_stats$DunnPostHoc

new_buty_dunn <- stat_plot_prep(but_long,
                                buty_dunn,
                                'class',
                                'taxon_function_abun')

buty_plot <- stat_plot(new_buty_dunn,
                       'class',
                       buty_ko_labs,
                       'Cecal Butyrate Functional Potential at Day 3')

## bile acid statistical analysis
bile_stats <- stat_calc(bile_long,
                        'genus')

bile_lm <- bile_stats$LinearModel
bile_kruskal <- bile_stats$KruskalTest
bile_dunn <- bile_stats$DunnPostHoc

new_bile_dunn <- stat_plot_prep(bile_long,
                                bile_dunn,
                                'genus',
                                'taxon_function_abun')

bile_plot <- stat_plot(new_bile_dunn,
                       'genus',
                       bile_ko_labs,
                       'Cecal Bile Acid Functional Potential at Day 3')

## saving my outputs as a .tsv
write_tsv(buty_lm,
          args$butyrate_lm_FP)
write_tsv(new_buty_dunn,
          args$butyrate_dunn_FP)
write_tsv(bile_lm,
          args$bile_lm_FP)
write_tsv(new_bile_dunn,
          args$bile_dunn_FP)

ggsave(args$buty_stat_FP,
       plot = buty_plot, 
       width = 15, 
       height = 5)
ggsave(args$bile_stat_FP,
       plot = bile_plot, 
       width = 12, 
       height = 4.5)
