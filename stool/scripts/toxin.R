## 9-21-23
## C. diff toxin abundances from mouse cecum at day 3

## needed libraries
library(broom)
library(cowplot)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(glue)
library(argparse)

## argparse for input file paths
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-t",
                    "--toxin",
                    dest = "toxin_FP",
                    help = "Filepath to toxin data file in .tsv format.")
parser$add_argument("-n",
                    "--neat_plot",
                    dest = "neat_plot_FP",
                    help = "Filepath to where the neat toxin concentration plot will go in .pdf format.")
parser$add_argument("-d",
                    "--diluted_plot",
                    dest = "dil_plot_FP",
                    help = "Filepath to where the diluted toxin concentration plot will go in .pdf format.")
parser$add_argument("-nk",
                    "--neat_kruskal",
                    dest = "neat_kruskal_FP",
                    help = "Filepath to where the neat toxin concentration Kruskal Wallis test results will go in .tsv format.")
parser$add_argument("-nd",
                    "--neat_dunn",
                    dest = "neat_dunn_FP",
                    help = "Filepath to where the neat toxin concentration Dunn's Post Hoc test results will go in .tsv format.")
parser$add_argument("-dk",
                    "--diluted_kruskal",
                    dest = "dil_kruskal_FP",
                    help = "Filepath to where the diluted toxin concentration Kruskal Wallis test results will go in .tsv format.")
parser$add_argument("-dd",
                    "--diluted_dunn",
                    dest = "dil_dunn_FP",
                    help = "Filepath to where the diluted toxin concentration Dunn's Post Hoc test results will go in .tsv format.")

args <- parser$parse_args()


## needed input file paths
# metadata_FP <- './data/misc/processed_metadata.tsv'
# toxin_FP <- './data/misc/toxin_final_data.tsv'

## labeling lists
neat_labs <- c('TcdA', 'TcdB')
names(neat_labs) <- c('Total TcA Neat', 'Total TcB Neat')
neat_x_labs <- c('Chow', 
                 'HFt/\nHFb', 
                 'HFt/\nLFb',
                 'LFt/\nHFb', 
                 'LFt/\nLFb')
neat_title <- 'Toxin Neat Concentration by Mouse Diet'

dil_labs <- c('TcdA', 'TcdB')
names(dil_labs) <- c('Total TcA 1:10', 'Total TcB 1:10')
dil_x_labs <- c('HFt/\nHFb', 
                'HFt/\nLFb',
                'LFt/\nHFb', 
                'LFt/\nLFb')
dil_title <- 'Toxin Diluted Concentration (1:10) by Mouse Diet'

## needed functions (in order)
## 1
## for toxin file prep
file_prep <- function(metadata_fp,
                      toxin_fp){
  ## metadata file
  metadata <- read_tsv(metadata_fp)
  ## toxin file
  toxin <- read_tsv(toxin_fp)
  wanted_toxin_ids <- toxin$mouse_id
  metadata %>% 
    group_by(mouse_id) %>% 
    filter(mouse_id %in% wanted_toxin_ids) -> meta_filt
  toxin %>% 
    left_join(meta_filt, by = 'mouse_id') %>% 
    gather('Total TcA Neat', 'Total TcB Neat', 
           key = neat_toxin, value = neat_conc) %>% 
    gather('Total TcA 1:10', 'Total TcB 1:10',
           key = dil_toxin, value = dil_conc) -> pre_toxin
  
  pre_toxin$neat_conc[pre_toxin$neat_conc == 'BDL'] <- '0'
  pre_toxin$dil_conc[pre_toxin$dil_conc == 'BDL'] <- '0'
  pre_toxin$dil_conc[pre_toxin$dil_conc == 'Chow'] <- '0'
  
  pre_toxin %>% 
    filter(!is.na(dil_conc),
           !is.na(neat_conc)) %>%
    mutate(neat_conc = as.numeric(neat_conc),
           dil_conc = as.numeric(dil_conc)) %>% 
    select(-'Extra_Sample', -'Tube_Label', -'Collection Date',
           -'Sample_Type') -> big_toxin
  big_toxin %>% 
    select(-neat_toxin, -neat_conc) %>% 
    distinct(mouse_id, dil_conc, .keep_all = TRUE) %>% 
    filter(diet != 'Chow') -> dil_toxin
  
  big_toxin %>%
    select(-dil_toxin, -dil_conc) %>% 
    distinct(mouse_id, neat_conc, .keep_all = TRUE) -> neat_toxin 
  ## creating a list of my outputs
  my_list <- list(Metadata = metadata,
                  Toxin = big_toxin,
                  NeatToxin = neat_toxin,
                  DilToxin = dil_toxin)
  return(my_list)
}

## 2
## statistical analysis function
stats <- function(biom_table,
                  tox_col,
                  conc_col){
  ## kruskal test
  biom_table %>% 
    group_by(biom_table[tox_col]) %>% 
    do(tidy(kruskal.test(.data[[conc_col]] ~ diet,
                         data = .))) %>% 
    ungroup() %>%
    arrange(p.value) %>%
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(.data[[tox_col]])) -> kruskal
  ## dunns post hoc test
  diet_name <- 'diet'
  test <- reformulate(glue("{diet_name}"),glue("{conc_col}"))
  
  biom_table %>% 
    group_by(biom_table[tox_col]) %>% 
    mutate(test_id = paste(.data[[tox_col]])) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(test,
              p.adjust.method = 'BH',
              data = .) %>% 
    add_y_position(scales = 'free_y') -> dunn
  ## linear modeling
  biom_table %>% 
    group_by(biom_table[tox_col]) %>%
    do(glance(lm(.data[[conc_col]] ~ (purified_diet * seq_depth) + high_fat * high_fiber,
                 data = .))) %>% 
    ungroup() %>%
    na.omit() %>%
    mutate(adj.p = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(.data[[tox_col]])) -> lm_full
  biom_table %>% 
    group_by(biom_table[tox_col]) %>% 
    mutate(test_id = paste(.data[[tox_col]])) %>% 
    filter(test_id %in% lm_full$test_id) %>%
    do(tidy(lm(.data[[conc_col]] ~ (purified_diet * seq_depth) + high_fat * high_fiber,
               data = .))) %>% 
    filter(term != '(Intercept)') -> lm
  
  lm['signif'] <- symnum(lm$p.value,
                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("****", "***", "**", "*", "+", "ns"),
                         abbr.colnames = FALSE,
                         na = "")
  ## creating a list of my outputs
  my_list <- list(KruskalTest = kruskal,
                  DunnPostHoc = dunn,
                  LinearModel = lm)
  return(my_list)
}

## 3 
## toxin ggplot function
tox_plot <- function(biom_table,
                     tox_col,
                     conc_col,
                     x_labels,
                     facet_labs,
                     dunn,
                     title){
  biom_table %>% 
    na.omit() %>% 
    ggplot(aes(x = diet, y = .data[[conc_col]])) +
    geom_violin(aes(group = diet), draw_quantiles = c(0.5)) +
    geom_jitter(alpha = 0.4, width = 0.1, height = 0)+
    scale_x_discrete(labels = x_labels) +
    facet_wrap(~.data[[tox_col]],
               labeller = labeller(.cols = facet_labs),
               scales = 'free_y') +
    stat_pvalue_manual(dunn,
                       tip.length = 0.01,
                       label = 'p.adj.signif',
                       hide.ns = TRUE) +
    theme_bw(base_size = 14) +
    xlab('Diet') +
    ylab('Concentration') +
    ggtitle(title) -> plot
  return(plot)
}

## file prep
toxin_files <- file_prep(args$metadata_FP,
                         args$toxin_FP)

metadata <- toxin_files$Metadata
neat_tox_table <- toxin_files$NeatToxin
dil_tox_table <- toxin_files$DilToxin

## statistical analysis
## neat
neat_stats <- stats(neat_tox_table,
                    "neat_toxin",
                    "neat_conc")

neat_dunn <- neat_stats$DunnPostHoc
neat_kruskal <- neat_stats$KruskalTest
neat_lm <- neat_stats$LinearModel

## diluted
dil_stats <- stats(dil_tox_table,
                   "dil_toxin",
                   "dil_conc")

dil_dunn <- dil_stats$DunnPostHoc
dil_kruskal <- dil_stats$KruskalTest
dil_lm <- dil_stats$LinearModel

## plots
## neat 
neat_plot <- tox_plot(neat_tox_table,
                      "neat_toxin",
                      "neat_conc",
                      neat_x_labs,
                      neat_labs,
                      neat_dunn,
                      neat_title)
## diluted
dil_plot <- tox_plot(dil_tox_table,
                      "dil_toxin",
                      "dil_conc",
                      dil_x_labs,
                      dil_labs,
                      dil_dunn,
                      dil_title)
## saving my outputs
## plots
ggsave(args$neat_plot_FP,
       plot = neat_plot,
       width = 10,
       height = 5)

ggsave(args$dil_plot_FP,
       plot = dil_plot,
       width = 10,
       height = 5)

## stats
write_tsv(neat_kruskal,
          args$neat_kruskal_FP)
write_tsv(neat_dunn,
          args$neat_dunn_FP)
write_tsv(dil_kruskal,
          args$dil_kruskal_FP)
write_tsv(dil_dunn,
          args$dil_dunn_FP)