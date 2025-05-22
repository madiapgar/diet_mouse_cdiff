## 9-26-23
## this script takes metabolomics data and creates plots/does statistical analysis
## on the data 

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
parser$add_argument("-mb",
                    "--metab",
                    dest = "metab_FP",
                    help = "Filepath to metabolomics file in .csv format.")
parser$add_argument("-mp",
                    "--metab_plot",
                    dest = "metab_plot_FP",
                    help = "Filepath to metabolomics output plot in .pdf format.")
parser$add_argument("-mlm",
                    "--metab_lm",
                    dest = "metab_lm_FP",
                    help = "Filepath to metabolomics linear model output in .tsv format.")
parser$add_argument("-md",
                    "--metab_dunn",
                    dest = "metab_dunn_FP",
                    help = "Filepath to metabolomics Dunn's Post Hoc test output in .tsv format.")
parser$add_argument("-mk",
                    "--metab_kruskal",
                    dest = "metab_kruskal_FP",
                    help = "Filepath to metabolomics Kruskal-Wallis test output in .tsv format.")

args <- parser$parse_args()

## input file paths
# metab_FP <- './data/misc/processed_metabolomics.tsv'
metab_labs <- c('Acetic Acid',
                'Propanoic Acid',
                'n-Butanoic Acid')
names(metab_labs) <- c('Acetic Acid (ug/g)',
                       'Propanoic Acid (ug/g)',
                       'n-Butanoic Acid (ug/g)')

metab_x_labs <- c('Chow', 
                 'HFt/\nHFb', 
                 'HFt/\nLFb',
                 'LFt/\nHFb', 
                 'LFt/\nLFb')

metab_title <- 'Cecal Metabolite Concentration by Mouse Diet'

## needed functions
## 1
## statistical analysis function
stats <- function(biom_table,
                  metab_col,
                  conc_col){
  ## kruskal test
  biom_table %>% 
    group_by(biom_table[metab_col]) %>% 
    do(tidy(kruskal.test(.data[[conc_col]] ~ diet,
                         data = .))) %>% 
    ungroup() %>%
    arrange(p.value) %>%
    mutate(p.adj = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(.data[[metab_col]])) -> kruskal
  ## dunns post hoc test
  diet_name <- 'diet'
  test <- reformulate(glue("{diet_name}"),glue("{conc_col}"))
  
  biom_table %>% 
    group_by(biom_table[metab_col]) %>% 
    mutate(test_id = paste(.data[[metab_col]])) %>% 
    filter(test_id %in% kruskal$test_id) %>% 
    dunn_test(test,
              p.adjust.method = 'BH',
              data = .) %>% 
    add_y_position(scales = 'free_y') -> dunn
  ## linear modeling
  biom_table %>% 
    group_by(biom_table[metab_col]) %>%
    do(glance(lm(.data[[conc_col]] ~ (purified_diet * seq_depth) + high_fat * high_fiber,
                 data = .))) %>% 
    ungroup() %>%
    na.omit() %>%
    mutate(adj.p = p.adjust(p.value,
                            method = "BH"),
           test_id = paste(.data[[metab_col]])) -> lm_full
  biom_table %>% 
    group_by(biom_table[metab_col]) %>% 
    mutate(test_id = paste(.data[[metab_col]])) %>% 
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

## 2
## metabolite ggplot function
metab_plot <- function(biom_table,
                       metab_col,
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
    facet_wrap(~.data[[metab_col]],
               labeller = labeller(.cols = facet_labs),
               scales = 'free_y') +
    stat_pvalue_manual(dunn,
                       tip.length = 0.01,
                       label = 'p.adj.signif',
                       hide.ns = TRUE) +
    theme_bw(base_size = 16) +
    xlab('Diet') +
    ylab('Concentration (ug/g)') +
    ggtitle(title) -> plot
  return(plot)
}

## reading in metabolomics file
metab <- read_tsv(args$metab_FP)

## statistical analysis
metab_stats <- stats(metab,
                     "metabolite",
                     "concentration")

kruskal <- metab_stats$KruskalTest
dunn <- metab_stats$DunnPostHoc
linear_model <- metab_stats$LinearModel

## plot
metab_plot <- metab_plot(metab,
                         "metabolite",
                         "concentration",
                         metab_x_labs,
                         metab_labs,
                         dunn,
                         metab_title)

## saving my outputs
## plot
ggsave(args$metab_plot_FP,
       plot = metab_plot,
       width = 15,
       height = 5)

## statistical tests 
write_tsv(linear_model,
          args$metab_lm_FP)
write_tsv(dunn,
          args$metab_dunn_FP)
write_tsv(kruskal,
          args$metab_kruskal_FP)
