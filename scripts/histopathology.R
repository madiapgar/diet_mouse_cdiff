## 7-17-23
## this script takes the histopathology scores for all day 3 mice and creates a
## plot along with running a linear model on the results 

## needed libraries
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(rstatix)
library(ggpubr)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-hi",
                    "--histo",
                    dest = "histo_FP",
                    help = "Filepath to histopathology score file in .csv format.")
parser$add_argument("-hp",
                    "--histo_plot",
                    dest = "histo_plot_FP",
                    help = "Filepath to histopathology score plot in .pdf format.")
parser$add_argument("-lm",
                    "--linear_model",
                    dest = "lm_FP",
                    help = "Filepath to histopathology linear modeling results in .tsv format.")
parser$add_argument("-d",
                    "--dunn",
                    dest = "dunn_FP",
                    help = "Filepath to histopathology Dunns Post Hoc test results in .tsv format.")

args <- parser$parse_args()

## input file paths and others
# metadata_FP <- './data/misc/processed_metadata.tsv'
# histo_FP <- './data/misc/histo_data.csv'
tissue_labs <- c('Cecum',
                 'Colon')
names(tissue_labs) <- c('cecum',
                        'colon')

## functions in order of usage 
## 1
histo_file_prep <- function(metadata_fp,
                            histo_fp){
  ## reading in metadata file
  metadata <- read_tsv(metadata_fp)
  ## reading in histopathology scores
  histo <- read_csv(histo_fp) %>% 
    filter(!is.na(mouse_id))
  ## joining the two together 
  metadata %>% 
    merge(histo, by = 'mouse_id') %>% 
    group_by(mouse_id) %>% 
    filter(day_post_inf == max(day_post_inf)) %>% 
    ungroup() %>% 
    mutate(day_post_inf = as.factor(day_post_inf)) %>% 
    gather(cecum, colon, key = tissue, value = score) -> big_histo
  return(big_histo)
}

## 2
histo_stats <- function(big_histo){
  ## kruskal-wallis test
  big_histo %>% 
    group_by(tissue) %>% 
    do(tidy(kruskal.test(score ~ diet,
                         data = .))) -> kruskal
  ## dunn's post hoc test
  big_histo %>% 
    group_by(tissue) %>% 
    dunn_test(score ~ diet,
              p.adjust.method = 'BH',
              data =.) %>% 
    add_y_position(scales = 'free_y', step.increase = 0) -> dunn
  ## linear model
  big_histo %>% 
    group_by(tissue) %>% 
    do(tidy(lm(score ~ purified_diet + high_fat * high_fiber,
               data =.))) %>% 
    na.omit() %>% 
    filter(term != '(Intercept)') %>% 
    arrange(p.value) -> linear_model
  ## creating a list of my outputs
  my_list <- list(KruskalWallis = kruskal,
                  DunnsPostHoc = dunn,
                  LinearModel = linear_model)
  return(my_list)
}

## 3
histo_plot <- function(big_histo,
                       histo_dunn){
  big_histo %>% 
    mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
    ggplot(aes(x = diet, y = score)) +
    geom_violin(aes(group = diet),  draw_quantiles = c(0.5)) +
    geom_jitter(alpha = 0.4, width = 0.1, height = 0) +
    scale_x_discrete(labels = c('Chow', 'High Fat/\nHigh Fiber', 'High Fat/\nLow Fiber',
                                'Low Fat/\nHigh Fiber', 'Low Fat/\nLow Fiber')) +
    facet_wrap(~tissue, labeller = labeller(tissue = tissue_labs),
               scales = "free_y") +
    stat_pvalue_manual(histo_dunn,
                       tip.length = 0.01,
                       label = 'p.adj.signif',
                       hide.ns = TRUE,
                       step.increase = 0.1) +
    theme_bw(base_size = 14) +
    xlab('Diet') +
    ylab('Histopathology Score') +
    ggtitle("Histopathology Score by Diet") -> all_day_plot
  return(all_day_plot)
}


## file prep 
histo <- histo_file_prep(args$metadata_FP,
                         args$histo_FP)

## stats 
stats <- histo_stats(big_histo = histo)

kruskal <- stats$KruskalWallis
dunn <- stats$DunnsPostHoc
linear_model <- stats$LinearModel

## plot 
plot <- histo_plot(histo,
                   dunn)

## saving my plot and stats outputs 
ggsave(args$histo_plot_FP, 
       plot = plot,
       width = 10, 
       height = 5)

write_tsv(linear_model,
          args$lm_FP)
write_tsv(dunn,
          args$dunn_FP)
