## 6-26-23 
## Qiime2 core metrics alpha diversity analysis output plot construction
## Faith's PD plot
## Shannon Entropy plot

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(viridis)
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
parser$add_argument("-of",
                    "--output_faith",
                    dest = "output_faith_FP",
                    help = "Filepath to Faith's PD plot in .pdf format.")
parser$add_argument("-os",
                    "--output_shannon",
                    dest = "output_shannon_FP",
                    help = "Filepath to Shannon Entropy plot in .pdf format.")

args <- parser$parse_args()

## input file paths 
# metadata_FP <- './data/misc/processed_metadata.tsv'
# faith_pd_FP <- './data/qiime/core_outputs/faith_pd.tsv'
# shannon_FP <- './data/qiime/core_outputs/shannon_entropy.tsv'

diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')

diet_names_labels <- c('Chow', 
                       'HF/HF', 
                       'HF/LF', 
                       'LF/HF', 
                       'LF/LF')

## needed functions (in order)
## 1
## faith's pd plot 
## assumes that the files is a .tsv
faith_pd_plot <- function(faith_fp,
                          metadata_file,
                          labels,
                          names_labels,
                          title){
  faith <- read_tsv(faith_fp)
  names(faith)[names(faith) == '#SampleID'] <- 'sampleid'
  metadata_file %>% 
    left_join(faith, by = 'sampleid') -> faith_pd
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  faith_pd %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = faith_pd)) +
    geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
    geom_line(aes(group = mouse_id), alpha = 0.1) +
    geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_smooth(se = FALSE) +
    theme_bw(base_size = 16) +
    facet_grid(~diet, labeller = labeller(diet = labs) ) +
    ggtitle(title) +
    xlab('Days Relative to Infection') +
    ylab("Faith's PD") -> faith_plot
  return(faith_plot)
}

## 2
## shannon entropy plot
## assumes that the file is a .tsv
shannon_plot <- function(shannon_fp,
                         metadata_file,
                         labels,
                         names_labels,
                         title){
  shannon <- read_tsv(shannon_fp)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  metadata_file %>% 
    left_join(shannon, by = 'sampleid') -> shannon_entropy
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  shannon_entropy %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = shannon_entropy)) +
    geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
    geom_line(aes(group = mouse_id), alpha = 0.1) +
    geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_smooth(se = FALSE) +
    theme_bw(base_size = 16) +
    facet_grid(~diet, labeller = labeller(diet = labs) ) +
    ggtitle(title) +
    xlab('Days Relative to Infection') +
    ylab("Shannon Entropy") -> shannon_plot
  return(shannon_plot)
}

## core metrics file prep
## metadata file prep
metadata <- read_tsv(args$metadata_FP)

## faith's pd plot 
faith_title <- "Faith's Phylogenetic Diversity"

faith_plot <- faith_pd_plot(args$faith_pd_FP,
                            metadata,
                            diet_labs,
                            diet_names_labels,
                            faith_title)

## shannon entropy plot 
shannon_title <- "Shannon Entropy"

shannon_entropy_plot <- shannon_plot(args$shannon_FP,
                                     metadata,
                                     diet_labs,
                                     diet_names_labels,
                                     shannon_title)

## saving my plot outputs to the plots folder
ggsave(args$output_faith_FP,
       plot = faith_plot, 
       width = 14, 
       height = 4.5)

ggsave(args$output_shannon_FP,
       plot = shannon_entropy_plot, 
       width = 14, 
       height = 4.5)