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
## alpha diversity plot function
alpha_div_plot <- function(alpha_div_fp,
                           rename_me,
                           metadata_file,
                           labels,
                           names_labels,
                           x_axis,
                           y_axis,
                           box_grouped_by,
                           line_grouped_by,
                           facet_by,
                           title,
                           x_label,
                           y_label){
  alpha_table <- read_tsv(alpha_div_fp)
  names(alpha_table)[names(alpha_table) == rename_me] <- 'sampleid'
  metadata_file %>% 
    left_join(alpha_table, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> joint_table
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  joint_table %>%
    ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
    geom_boxplot(aes(group = .data[[box_grouped_by]]), outlier.shape = NA) +
    geom_line(aes(group = .data[[line_grouped_by]]), alpha = 0.1) +
    geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_smooth(se = FALSE) +
    theme_bw(base_size = 20) +
    facet_grid(~.data[[facet_by]], 
               labeller = labeller(.cols = labs) ) +
    ggtitle(title) +
    xlab(x_label) +
    ylab(y_label) -> plot
  return(plot)
}

## core metrics file prep
## metadata file prep
metadata <- read_tsv(args$metadata_FP)

## faith's pd plot 
faith_title <- "Faith's Phylogenetic Diversity"
faith_renamed <- '#SampleID'
faith_x_lab <- 'Days Relative to Infection'
faith_y_lab <- "Faith's PD"

faith_plot <- alpha_div_plot(alpha_div_fp = args$faith_pd_FP,
                             rename_me = faith_renamed,
                             metadata_file = metadata,
                             labels = diet_labs,
                             names_labels = diet_names_labels,
                             x_axis = 'day_post_inf',
                             y_axis = 'faith_pd',
                             box_grouped_by = 'day_post_inf',
                             line_grouped_by = 'mouse_id',
                             facet_by = 'diet',
                             title = faith_title,
                             x_label = faith_x_lab,
                             y_label = faith_y_lab)

## shannon entropy plot 
shannon_title <- "Shannon Entropy"
shannon_renamed <- '...1'
shannon_x_lab <- 'Days Relative to Infection'
shannon_y_lab <- 'Shannon Entropy'

shannon_entropy_plot <- alpha_div_plot(alpha_div_fp = args$shannon_FP,
                                       rename_me = shannon_renamed,
                                       metadata_file = metadata,
                                       labels = diet_labs,
                                       names_labels = diet_names_labels,
                                       x_axis = 'day_post_inf',
                                       y_axis = 'shannon_entropy',
                                       box_grouped_by = 'day_post_inf',
                                       line_grouped_by = 'mouse_id',
                                       facet_by = 'diet',
                                       title = shannon_title,
                                       x_label = shannon_x_lab,
                                       y_label = shannon_y_lab)

## saving my plot outputs to the plots folder
ggsave(args$output_faith_FP,
       plot = faith_plot, 
       width = 14, 
       height = 4.5)

ggsave(args$output_shannon_FP,
       plot = shannon_entropy_plot, 
       width = 14, 
       height = 4.5)