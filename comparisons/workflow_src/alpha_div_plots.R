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
library(vegan)
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


## needed functions (in order)
## 1
make_plot <- function(input_table,
                      metadata_file,
                      x_axis,
                      y_axis,
                      x_group_by,
                      x_labels,
                      fill_by,
                      x_name,
                      y_name,
                      title){
  metadata_file %>% 
    left_join(input_table, by = 'sampleid') -> plot_table
  
  plot_table %>%
    ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
    geom_boxplot(aes(group = .data[[x_group_by]]), outlier.shape = NA) +
    geom_jitter(aes(fill = .data[[fill_by]]), width = 0.1, height = 0, alpha = 0.7, size = 2, pch = 21) +
    scale_fill_brewer(palette = 'Pastel1', 
                      name = 'Vendor',
                      labels = c('Charles River',
                                 'Taconic')) +
    scale_x_discrete(labels = x_labels) +
    theme_bw(base_size = 20) +
    ggtitle(title) +
    xlab(x_name) +
    ylab(y_name) -> plot
  return(plot)
}

## core metrics file prep
## reading in files
metadata <- read_tsv(args$metadata_FP)
names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'

faith <- read_tsv(args$faith_pd_FP)
names(faith)[names(faith) == '#SampleID'] <- 'sampleid'

shannon <- read_tsv(args$shannon_FP)
names(shannon)[names(shannon) == '...1'] <- 'sampleid'

exp_x_labs <- c('Old Anschutz (2020)',
                'New Anschutz (2024)',
                'U of Arizona')

## faith's pd plot 
faith_plot <- make_plot(input_table = faith,
                        metadata_file = metadata,
                        x_axis = 'experiment_set',
                        y_axis = 'faith_pd',
                        x_group_by = 'experiment_set',
                        x_labels = exp_x_labs,
                        fill_by = 'vendor',
                        x_name = 'Experiment',
                        y_name = "Faith's PD",
                        title = "All Exp Faith's Phylogenetic Diversity")

## shannon entropy plot 
shannon_plot <- make_plot(input_table = shannon,
                          metadata_file = metadata,
                          x_axis = 'experiment_set',
                          y_axis = 'shannon_entropy',
                          x_group_by = 'experiment_set',
                          x_labels = exp_x_labs,
                          fill_by = 'vendor',
                          x_name = 'Experiment',
                          y_name = "Shannon Entropy",
                          title = "All Exp Shannon Entropy")

## saving my plot outputs to the plots folder
ggsave(args$output_faith_FP,
       plot = faith_plot, 
       width = 12, 
       height = 7)

ggsave(args$output_shannon_FP,
       plot = shannon_plot, 
       width = 12, 
       height = 7)