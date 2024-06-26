## 7-14-23
## putting together my microbe family abundance plots based on the 16S rDNA data
## a different way of visualizing the taxa barplot (might be easier)
## you have two options

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
parser$add_argument("-o",
                    "--otu",
                    dest = "otu_table_FP",
                    help = "Filepath to OTU Table file in .qza format.")
parser$add_argument("-t",
                    "--taxonomy",
                    dest = "tax_FP",
                    help = "Filepath to taxonomy file in .qza format.")
parser$add_argument("-p1",
                    "--plot1",
                    dest = "plot1_FP",
                    help = "Filepath to first family abundance plot in .pdf format.")
parser$add_argument("-p2",
                    "--plot2",
                    dest = "plot2_FP",
                    help = "Filepath to second family abundance plot in .pdf format.")

args <- parser$parse_args()



wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Bacteroidaceae', 'Ruminococcaceae')

## functions in order of usage 
## 1
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
             seq_depth, .data[[tax_level]]) %>% 
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
abun_plots <- function(input_table,
                       x_axis,
                       y_axis,
                       x_group_by,
                       fill_by,
                       facet_by,
                       title,
                       x_name,
                       y_name){
  ## first plot
  input_table %>%
    ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
    scale_y_continuous(trans = 'log10') +
    geom_boxplot(aes(group = .data[[x_group_by]]), outlier.shape = NA) +
    geom_jitter(aes(fill = .data[[fill_by]]), width = 0.1, height = 0, alpha = 0.7, pch = 21, size = 2) +
    scale_fill_brewer(palette = 'Pastel1', 
                      name = 'Vendor',
                      labels = c('Charles River',
                                 'Taconic')) +
    theme_bw(base_size = 20) +
    facet_grid(~.data[[facet_by]]) +
    ggtitle(title) +
    ylab(y_name) +
    xlab(x_name) -> plot
  
  return(plot)
}


## family abundance table prep 
abun_files <- family_abun_file_prep(args$metadata_FP,
                                    args$tax_FP,
                                    args$otu_table_FP,
                                    wanted_level,
                                    wanted_family)

exp_x_labs <- c('New Anschutz (2024)',
                'U of Arizona')

exp_labs <- c('New Anschutz (2024)',
              'U of Arizona')

names(exp_labs) <- c('new_exp_anschutz',
                     'second_set_arizona')

## pulling the abundance table out, you can also take metadata, otu table, and taxonomic info
## out of the list output 
abun_filt <- abun_files$AbundanceTable

## generating the plots
abun1 <- abun_plots(input_table = abun_filt,
                    x_axis = 'Family',
                    y_axis = 'rel_abund',
                    x_group_by = 'Family',
                    fill_by = 'vendor',
                    facet_by = 'experiment_set',
                    title = 'New Exp v AZ Exp Microbes at Day -15',
                    x_name = 'Family',
                    y_name = 'Relative Abundance (log10)')

abun1 <- abun1 +
          theme(strip.text.y = element_text(angle = 0),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(~experiment_set,
                     labeller = labeller(experiment_set = exp_labs))


## scatter plot faceted by wanted microbe families with diet on the x-axis for visual comparisons by diet
abun2 <- abun_plots(input_table = abun_filt,
                    x_axis = 'experiment_set',
                    y_axis = 'rel_abund',
                    x_group_by = 'experiment_set',
                    fill_by = 'vendor',
                    facet_by = 'Family',
                    title = 'New Exp v AZ Exp Microbes at Day -15',
                    x_name = 'Experiment',
                    y_name = 'Relative Abundance (log10)')

abun2 <- abun2 +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = exp_x_labs)


## saving my plot outputs
## option #1
ggsave(args$plot1_FP,
       plot = abun1, 
       width = 18, 
       height = 7)
## option #2
ggsave(args$plot2_FP,
       plot = abun2, 
       width = 20, 
       height = 7)
