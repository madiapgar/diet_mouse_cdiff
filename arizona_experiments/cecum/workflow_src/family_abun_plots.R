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

## input file paths and others
# otu_table_FP <- './cecum/data/cecal_qiime/tax_filt_actual.qza'
# tax_FP <- './cecum/data/cecal_qiime/taxonomy.qza'
# metadata_FP <- './cecum/data/misc/cecal_processed_meta.tsv'
diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')
names(diet_labs) <- c('Chow', 
                      'HF/HF', 
                      'HF/LF', 
                      'LF/HF', 
                      'LF/LF')
wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Tannerellaceae', 'Muribaculaceae', 'Bacteroidaceae', 
                   'Marinifilaceae', 'Ruminococcaceae')

## functions in order of usage 
## 1
family_abun_file_prep <- function(metadata_fp,
                                  tax_fp,
                                  otu_table_fp,
                                  tax_level,
                                  wanted_tax){
  ## metadata
  metadata <- read_tsv(metadata_fp)
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
    group_by(sampleid, day_post_inf, diet, mouse_id, 
             purified_diet, high_fat, high_fiber, 
             seq_depth, .data[[tax_level]]) %>% 
    summarise(rel_abund = sum(rel_abun)) %>% 
    filter(.data[[tax_level]] %in% wanted_tax) %>% 
    mutate(mouse_fact = as.factor(mouse_id),
           day_fact = as.factor(day_post_inf)) -> abun_filt
  ## creating a list for my outputs
  my_list <- list(Metadata = metadata,
                  Taxonomy = taxonomy,
                  OTUTable = otu_table,
                  AbundanceTable = abun_filt)
  return(my_list)
}

## 2 
abun_plots <- function(abundance_table){
  ## first plot
  abundance_table %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = Family, y = rel_abund)) +
    scale_y_continuous(trans = 'log10') +
    geom_boxplot(aes(group = Family), outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    theme_bw(base_size = 16) +
    facet_grid(~diet, labeller = labeller(diet = diet_labs)) +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Cecal Microbe Family Relative Abundance at Day 3") +
    ylab("Relative Abundance") +
    xlab("Family") -> family_abun1
  ## second plot
  abundance_table %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = diet, y = rel_abund)) +
    scale_y_continuous(trans = 'log10') +
    geom_boxplot(aes(group = diet), outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_discrete(labels = c("Chow",
                                "HFt/HFb",
                                "HFt/LFb",
                                "LFt/HFb",
                                "LFt/LFb")) +
    theme_bw(base_size = 16) +
    facet_wrap(~Family,
               ncol = 1,
               strip.position = "right") +
    theme(strip.text.y = element_text(angle = 0)) +
    ggtitle("Cecal Microbe Family Relative Abundance at Day 3") +
    ylab("Relative Abundance") +
    xlab("Diet") -> family_abun2
  ## creating a list of my two plots
  my_list <- list(FamilyAbundance1 = family_abun1,
                  FamilyAbundance2 = family_abun2)
  return(my_list)
}

## family abundance table prep 
abun_files <- family_abun_file_prep(args$metadata_FP,
                                    args$tax_FP,
                                    args$otu_table_FP,
                                    wanted_level,
                                    wanted_family)

## pulling the abundance table out, you can also take metadata, otu table, and taxonomic info
## out of the list output 
abun_filt <- abun_files$AbundanceTable

## generating the plots
family_abun_plots <- abun_plots(abun_filt)

## scatter plot by wanted microbe families and diet at day 3 in the cecum 
abun1 <- family_abun_plots$FamilyAbundance1
## scatter plot faceted by wanted microbe families with diet on the x-axis for visual comparisons by diet 
abun2 <- family_abun_plots$FamilyAbundance2

## saving my plot outputs
## option #1
ggsave(args$plot1_FP,
       plot = abun1, 
       width = 18, 
       height = 7)
## option #2
ggsave(args$plot2_FP,
       plot = abun2, 
       width = 12, 
       height = 17)
