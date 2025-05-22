## 6-26-23 
## Qiime2 core metrics diversity analysis output plot construction
## unweighted and weighted UniFrac PCoA plots

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(rstatix)
library(ape)
library(vegan)
library(ggh4x)
library(viridis)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-uu",
                    "--unweighted_unifrac",
                    dest = "unweighted_FP",
                    help = "Filepath to Unweighted UniFrac PCoA in .qza format.")
parser$add_argument("-wu",
                    "--weighted_unifrac",
                    dest = "weighted_FP",
                    help = "Filepath to Weighted UniFrac PCoA in .qza format")
parser$add_argument("-f",
                    "--faith_pd",
                    dest = "faith_pd_FP",
                    help = "Filepath to Faith's PD file in .tsv format.")
parser$add_argument("-s",
                    "--shannon",
                    dest = "shannon_FP",
                    help = "Filepath to Shannon Entropy file in .tsv format.")
parser$add_argument("-ou",
                    "--output_uu",
                    dest = "output_uu_FP",
                    help = "Filepath to Unweighted UniFrac PCoA plot in .pdf format.")
parser$add_argument("-ow",
                    "--output_wu",
                    dest = "output_wu_FP",
                    help = "Filepath to Weighted UniFrac PCoA plot in .pdf format.")

args <- parser$parse_args()


## input file paths
# metadata_FP <- './data/misc/processed_metadata.tsv'
# unweighted_FP <- './data/qiime/core_outputs/unweighted_unifrac_pcoa_results.qza'
# weighted_FP <- './data/qiime/core_outputs/weighted_unifrac_pcoa_results.qza'
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

## functions in order of usage
## 1 
## unweighted/weighted unifrac pcoa result, faith's pd, and shannon entropy file prep 
## going to attempt to return multiple outputs so I can just have one function for file prep
biom_table_prep <- function(unweighted_fp,
                            weighted_fp,
                            faith_fp,
                            shannon_fp,
                            metadata_file){
  ## unweighted pcoa
  unweighted <- read_qza(unweighted_fp)$data
  unweighted_var <- unweighted$ProportionExplained
  unweighted_pcoa <- unweighted$Vectors ##used for pcoa plot
  names(unweighted_pcoa)[names(unweighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## weighted pcoa
  weighted <- read_qza(weighted_fp)$data
  weighted_var <- weighted$ProportionExplained
  weighted_pcoa <- weighted$Vectors
  names(weighted_pcoa)[names(weighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## faith's 
  faith <- read_tsv(faith_fp)
  names(faith)[names(faith) == '#SampleID'] <- 'sampleid'
  ## shannon 
  shannon <- read_tsv(shannon_fp)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  ## unweighted biom 
  unweighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> unweighted_biom
  ## weighted biom
  weighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> weighted_biom
  ## creating a list to return multiple outputs 
  my_list <- list(UnweightedVar = unweighted_var, 
                  WeightedVar = weighted_var,
                  UnweightedBiom = unweighted_biom,
                  WeightedBiom = weighted_biom)
  return(my_list)
}

## 2
## this function will pull out the percent variations from a specified column so you can add it to your pcoa plots 
pcoa_ax_lab <- function(unifrac_var, col_name){
  uni_lab <- as.character(round(unifrac_var[col_name] * 100, 2))
  uni_lab <- paste0(col_name, ' - ', uni_lab, '%')
  return(uni_lab)
}

## 3
## pcoa plot function
## xlab and ylab are outputs from pcoa_ax_lab function
pcoa_plot <- function(biom_file,
                      labels,
                      names_labels,
                      xlab,
                      ylab,
                      title){
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  biom_file %>% 
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = faith_pd), pch = 21, alpha = 0.7, size = 2) +
    theme_bw(base_size = 14) +
    #scale_fill_brewer(palette = 'Dark2', name = 'Vendor') +
    scale_fill_distiller(palette = 'Spectral', name = "Faith's PD") +
    facet_grid(day_post_inf~diet, 
               labeller = labeller(diet = labs)) +
    theme(legend.text = element_text(size = 8.5),
          strip.text.y = element_text(angle = 0)) +
    ggtitle(title) +
    labs(x = xlab, y = ylab) -> pcoa
  return(pcoa)
}

## core metrics file prep
## metadata file prep
metadata <- read_tsv(args$metadata_FP)
names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'

## preparing core beta diversity files for ggplot
core_files <- biom_table_prep(args$unweighted_FP,
                              args$weighted_FP,
                              args$faith_pd_FP,
                              args$shannon_FP,
                              metadata)
## extracting core beta diversity files from named list 
uw_var <- core_files$UnweightedVar
w_var <- core_files$WeightedVar
unweighted_biom <- core_files$UnweightedBiom
weighted_biom <- core_files$WeightedBiom

## unweighted pcoa
uw_uni_xlab <- pcoa_ax_lab(uw_var, 'PC1')
uw_uni_ylab <- pcoa_ax_lab(uw_var, 'PC2')

uw_title <- 'Unweighted UniFrac PCoA Plot'

unweighted_pcoa <- pcoa_plot(unweighted_biom,
                             diet_labs,
                             diet_names_labels,
                             uw_uni_xlab,
                             uw_uni_ylab,
                             uw_title)
## weighted pcoa 
w_uni_xlab <- pcoa_ax_lab(w_var, 'PC1')
w_uni_ylab <- pcoa_ax_lab(w_var, 'PC2')

w_title <- 'Weighted UniFrac PCoA Plot'

weighted_pcoa <- pcoa_plot(weighted_biom,
                           diet_labs,
                           diet_names_labels,
                           w_uni_xlab,
                           w_uni_ylab,
                           w_title)

## saving my plot outputs to the plots folder
ggsave(args$output_uu_FP,
       plot = unweighted_pcoa, 
       width = 12, 
       height = 4.5)

ggsave(args$output_wu_FP,
       plot = weighted_pcoa, 
       width = 12, 
       height = 4.5)
