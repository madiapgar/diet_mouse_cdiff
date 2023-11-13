## 7-13-23
## script that contains all of the needed functions for baiH and baiI plots 

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
parser$add_argument("-t",
                    "--taxonomy",
                    dest = "tax_FP",
                    help = "Filepath to taxonomy file in .qza format.")
parser$add_argument("-k",
                    "--ko",
                    dest = "ko_contrib_FP",
                    help = "Filepath to KO metagenome contrib file in .tsv format.")
parser$add_argument("-bh",
                    "--baiH_plot",
                    dest = "baiH_plot_FP",
                    help = "Filepath to baiH plot in .pdf format.")
parser$add_argument("-bi",
                    "--baiI_plot",
                    dest = "baiI_plot_FP",
                    help = "Filepath to baiI plot in .pdf format.")

args <- parser$parse_args()

## input file paths 
# metadata_FP <- './cecum/data/misc/cecal_processed_metadata.tsv'
# tax_FP <- './cecum/data/cecal_qiime/taxonomy.qza'
# ko_contrib_FP <- './cecum/data/picrust/meta_contrib.tsv'

## other needed inputs 
diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')
diet_names_labs <- c('Chow',
                     'HF/HF',
                     'HF/LF',
                     'LF/HF',
                     'LF/LF')
## baiH specs
bile_tax_level <- 'genus'
baiH_ko <- 'K15873'
baiH_title <- 'Cecal baiH Potential at Day 3'

## baiI specs
baiI_ko <- 'K15874'
baiI_title <- 'Cecal baiI Potential at Day 3'

## functions in order of usage 
## 1 
## input file prep for bile acid plots 
bile_file_prep <- function(tax_fp,
                           ko_contrib_fp,
                           metadata_fp,
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
    left_join(metadata, by = 'sample') -> kometa_contrib_big
  ## filtering for wanted kos and taxonomic level
  kometa_contrib_big %>% 
    select(sample, ko, taxon_function_abun, study, diet, facility, 
           any_of(taxonomy_level)) %>% 
    filter(ko %in% ko_list) -> filtered_biom
  ## summing the wanted taxonomic level's abundance for a particular ko
  filtered_biom %>% 
    group_by(sample, ko, study, diet, facility, .data[[taxonomy_level]]) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    ungroup() -> filtered_biom_sum
  ## prep for bile acid table
  filtered_biom_sum %>% 
    spread(diet, taxon_function_abun, fill = 0.01) %>% 
    gather(-sample, -ko, -study, -facility, -genus, 
           key = diet, value = taxon_function_abun) -> bile_sum
  ## creating the 'total' facet table to rbind to bile_sum
  filtered_biom_sum %>% 
    group_by(sample, facility, ko, diet, study) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    mutate(genus = 'Total') %>% 
    spread(diet, taxon_function_abun, fill = 0.01) %>% 
    gather(-sample, -ko, -study, -facility, -genus, 
           key = diet, value = taxon_function_abun) -> bile_total
  ## rbinding them together 
  bile_total %>% 
    rbind(bile_sum) -> bile_filtered_sum
  ## creating a list of my outputs
  my_list <- list(Taxonomy = taxonomy,
                  Metadata = metadata,
                  KOContrib = kometa_contrib,
                  ProcessedKOBiom = bile_filtered_sum)
  return(my_list)
}

## 2
## bile acid plot 
bile_plot <- function(processed_ko_biom,
                      labels,
                      names_labels,
                      title){
  labs <- labels
  names(labs) <- names_labels
  processed_ko_biom %>% 
    filter(!is.na(diet)) %>% 
    filter(!is.na(genus)) %>% 
    ggplot(aes(x = diet, y = taxon_function_abun)) +
    scale_y_continuous(trans = 'log10') +
    scale_x_discrete(labels = c('Chow',
                                'HFt/\nHFb',
                                'HFt/\nLFb',
                                'LFt/\nHFb',
                                'LFt/\nLFb')) +
    geom_violin(aes(group = diet)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    geom_smooth(se = FALSE, color = 'blue') +
    facet_grid(~genus) +
    theme_bw(base_size = 16) +
    theme(legend.text = element_text(size = 8.5),
          strip.text.y = element_text(angle = 0)) +
    guides(fill = guide_legend(override.aes = list(size = 2.5))) + 
    ggtitle(title) +
    xlab('Diet') +
    ylab('KO Counts') -> bile_plot
  return(bile_plot)
}

## baiH
## file prep
baiH_files <- bile_file_prep(args$tax_FP,
                             args$ko_contrib_FP,
                             args$metadata_FP,
                             baiH_ko,
                             bile_tax_level)

for_baiH_plot <- baiH_files$ProcessedKOBiom

## plot
baiH <- bile_plot(for_baiH_plot,
                  diet_labs,
                  diet_names_labs,
                  baiH_title)

## baiI
## file prep
baiI_files <- bile_file_prep(args$tax_FP,
                             args$ko_contrib_FP,
                             args$metadata_FP,
                             baiI_ko,
                             bile_tax_level)

for_baiI_plot <- baiI_files$ProcessedKOBiom

## plot
baiI <- bile_plot(for_baiI_plot,
                  diet_labs,
                  diet_names_labs,
                  baiI_title)

## saving my plot outputs
ggsave(args$baiH_plot_FP,
       plot = baiH, 
       width = 12, 
       height = 4.5)

ggsave(args$baiI_plot_FP,
       plot = baiI, 
       width = 8, 
       height = 4.5)
