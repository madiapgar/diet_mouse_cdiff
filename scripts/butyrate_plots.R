## 7-13-23
## script that contains all of the needed functions for butyrate kinase and
## butyryl coa transferase plots 

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
parser$add_argument("-bk",
                    "--buk_plot",
                    dest = "buk_plot_FP",
                    help = "Filepath to butyrate kinase plot in .pdf format.")
parser$add_argument("-bt",
                    "--but_plot",
                    dest = "but_plot_FP",
                    help = "Filepath to butyryl-coa transferase plot in .pdf format.")

args <- parser$parse_args()

## input file paths 
# metadata_FP <- './data/misc/processed_metadata.tsv'
# tax_FP <- './data/qiime/taxonomy.qza'
# ko_contrib_FP <- './data/picrust/tss3_meta_contrib.tsv'

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
## butyrate kinase specs
buk_tax_level <- 'class'
buk_thresh_level <- 200
buk_ko <- 'K00929'
buk_title <- 'Butyrate Kinase Potential Over Time'
## butyryl coa transferase specs
but_tax_level <- 'class'
but_thresh_level <- 50
but_ko <- 'K01034'
but_title <- 'Butyryl-CoA Transferase Potential Over Time'

## functions in order of usage 
## 1 
## input file prep (metadata, taxonomy, ko meta contrib prep) for butyrate plots
buty_file_prep <- function(tax_fp,
                           ko_contrib_fp,
                           metadata_fp,
                           ko_list,
                           taxonomy_level,
                           threshold){
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
    select(sample, ko, taxon_function_abun, study, diet, day_post_inf, 
           any_of(taxonomy_level)) %>% 
    filter(ko %in% ko_list) -> filtered_biom
  ## summing the wanted taxonomic level's abundance for a particular ko
  filtered_biom %>% 
    group_by(sample, ko, study, diet, day_post_inf, .data[[taxonomy_level]]) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    ungroup() -> filtered_biom_sum
  ## creating a facet that will contain the total plot to reference
  ## fill is set at 0.01 since its less than the lowest abundance which is 0.88
  THRESHOLD = threshold
  filtered_biom_sum %>% 
    group_by(sample, day_post_inf, ko, diet, study) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    mutate(class = 'Total') %>% 
    filter(!is.na(day_post_inf)) %>%
    spread(day_post_inf, taxon_function_abun, fill = 0.01) %>% 
    gather(-sample, -ko, -diet, -study, -class,
           key = day_post_inf, value = taxon_function_abun) %>% 
    mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
    group_by(class) %>% 
    mutate(other_col = mean(taxon_function_abun),
           class = if_else(other_col < THRESHOLD, 'Other', class)) %>% 
    arrange(other_col) -> filtered_sample_abun
  ## classifying all taxonomic level abundances below determined threshold as 'Other'
  ## needs to be processed the exact same way as above so rbind() will work
  filtered_biom_sum %>% 
    filter(!is.na(day_post_inf)) %>%
    spread(day_post_inf, taxon_function_abun, fill = 0.01) %>% 
    gather(-sample, -ko, -diet, -study, -class,
           key = day_post_inf, value = taxon_function_abun) %>% 
    mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
    group_by(class) %>% 
    mutate(other_col = mean(taxon_function_abun),
           class = if_else(other_col < THRESHOLD, 'Other', class)) %>% 
    arrange(other_col) -> filtered_sum_other
  ## rbinding filtered_sample_abun and filtered_sum_other together into the same table
  filtered_sum_other %>% 
    rbind(filtered_sample_abun) -> big_filtered_sum
  ## creating a list of my outputs
  my_list <- list(Taxonomy = taxonomy,
                  Metadata = metadata,
                  KOContrib = kometa_contrib,
                  ProcessedKOBiom = big_filtered_sum)
  return(my_list)
}

## 2 
## buty plot construction 
butyrate_plot <- function(processed_ko_biom,
                          labels,
                          names_labels,
                          title){
  labs <- labels
  names(labs) <- names_labels
  processed_ko_biom %>% 
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = taxon_function_abun)) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_violin(aes(group = day_post_inf), outlier.shape = NA) +
    geom_smooth(se = FALSE, size = 0.5) +
    geom_jitter(width = 0.1, height = 0, 
                alpha = 0.4) +
    geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
    facet_grid(class~diet, labeller = labeller(diet = labs)) +
    theme_bw(base_size = 14) +
    theme(legend.text = element_text(size = 8.5),
          strip.text.y = element_text(angle = 0)) +
    guides(color = guide_legend(override.aes = list(size = 0.9))) +
    ggtitle(title) +
    xlab('Days Relative to Infection') +
    ylab('KO Counts') -> butyrate
  return(butyrate)
}

## butyrate kinase
## file prep
buk_files <- buty_file_prep(args$tax_FP,
                            args$ko_contrib_FP,
                            args$metadata_FP,
                            buk_ko,
                            buk_tax_level,
                            buk_thresh_level)
## metadata, taxonomy, and kometa_contrib are universal for all plots 
metadata <- buk_files$Metadata
taxonomy <- buk_files$Taxonomy
kometa_contrib <- buk_files$KOContrib
for_buk_plot <- buk_files$ProcessedKOBiom

## plot 
butyrate_kinase <- butyrate_plot(for_buk_plot,
                                 diet_labs,
                                 diet_names_labs,
                                 buk_title)

## butyryl coa transferase 
## file prep
but_files <- buty_file_prep(args$tax_FP,
                            args$ko_contrib_FP,
                            args$metadata_FP,
                            but_ko,
                            but_tax_level,
                            but_thresh_level)

for_but_plot <- but_files$ProcessedKOBiom

## plot
butyryl_coa_transferase <- butyrate_plot(for_but_plot,
                                         diet_labs,
                                         diet_names_labs,
                                         but_title)

## saving my plot outputs
ggsave(args$buk_plot_FP, 
       plot = butyrate_kinase,
       width = 11, 
       height = 7)

ggsave(args$but_plot_FP,
       plot = butyryl_coa_transferase, 
       width = 11, 
       height = 7)