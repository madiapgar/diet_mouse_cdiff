## 6-29-23
## creating my processed metadata file 

## needed libraries
library(ggpubr)
library(magrittr)
library(tidyverse)
library(broom)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-c",
                    "--cecal_key",
                    dest = "cecal_key_FP",
                    help = "Filepath to cecal mouse id key file in .tsv/.txt format.")
parser$add_argument("-d",
                    "--seq_depth",
                    dest = "seq_depth_FP",
                    help = "Filepath to sequencing depth file in .tsv format.")
parser$add_argument("-i",
                    "--mouse_id_facil",
                    dest = "id_facil_FP",
                    help = "Filepath to mouse ID with barrier facility file in .tsv format.")
parser$add_argument("-o",
                    "--output",
                    dest = "output_fp",
                    help = "Filepath to location for output file(s).")

args <- parser$parse_args()

## input file paths
# metadata_FP <- './cecum/data/misc/updated_cecal_metadata.tsv'
# cecal_key_FP <- '~/projects/diet_mouse_cdiff_background/cecal_key.txt'
# seq_depth_FP <- './cecum/data/misc/seq_depth.tsv'
# id_facil_FP <- './cecum/data/misc/mouseID_facil.tsv'
# output_fp <- './cecum/data/misc/cecal_processed_metadata.tsv'

## needed functions 
## 1 
## general function to prep the metadata file for further data analyses 
metadata_fixer <- function(metadata_fp,
                           cecal_key_fp) {
  metadata <- read_tsv(metadata_fp)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'
  metadata %>% 
    separate_wider_delim(cols = 'corr_sample_num',
                         delim = '.',
                         names = c('diet',
                                   'tube_num',
                                   'date'),
                         cols_remove = FALSE) %>% 
    mutate(diet = ifelse(diet == 'Ctrl', 'LF/LF', diet),
           diet = ifelse(diet == 'CtrlF', 'LF/HF', diet),
           diet = ifelse(diet == 'WD', 'HF/LF', diet),
           diet = ifelse(diet == 'WDF', 'HF/HF', diet),
           day_post_inf = 3,
           sample_type = 'cecum') %>% 
    select(diet, tube_num, date, sampleid, day_post_inf, sample_type, corr_sample_num) -> metadata
  metadata[[1]][[11]] <- 'LF/LF'
  metadata[[1]][[10]] <- 'LF/HF'
  ## cecal key
  cecal_key <- read_tsv(cecal_key_fp)
  cecal_key %>% 
    mutate(date = paste0(0, date)) -> cecal_key
  names(cecal_key)[names(cecal_key) == 'tube_numb'] <- 'tube_num'
  
  ## joining them together
  merge(metadata, cecal_key,
        by.x = c('diet',
                 'tube_num',
                 'date'),
        by.y = c('diet',
                 'tube_num',
                 'date'),
        all.x = TRUE,
        all.y = FALSE) -> actual_metadata
  return(actual_metadata)
}

## 2 
## for editing my metadata file post metadata fixer 
meta_diet_fixer <- function(metadata_file,
                            seq_depth_fp,
                            id_facil_fp){
  seq_depths <- read_tsv(seq_depth_fp)
  id_facil <- read_tsv(id_facil_fp)
  metadata_file %>% 
    mutate(high_fat = case_when(
             diet == 'HF/HF' ~ 1,
             diet == 'HF/LF' ~ 1,
             .default = 0
           ), 
           high_fiber = case_when(
             diet == 'HF/HF' ~ 1,
             diet == 'LF/HF' ~ 1,
             .default = 0
           ), 
           purified_diet = case_when(
             diet == 'Chow' ~ 0,
             .default = 1
           )) %>% 
    left_join(seq_depths) %>% 
    left_join(id_facil, by = 'mouse_id') -> metadata
  return(metadata)
}


## metadata processing and adding in the sequencing depth 
metadata_pre <- metadata_fixer(args$metadata_FP,
                               args$cecal_key_FP)

metadata <- meta_diet_fixer(metadata_pre,
                            args$seq_depth_FP,
                            args$id_facil_FP)

## writing out processed metadata file to the data/misc directory
write_tsv(metadata,
          args$output_fp)
