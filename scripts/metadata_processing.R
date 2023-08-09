## 6-29-23
## creating my processed metadata file 

## needed libraries
library(ggpubr)
library(magrittr)
library(tidyverse)
library(broom)

## input file paths
metadata_FP <- './data/misc/merged_metadata1.tsv'
seq_depth_FP <- './data/misc/tss_seq_depth.tsv'

## needed functions 
## 1 
## general function to prep the metadata file for further data analyses 
metadata_fixer <- function(metadata_fp) {
  tmpMeta <- read_tsv(metadata_fp, n_max = 2)
  mycols <- colnames(tmpMeta)
  metadata <- read_tsv(metadata_fp, skip = 2, col_names = mycols)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'
  metadata %>% 
    filter(!is.na(diet)) %>% 
    mutate(day_post_inf = if_else(day_post_inf == 2, 3, day_post_inf)) %>% 
    mutate(diet = as.factor(diet)) -> metadata
  return(metadata)
}

## 2 
## for editing my metadata file post metadata fixer 
meta_diet_fixer <- function(metadata_file,
                            seq_depth_fp){
  seq_depths <- read_tsv(seq_depth_fp)
  metadata_file %>% 
    select(sampleid, diet, day_post_inf, mouse_id, study) %>% 
    mutate(diet_true = diet,
           diet_true = if_else(day_post_inf == -15, "Chow", diet_true),
           high_fat = case_when(
             diet_true == 'HF/HF' ~ 1,
             diet_true == 'HF/LF' ~ 1,
             .default = 0
           ), 
           high_fiber = case_when(
             diet_true == 'HF/HF' ~ 1,
             diet_true == 'LF/HF' ~ 1,
             diet_true == 'Chow' ~ 1,
             .default = 0
           ), 
           purified_diet = case_when(
             diet_true == 'Chow' ~ 0,
             .default = 1
           )
    ) %>% 
    left_join(seq_depths) -> metadata
  return(metadata)
}


## metadata processing and adding in the sequencing depth 
metadata_pre <- metadata_fixer(metadata_FP)
metadata <- meta_diet_fixer(metadata_pre,
                            seq_depth_FP)

## writing out processed metadata file to the data/misc directory
write_tsv(metadata,
          './data/misc/processed_metadata.tsv')
