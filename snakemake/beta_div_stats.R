## 6-26-23
## Qiime2 core diversity analysis statistical analysis for beta diversity metrics

## needed libraries
library(qiime2R)
library(tidyverse)
library(cowplot)
library(magrittr)
library(vegan)
library(viridis)
library(microshades)
library(phyloseq)
library(ggh4x)
library(broom)
library(rstatix)
library(dunn.test)

## input file paths
metadata_FP <- '../data/misc/merged_metadata1.tsv'
seq_depth_FP <- '../data/misc/tss_seq_depth.tsv'
uw_dist_fp <- '../data/qiime/core_outputs/uw_dist_matrix.tsv'
w_dist_fp <- '../data/qiime/core_outputs/w_dist_matrix.tsv'
unwanted_samples <- c('Mock20220615A', 'Mock_1A', 'Mock_2A',
                      'Mock_3A', 'Mock_4A', 'Mock_5A', 'Mock_6A',
                      'Mock_7A', 'PCR Blank0',
                      'PCR Blank1', 'Mock_7', 'Mock_6',
                      'Mock_5', 'Mock_4', 'PCR blank')

## functions in order of useage 
##1 
## initial metadata fixer 
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

## 3 
## for distance matrix processing
## for beta diversity statistical analysis 
dist_matrix_prep <- function(metadata_file,
                             dist_matrix_fp,
                             sample_filter){ 
  ## metadata filtering
  metadata_file %>% 
    filter(!(sampleid %in% sample_filter)) -> metadata
  ## distance matrix
  dist <- read_tsv(dist_matrix_fp)
  names(dist)[names(dist) == '...1'] <- 'sampleid'
  dist %>% 
    gather(-sampleid, key = sample_col, value = dist) %>% 
    filter(sampleid %in% metadata$sampleid) %>% 
    filter(sample_col %in% metadata$sampleid) %>% 
    spread(sample_col, dist) -> dist_long
  dist_long %>% 
    select(-sampleid) -> dist_proc
  metadata %>% 
    arrange(sampleid) -> metadata
  metadata %>% 
    filter(sampleid %in% dist_long$sampleid) -> filt_meta
  dist_proc <- as.matrix(dist_proc)
  row.names(dist_proc) <- colnames(dist_proc)
  filt_meta <- filt_meta[order(filt_meta$sampleid),]
  ## list of outputs
  my_list <- list(Metadata = filt_meta,
                  DistanceMatrix = dist_proc)
  return(my_list)
}


## 4 
## beta diversity adonis2 testing function
adonis_test <- function(dist_matrix,
                        metadata_file){
  adonis_results <- adonis2(as.dist(dist_matrix) ~ purified_diet * seq_depth + high_fat + high_fiber +
                              day_post_inf + study,
                            data = metadata_file,
                            permutations = 999, 
                            parallel = 4)
  adonis_results <- tidy(adonis_results)
  return(adonis_results)
}

## actually using the functions
## weighted unifrac 
stat_meta <- metadata_fixer(metadata_FP)
w_dist_files <- dist_matrix_prep(stat_meta,
                                 w_dist_fp,
                                 unwanted_samples)

w_dist <- w_dist_files$DistanceMatrix
stat_meta <- w_dist_files$Metadata


filt_stat_meta <- meta_diet_fixer(stat_meta,
                                  seq_depth_FP)

w_adonis <- adonis_test(w_dist,
                        filt_stat_meta)


## unweighted unifrac
uw_dist_files <- dist_matrix_prep(stat_meta,
                                  uw_dist_fp,
                                  unwanted_samples)

uw_dist <- uw_dist_files$DistanceMatrix
stat_meta <- uw_dist_files$Metadata


filt_stat_meta <- meta_diet_fixer(stat_meta,
                                  seq_depth_FP)

uw_adonis <- adonis_test(uw_dist,
                         filt_stat_meta)

## writing results out as a .tsv file 
write_tsv(w_adonis, '../../stats/w_adonis_results.tsv')
write_tsv(uw_adonis, '../../stats/uw_adonis_results.tsv')