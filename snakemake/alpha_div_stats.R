## 6-26-23
## Qiime2 core diversity analysis statistical analysis for alpha diversity metrics

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
metadata_FP <- '../../data/misc/merged_metadata1.tsv'
seq_depth_FP <- '../../data/misc/tss_seq_depth.tsv'
faith_pd_fp <- '../../data/qiime/core_outputs/faith_pd.tsv'
shannon_fp <- '../../data/qiime/core_outputs/shannon_entropy.tsv'
unwanted_samples <- c('Mock20220615A', 'Mock_1A', 'Mock_2A',
                      'Mock_3A', 'Mock_4A', 'Mock_5A', 'Mock_6A',
                      'Mock_7A', 'PCR Blank0',
                      'PCR Blank1', 'Mock_7', 'Mock_6',
                      'Mock_5', 'Mock_4', 'PCR blank')

## functions in order that they're used
## 1 
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
## alpha diversity file prep 
alpha_div_prep <- function(file_path1,
                           file_path2,
                           sample_filter,
                           metadata_fp,
                           seq_depth_fp){
  ## faith's pd 
  alpha_faith <- read_tsv(file_path1)
  names(alpha_faith)[names(alpha_faith) == '#SampleID'] <- 'sampleid'
  alpha_faith %>% 
    filter(!(sampleid %in% sample_filter)) -> faith_pd
  ## metadata file for both
  stat_meta <- metadata_fixer(metadata_fp)
  stat_meta %>% 
    filter(!(sampleid %in% sample_filter)) -> stat_meta
  ## joining faith's pd and metadata file together into one table
  faith_stat_meta <- meta_diet_fixer(stat_meta,
                                     seq_depth_fp)
  faith_stat_meta %>% 
    filter(sampleid %in% faith_pd$sampleid) %>% 
    left_join(faith_pd, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> faith_biom
  ## shannon entropy
  alpha_shannon <- read_tsv(file_path2)
  names(alpha_shannon)[names(alpha_shannon) == '...1'] <- 'sampleid'
  alpha_shannon %>% 
    filter(!(sampleid %in% sample_filter)) -> shannon
  ## joining shannon and metadata file together into one table 
  shannon_stat_meta <- meta_diet_fixer(stat_meta,
                                       seq_depth_fp)
  shannon_stat_meta %>% 
    filter(sampleid %in% shannon$sampleid) %>% 
    left_join(shannon, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> shannon_biom
  ## creating a list for outputs 
  my_list <- list(FaithPD = faith_biom,
                  Shannon = shannon_biom, 
                  Metadata = stat_meta)
  return(my_list)
}

## 4 
## stats calculations
## faith's pd 
faith_div_stats <- function(biom_table){
  ## alpha_cat is what the alpha div column is called (faith_pd or shannon_entropy)
  ## sectioned out by diet 
  biom_table %>% 
    group_by(day_post_inf) %>% 
    do(tidy(lm(faith_pd ~ (purified_diet * seq_depth) + high_fat + high_fiber + study,
               data = .))) -> sectioned_lm
  sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> sectioned_lm
  ## not sectioned out by diet 
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(faith_pd ~ diet * seq_depth,
               data = .))) -> not_sectioned_lm
  not_sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> not_sectioned_lm
  ## creating a list 
  my_list <- list(DietSpecific = sectioned_lm,
                  OverallDiet = not_sectioned_lm)
  return(my_list)
}

## 5 
## shannon entropy 
shannon_div_stats <- function(biom_table){
  ## alpha_cat is what the alpha div column is called (faith_pd or shannon_entropy)
  ## sectioned out by diet 
  biom_table %>% 
    group_by(day_post_inf) %>% 
    do(tidy(lm(shannon_entropy ~ (purified_diet * seq_depth) + high_fat + high_fiber + study,
               data = .))) -> sectioned_lm
  sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> sectioned_lm
  ## not sectioned out by diet 
  biom_table %>%
    group_by(day_post_inf) %>% 
    do(tidy(lm(shannon_entropy ~ diet * seq_depth,
               data = .))) -> not_sectioned_lm
  not_sectioned_lm %>% 
    filter(day_post_inf != -15) %>% 
    filter(p.value <= 0.05) -> not_sectioned_lm
  ## creating a list 
  my_list <- list(DietSpecific = sectioned_lm,
                  OverallDiet = not_sectioned_lm)
  return(my_list)
}

## use of functions 
## alpha diversity analysis  
alpha_files <- alpha_div_prep(faith_pd_fp,
                              shannon_fp,
                              unwanted_samples,
                              metadata_FP,
                              seq_depth_FP)

faith <- alpha_files$FaithPD
shannon <- alpha_files$Shannon
metadata <- alpha_files$Metadata

## faith's pd stats
faith_stats <- faith_div_stats(faith)
sectioned_faith_lm <- faith_stats$DietSpecific
faith_lm <- faith_stats$OverallDiet

## shannon entropy stats
shannon_stats <- shannon_div_stats(shannon)
sectioned_shannon_lm <- shannon_stats$DietSpecific
shannon_lm <- shannon_stats$OverallDiet

## writing out results as a .tsv file 
write_tsv(faith_lm, '../../stats/faith_total_results.tsv')
write_tsv(sectioned_faith_lm, '../../stats/faith_diet_results.tsv')
write_tsv(shannon_lm, '../../stats/shannon_total_results.tsv')
write_tsv(sectioned_shannon_lm, '../../stats/shannon_diet_results.tsv')