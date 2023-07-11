## 6-26-23
## Qiime2 core diversity analysis statistical analysis for alpha diversity metrics

## needed libraries
library(qiime2R)
library(tidyverse)
library(cowplot)
library(magrittr)
library(vegan)
library(ggh4x)
library(broom)
library(rstatix)
library(dunn.test)

## input file paths
metadata_FP <- './data/misc/processed_metadata.tsv'
faith_pd_fp <- './data/qiime/core_outputs/faith_pd.tsv'
shannon_fp <- './data/qiime/core_outputs/shannon_entropy.tsv'
unwanted_samples <- c('Mock20220615A', 'Mock_1A', 'Mock_2A',
                      'Mock_3A', 'Mock_4A', 'Mock_5A', 'Mock_6A',
                      'Mock_7A', 'PCR Blank0',
                      'PCR Blank1', 'Mock_7', 'Mock_6',
                      'Mock_5', 'Mock_4', 'PCR blank')

## functions in order that they're used
## 1
## alpha diversity file prep 
alpha_div_prep <- function(file_path1,
                           file_path2,
                           sample_filter,
                           metadata_fp){
  ## faith's pd 
  alpha_faith <- read_tsv(file_path1)
  names(alpha_faith)[names(alpha_faith) == '#SampleID'] <- 'sampleid'
  alpha_faith %>% 
    filter(!(sampleid %in% sample_filter)) -> faith_pd
  ## metadata file for both
  stat_meta <- read_tsv(metadata_fp)
  stat_meta %>% 
    filter(!(sampleid %in% sample_filter)) -> stat_meta
  ## joining faith's pd and metadata file together into one table
  stat_meta %>% 
    filter(sampleid %in% faith_pd$sampleid) %>% 
    left_join(faith_pd, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> faith_biom
  ## shannon entropy
  alpha_shannon <- read_tsv(file_path2)
  names(alpha_shannon)[names(alpha_shannon) == '...1'] <- 'sampleid'
  alpha_shannon %>% 
    filter(!(sampleid %in% sample_filter)) -> shannon
  ## joining shannon and metadata file together into one table 
  stat_meta %>% 
    filter(sampleid %in% shannon$sampleid) %>% 
    left_join(shannon, by = 'sampleid') %>% 
    filter(!is.na(diet)) -> shannon_biom
  ## creating a list for outputs 
  my_list <- list(FaithPD = faith_biom,
                  Shannon = shannon_biom, 
                  Metadata = stat_meta)
  return(my_list)
}

## 2
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

## 3
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
                              metadata_FP)

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
write_tsv(faith_lm, './stats/faith_total_results.tsv')
write_tsv(sectioned_faith_lm, './stats/faith_diet_results.tsv')
write_tsv(shannon_lm, './stats/shannon_total_results.tsv')
write_tsv(sectioned_shannon_lm, './stats/shannon_diet_results.tsv')