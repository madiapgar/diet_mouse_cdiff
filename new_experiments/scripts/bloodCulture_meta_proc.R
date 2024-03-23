## 3-13-24
## this script processes the blood culture metadata file

## needed libraries
library(tidyverse)
library(magrittr)
library(broom)

## input file path
pre_meta_FP <- './new_experiments/data/misc/blood_culture_metadata.tsv'
pre_colony_FP <- './new_experiments/data/misc/blood_colony_count.txt'

## reading in metadata and colony count files
pre_meta <- read_tsv(pre_meta_FP)
pre_colony <- read_tsv(pre_colony_FP)

## can't combine colony count and metadata files since the sampleids are different - not sure how to 
## approach that problem right now

## processing colony count file 
pre_colony %>% 
  gather(-sampleid, key = 'separate_me', value = 'colony_count') %>% 
  separate_wider_delim(cols = 'separate_me',
                       delim = '.',
                       names = c('media',
                                 'oxygen_tol',
                                 'cultured_from'),
                       cols_remove = FALSE) %>% 
  select(-separate_me) %>% 
  mutate(another_sep = sampleid) %>% 
  separate_wider_delim(cols = 'another_sep',
                       delim = '.',
                       names = c('experiment',
                                 'vendor',
                                 'diet',
                                 'num'),
                       cols_remove = FALSE) %>% 
  select(!c('num', 'another_sep')) %>% 
  mutate(diet = ifelse(diet == 'LFLF', 'LF/LF', diet),
         vendor = ifelse(vendor == 'CR', 'charles_river', 'taconic')) -> colony_count

## processing metadata file 
names(pre_meta)[names(pre_meta) == '#SampleID'] <- 'sampleid'

pre_meta %>% 
  select(sampleid) %>% 
  mutate(separate_me = sampleid) %>% 
  separate_wider_delim(cols = 'separate_me',
                       delim = '.',
                       names = c('experiment',
                                 'vendor',
                                 'diet',
                                 'num',
                                 'cultured_from'),
                       cols_remove = FALSE) %>% 
  select(-separate_me) %>% 
  mutate(vendor = ifelse(vendor == 'CR', 'charles_river', 'taconic'),
         diet = ifelse(diet == 'LFLF', 'LF/LF', diet),
         cultured_from = ifelse(cultured_from == 'aer', 'aerobic', cultured_from),
         cultured_from = ifelse(cultured_from == 'an', 'anaerobic', cultured_from),
         oxygen_tol = cultured_from,
         cultured_from = ifelse(cultured_from == 'aerobic' | cultured_from == 'anaerobic', 
                                'blood', cultured_from),
         oxygen_tol = ifelse(oxygen_tol == 'spleen', 'idk', oxygen_tol),
         sample_type = 'culture') -> meta

## writing output file as a .tsv
write_tsv(meta,
          './new_experiments/data/misc/proc_blood_culture_meta.tsv')
write_tsv(colony_count,
          './new_experiments/data/misc/proc_colony_count.tsv')

