## 5-28-24
## new experiment metadata processing 

library(tidyverse)
library(readr)
library(magrittr)

## file paths
meta_fp <- './new_experiments/data/misc/newExp_d15-d3_metadata.txt'
seq_depth_fp <- './new_experiments/data/misc/newExp_d15-d3_seq_depth.tsv'

## reading in metadata file and sequencing depth
pre_meta <- read_tsv(meta_fp)
seq_depth <- read_tsv(seq_depth_fp)

## data wrangling/processing
pre_meta %>% 
  select(SampleID, Study, MouseSex, Diet, Vendor, ExperimentID, ExperimentDay) -> pre_meta

pre_meta <- pre_meta[1:98,]

new_colnames <- c('sampleid',
                  'study',
                  'mouse_sex',
                  'diet',
                  'vendor',
                  'experiment',
                  'experiment_day')

colnames(pre_meta) <- new_colnames

pre_meta %>% 
  mutate(separate_me = sampleid,
         experiment_set = 'new_exp_anschutz',
         mouse_sex = 'F',
         study = 2,
         diet = ifelse(diet == 'HFLF', 'HF/LF', diet),
         diet = ifelse(diet == 'LFLF', 'LF/LF', diet),
         diet = ifelse(diet == 'HFHF', 'HF/HF', diet),
         diet = ifelse(diet == 'LFHF', 'LF/HF', diet),
         vendor = ifelse(vendor == 'Taconic', 'taconic', 'charles_river'),
         day_post_inf = ifelse(experiment_day == 0, -15, 3),
         sample_type = 'colon') %>% 
  separate_wider_delim(cols = 'separate_me',
                       delim = '.',
                       names = c('take_out1',
                                 'mouse1',
                                 'mouse2',
                                 'mouse3',
                                 'take_out2'),
                       cols_remove = FALSE) %>% 
  mutate(mouse_id = paste(mouse1, mouse2, mouse3, sep = '.')) %>% 
  select(!c('take_out1', 'take_out2', 'mouse1', 'mouse2', 'mouse3', 'separate_me')) %>% 
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
  left_join(seq_depth,
            by = 'sampleid') -> metadata

## saving processed metadata file as a .tsv
write_tsv(metadata,
          './new_experiments/data/misc/proc_newExp_d15-d3_metadata.tsv')
