## 5-28-24
## new experiment metadata processing 

library(tidyverse)
library(readr)
library(magrittr)

## file paths
cdd02_meta_fp <- './new_experiments/data/misc/newExp_d15-d3_metadata.txt'
cdd03_meta_fp <- './new_experiments/data/SEQ096/new_stool_barcodes.txt'
seq_depth_fp <- './new_experiments/data/misc/newExp_d15-d3_seq_depth.tsv'

## reading in metadata file and sequencing depth
pre_cdd02_meta <- read_tsv(cdd02_meta_fp)
pre_cdd03_meta <- read_tsv(cdd03_meta_fp)
seq_depth <- read_tsv(seq_depth_fp) %>% 
  rename(`#SampleID` = sampleid)

## data wrangling/processing
## CDD02 batch (since we got the data for the batches at different times)
pre_cdd02_meta %>% 
  select(SampleID, Study, MouseSex, Diet, Vendor, ExperimentID, ExperimentDay) -> pre_cdd02_meta

pre_cdd02_meta <- pre_cdd02_meta[1:98,]

cdd02_colnames <- c('#SampleID',
                  'study',
                  'mouse_sex',
                  'diet',
                  'vendor',
                  'experiment',
                  'experiment_day')

colnames(pre_cdd02_meta) <- cdd02_colnames

pre_cdd02_meta %>% 
  mutate(separate_me = `#SampleID`,
         study = 2) %>% 
  separate_wider_delim(cols = 'separate_me',
                       delim = '.',
                       names = c('mouse1',
                                 'mouse2',
                                 'mouse3',
                                 'mouse4',
                                 'take_out'),
                       cols_remove = FALSE) %>% 
  mutate(mouse_id = paste(mouse1, mouse2, mouse3, mouse4, sep = '.')) %>% 
  select(!c('take_out', 'mouse1', 'mouse2', 'mouse3', 'mouse4', 'separate_me')) -> pre_cdd02_meta

## CDD03 batch (came later and looks slightly different)
pre_cdd03_meta %>% 
  select(`#SampleID`, MouseID, Diet, Vendor, ExperimentID, ExperimentDay) -> pre_cdd03_meta

cdd03_colnames <- c('#SampleID',
                    'mouse_id',
                    'diet',
                    'vendor',
                    'experiment',
                    'experiment_day')

colnames(pre_cdd03_meta) <- cdd03_colnames

pre_cdd03_meta %>% 
  mutate(study = 3,
         mouse_sex = 'Female') -> pre_cdd03_meta

## combining the metadata for both of the batches
comb_pre_meta <- rbind(pre_cdd02_meta,
                       pre_cdd03_meta)

## adding columns and fixing column values so they're consistent with the metadata from other experiments
comb_pre_meta %>% 
  mutate(experiment_set = 'new_exp_anschutz',
         diet = ifelse(diet == 'HFLF', 'HF/LF', diet),
         diet = ifelse(diet == 'LFLF', 'LF/LF', diet),
         diet = ifelse(diet == 'HFHF', 'HF/HF', diet),
         diet = ifelse(diet == 'LFHF', 'LF/HF', diet),
         vendor = ifelse(vendor == 'Taconic', 'taconic', 'charles_river'),
         day_post_inf = ifelse(experiment_day == 0, -15, 3),
         sample_type = 'colon') %>% 
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
    left_join(seq_depth, by = '#SampleID') -> comb_metadata


## saving processed metadata file as a .tsv
write_tsv(comb_metadata,
          './new_experiments/data/misc/proc_newExp_d15-d3_metadata.tsv')
