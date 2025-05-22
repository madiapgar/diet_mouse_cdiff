## new experiment CDD02/CDD03 culture samples metadata processing
## 9-27-24

## needed libraries
library(tidyverse)
library(magrittr)
library(broom)

## file paths
metadata_FP <- './new_experiments/data/new_culture_seqs/culture_mock_barcodes.txt'
colony_FP <- './new_experiments/data/misc/newExp_cfus.txt'

## metadata pre-processing
mock_list <- c('EW.Mock.G20',
               'EW.Mock.H20',
               'EW.Mock.I20',
               'EW.Mock.A20',
               'EW.Mock.B20',
               'EW.Mock.C20',
               'EW.Mock.D20',
               'EW.Mock.E20',
               'EW.Mock.F20')
new_col_names <- c('sampleid',
                   'sample_type',
                   'mouse_id',
                   'diet',
                   'vendor',
                   'experiment')

## filtering out extraction mock samples
pre_meta %>% 
  select(`#SampleID`, sampleType, MouseID, Diet, Vendor, ExperimentID) %>% 
  filter(!(`#SampleID` %in% mock_list)) -> pre_meta

colnames(pre_meta) <- new_col_names

pre_meta %>%
  mutate(diet = case_when(
    diet == 'Chow' ~ 'Chow',
    diet == 'HFLF' ~ 'HF/LF',
    diet == 'HFHF' ~ 'HF/HF',
    diet == 'LFHF' ~ 'LF/HF',
    diet == 'LFLF' ~ 'LF/LF'
  ),
  vendor = case_when(
    vendor == 'CharlesRiver' ~ 'charles_river',
    vendor == 'Taconic' ~ 'taconic'
  ),
  separate_me = paste(sampleid)) %>% 
  separate_wider_delim(cols = 'separate_me',
                       delim = '.',
                       names = c('extra1',
                                 'extra2',
                                 'extra3',
                                 'extra4',
                                 'location'),
                       cols_remove = FALSE) %>% 
  select(!c('extra1', 'extra2', 'extra3', 'extra4', 'separate_me')) %>% 
  mutate(location = tolower(location)) -> unfixed_meta

## fixing a sample where the experiment, vendor, and diet were incorrect in the metadata
## and rejoining it with the metadata file
unfixed_meta %>% 
  filter(sampleid == "CDD02.CR.LFHF.2.Liver") %>% 
  mutate(mouse_id = paste('CDD02.CR.LFHF.2'),
         diet = paste('LF/HF'),
         vendor = paste('charles_river'),
         experiment = paste('CDD02')) -> weird_sample

unfixed_meta %>% 
  filter(sampleid != "CDD02.CR.LFHF.2.Liver") -> int_meta

proc_meta <- rbind(int_meta,
                   weird_sample)

## cfu file pre-processing to be joined with the metadata
pre_colony <- read_tsv(colony_FP)

pre_colony %>% 
  rename('sampleid' = 'mouse_id.tissue') %>% 
  mutate(colony_count = ifelse(colony_count == 'BLD', 0, colony_count),
         colony_count = ifelse(colony_count == 'TNTC', 2000, colony_count),
         colony_count = as.numeric(colony_count)) %>% 
  select(sampleid, colony_count) -> proc_colony

## combining the metadata and cfu tables 
proc_meta %>% 
  left_join(proc_colony, by = 'sampleid') -> meta_cfu_table

## saving the processed metadata file
write_tsv(meta_cfu_table,
          './new_experiments/data/misc/proc_cdd02-3_culture_metadata.tsv')
