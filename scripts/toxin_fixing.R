## attempting to find a faster way to match mouse_ids with the numbering system on the toxin data

# change to scripts directory if not there already
curr_dir <- getwd()
curr_dir <- str_split(curr_dir, '\\/')
if (curr_dir[length(curr_dir)] != 'scripts'){
  setwd('./scripts')
}

## needed libraries 
library(broom)
library(cowplot)
library(magrittr)
library(qiime2R)
library(tidyverse)

## file paths 
metab_FP <- '../data/misc/metabolomics.csv'
toxin_FP <- '../data/misc/toxin.csv'

## reading in metabolomics and toxin .csvs
metab <- read_csv(metab_FP)

metab %>% 
  select('Tube_Label', 'Collection Date', 'Sample_Type', 'mouse_id') -> metab_mini



toxin_pre <- read_csv(toxin_FP)

toxin_pre %>% 
  filter(Tube_Label != 'Non-standard dilution factor') %>% 
  filter(!is.na(Tube_Label)) -> toxin_pre

## join_by allows me to match up the tube label names with the dates they belong to 
## since there are repeating tube label names that are dependent on collection date 
toxin_pre %>%
  left_join(metab_mini, join_by('Tube_Label', 'Collection Date')) -> toxin

toxin %>% 
  select('Tube_Label', 'Collection Date', 
         'Sample_Type.x', 'Total TcA Neat', 
         'Total TcB Neat', 'Total TcA 1:10',
         'Total TcB 1:10', 'Extra_Sample', 
         'mouse_id') -> toxin_final

names(toxin_final)[names(toxin_final) == 'Sample_Type.x'] <- 'Sample_Type'

write_tsv(toxin_final,
          '../data/misc/toxin_final_data.tsv')
