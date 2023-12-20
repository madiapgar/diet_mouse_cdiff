## 12-20-23
## creating my processed bile acid file

## needed libraries
library(ggpubr)
library(magrittr)
library(tidyverse)
library(broom)

## input file paths 
bile_acid_FP <- '../data/misc/bile_acid.txt'
mouseID_facil_FP <- '../data/misc/mouseID_facil.tsv'
cecal_key_FP <- '../data/misc/cecal_key.txt'

## needed functions
## 1
bile_acid_processing <- function(bile_acid_fp,
                                 mouseID_facil_fp,
                                 cecal_key_fp){
  ## reading in bile acid file
  bile_acid <- read_tsv(bile_acid_fp)
  ## reading in cecal key file
  cecal_key <- read_tsv(cecal_key_fp)
  ## reading in mouse ID with facility file
  mouseID_facil <- read_tsv(mouseID_facil_fp)
  ## prepping bile acid table for joining with other data
  bile_acid %>% 
    separate_wider_delim(cols = 'tube_label',
                         delim = '.',
                         names = c('diet',
                                   'tube_numb'),
                         cols_remove = FALSE) %>% 
    mutate(diet = ifelse(diet == 'Ctrl', 'LF/LF', diet),
           diet = ifelse(diet == 'Ctrl+F', 'LF/HF', diet),
           diet = ifelse(diet == 'WD', 'HF/LF', diet),
           diet = ifelse(diet == 'WD+F', 'HF/HF', diet)) %>% 
    select(!c('reisdorph_lab_id', 'sample_type')) -> bile_acid
  
  names(bile_acid)[names(bile_acid) == 'collection_date'] <- 'date'
  ## correcting tubes that were mislabeled under 50622
  bile_acid[[1]][[20]] <- 'LF/LF'
  bile_acid[[1]][[21]] <- 'LF/LF'
  bile_acid[[1]][[22]] <- 'HF/LF'
  ## putting bile acid table and cecal key together
  merge(bile_acid, cecal_key,
        by.x = c('diet',
                 'tube_numb',
                 'date'),
        by.y = c('diet',
                 'tube_numb',
                 'date'),
        all.x = TRUE,
        all.y = FALSE) -> actual_bile_acid
  ## putting mouse ID with facility information on the table
  actual_bile_acid %>% 
    left_join(mouseID_facil, by = 'mouse_id') -> actual_bile_acid
  return(actual_bile_acid)
}

## using the function 
bile_acid <- bile_acid_processing(bile_acid_FP,
                                  mouseID_facil_FP,
                                  cecal_key_FP)

## writing results out as a .tsv
write_tsv(bile_acid,
          '../data/misc/corrected_bile_acid.tsv')