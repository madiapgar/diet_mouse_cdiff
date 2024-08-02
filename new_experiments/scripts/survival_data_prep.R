## 7-22-24
## putting together new experiment survival results for plots and stats

## needed libraries
library(tidyverse)
library(readr)
library(magrittr)

## file path
cfu_count_FP <- './new_experiments/data/misc/proc_combined_CDD_cfuCounts.tsv'

## reading in file 
cfuCount_table <- read_tsv(cfu_count_FP)

## data wrangling
## 0: alive (or alive until the end)
## 1: dead
surv <- cfuCount_table %>% 
              spread(key = 'location', value = 'microbe_presence') %>% 
              select(!c('spleen', 'liver', 'blood')) %>% 
              mutate(survival = ifelse(sac_exptDay == 30, 0, 1),
                     day_post_inf = case_when(
                       sac_exptDay == 17 ~ 2,
                       sac_exptDay == 18 ~ 3,
                       sac_exptDay == 19 ~ 4,
                       sac_exptDay == 20 ~ 5,
                       sac_exptDay == 21 ~ 6,
                       sac_exptDay == 22 ~ 7,
                       sac_exptDay == 23 ~ 8,
                       sac_exptDay == 29 ~ 14,
                       sac_exptDay == 30 ~ 15
                     ),
                     actual_surv = ifelse(day_post_inf == 14, 0, survival))

## mini version bc idek what I'm doing with this its humbling
write_tsv(surv,
          './new_experiments/data/misc/survival_data.tsv')
