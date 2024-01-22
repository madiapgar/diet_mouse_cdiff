library(tidyverse)

metadata_FP <- './comparisons/data/misc/d3_filt_comp_metadata.tsv'
survival_FP <- './colon/data/misc/aim1a_survival.csv'

metadata <- read_tsv(metadata_FP)
survival <- read_csv(survival_FP)

## all mice, without study one filtered out 
## day 3 is also filtered out of the stool metadata 
## should I reflect the number of mice used for the analysis or the actual total number of mice?
metadata %>% 
  filter(day_post_inf == -15) %>% 
  count(diet) -> n_sample_stool

metadata %>% 
  filter(sample_type == 'cecum') %>% 
  count(diet) -> n_sample_cecum


metadata %>% 
  ungroup() %>% 
  distinct(mouse_id, .keep_all = TRUE) %>% 
  count(diet) -> n_mouse_day3


## survival cohort
survival %>% 
  count(diet) -> n_mouse_survival 
