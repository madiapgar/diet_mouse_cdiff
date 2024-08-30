## wrangling the metadata for the first set of experiments 
## and then combining it with the metadata I currently have for comparisons

library(tidyverse)
library(magrittr)
library(broom)

## file paths 
old_samples_fp <- './comparisons/data/misc/sample_ids.csv'
old_set4_fp <- './comparisons/data/misc/mouse_set_4_mapping.txt'
old_set3_fp <- './comparisons/data/misc/mouse_set_3_mapping.txt'
old_set5_fp <- './comparisons/data/misc/mouse_CF_set_5_mapping.txt'
new_meta_comb_fp <- './comparisons/data/misc/newExp_comp_d15_metadata.tsv'

## reading in files
old_sample_list <- read_csv(old_samples_fp)
old_set4 <- read_tsv(old_set4_fp)
old_set3 <- read_tsv(old_set3_fp)
old_set5 <- read_tsv(old_set5_fp)

new_meta_comb <- read_tsv(new_meta_comb_fp) %>% 
  select(-seq_depth) %>% 
  mutate(mouse_sex = paste('Female'))

wanted_cols <- c('#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 
                 'mouse_PID', 'abx_treatment', 'diet', 
                 'day_post_abx', 'day_post_infection', 'Description')

## first set of experiments ("old") metadata wrangling 
## fixing sample ids 
old_sample_list <- old_sample_list %>% 
                      separate_wider_delim(cols = 'sample_name',
                                           delim = '.',
                                           names = c('extra',
                                                     '#SampleID'),
                                           too_many = 'merge',
                                           cols_remove = FALSE) %>% 
                      select(-extra)

## finding which samples i need from the barcode files
## set 3
filt_old_set3 <- old_set3 %>% 
                  filter(`#SampleID` %in% old_sample_list$`#SampleID`) %>% 
  select(all_of(wanted_cols)) %>% 
  mutate(diet = ifelse(diet == 'HFD', 'HF/LF', diet),
         study = paste(3)) %>% 
  rename(mouse_id = mouse_PID) %>% 
  rename(day_post_inf = day_post_infection)

## set 5
filt_old_set5 <- old_set5 %>% 
  filter(`#SampleID` %in% old_sample_list$`#SampleID`) %>% 
  select(all_of(wanted_cols)) %>% 
  separate_wider_delim(cols = `#SampleID`,
                       delim = '.',
                       names = c('extra',
                                 'diet_actual',
                                 'others'),
                       too_many = 'merge',
                       cols_remove = FALSE) %>% 
  select(!c('extra', 'others')) %>% 
  mutate(diet = paste(diet_actual)) %>% 
  select(-diet_actual) %>% 
  mutate(diet = case_when(
          diet == 'LF' ~ 'LF/LF',
          diet == 'hfd' ~ 'HF/LF',
          diet == 'chow' ~ 'Chow'),
         study = paste(5)) %>% 
  rename(mouse_id = mouse_PID) %>% 
  rename(day_post_inf = day_post_infection)

         

correct_samples <- c(filt_old_set3$`#SampleID`,
                     filt_old_set5$`#SampleID`)

## there are none from this file in the sample list I was sent...I am confused
## it should be this one, there's just a difference in how the sample ids are written 
mini_sample_list <- old_sample_list %>% 
  filter(!(`#SampleID` %in% correct_samples))

mini_samples <- mini_sample_list %>% 
  separate_wider_delim(cols = '#SampleID',
                       delim = '.',
                       names = c('1.1',
                                 '1.2',
                                 '1.3',
                                 '2.1',
                                 '2.2',
                                 '2.3'),
                       cols_remove = TRUE) %>% 
  mutate(sample1 = paste(`1.1`, `1.2`, `1.3`, sep = '.'),
         sample2 = paste(`2.1`, `2.2`, `2.3`, sep = '.'),
         `#SampleID` = paste(sample1, sample2, sep = '_')) %>% 
  select(`#SampleID`, sample_name)

proc_sample_list <- old_sample_list %>% 
  filter(`#SampleID` %in% correct_samples)

proc_sample_list <- rbind(proc_sample_list,
                          mini_samples)

## set 4
filt_old_set4 <- old_set4 %>% 
                   filter(`#SampleID` %in% proc_sample_list$`#SampleID`) %>% 
  select(all_of(wanted_cols)) %>% 
  mutate(diet = ifelse(diet == 'chow', 'Chow', 'HF/LF'),
         study = paste(4)) %>% 
  rename(mouse_id = mouse_PID) %>% 
  rename(day_post_inf = day_post_infection)


## putting together all of the barcode files 
## associated sequencing runs:
## set 5 = SEQ024
## set 4 = SEQ021
## set 3 = SEQ016
old_metadata <- rbind(filt_old_set3,
                      filt_old_set4,
                      filt_old_set5)

baseline_old_meta <- old_metadata %>% 
                      filter(day_post_inf == -17 | day_post_inf == -15 | day_post_inf == -13)

proc_baseline_oldMeta <- baseline_old_meta %>% 
                            select(`#SampleID`, mouse_id, diet, day_post_inf, study) %>% 
                            mutate(sample_type = paste('colon'),
                                   mouse_sex = paste('Female'),
                                   experiment_set = paste('first_set_anschutz'),
                                   vendor = paste('taconic'),
                                   high_fat = ifelse(diet == 'HF/LF', 1, 0),
                                   high_fiber = paste(0),
                                   purified_diet = ifelse(diet == 'Chow', 0, 1))

## combining processed old metadata file with the new experiments combined metadata   
oldNew_meta_comb <- rbind(new_meta_comb,
                          proc_baseline_oldMeta)

## saving my outputs
write_tsv(oldNew_meta_comb,
          './comparisons/data/misc/oldNew_comp_d15_metadata.tsv')
write_tsv(filt_old_set3,
          './comparisons/data/first_set_qiime/SEQ016/oldExp_s3-016_barcodes.txt')
write_tsv(filt_old_set4,
          './comparisons/data/first_set_qiime/SEQ021/oldExp_s4-021_barcodes.txt')
write_tsv(filt_old_set5,
          './comparisons/data/first_set_qiime/SEQ024/oldExp_s5-024_barcodes.txt')
