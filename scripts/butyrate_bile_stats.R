## 7-14-23
## creating a script that runs linear modeling statistics on the picrust 
## taxon functional abundance outputs for butyrate and secondary bile acid enzymes

## needed libraries 
packages <- c("ggpubr", 
              "magrittr", 
              "qiime2R", 
              "tidyverse", 
              "broom",
              "cowplot",
              "viridis",
              "ape",
              "AICcmodavg",
              "rstatix")

for(package in packages){
  if(!require(package, character.only = T)){
    install.packages(package)
  }
  library(package)
}

## input file paths and KOs
metadata_FP <- './data/misc/processed_metadata.tsv'
ko_contrib_FP <- './data/picrust/tss3_meta_contrib.tsv'
but_kos <- c('K00929','K01034')
bile_kos <- c('K15873', 'K15874')

## needed function
stat_file_prep <- function(metadata_fp,
                           ko_contrib_fp,
                           ko_list){
  ## metadata
  metadata <- read_tsv(metadata_fp)
  names(metadata)[names(metadata) == 'sampleid'] <- 'sample'
  ## ko meta contrib 
  ko_contrib <- read_tsv(ko_contrib_fp)
  ko_contrib %>% 
    left_join(metadata, by = 'sample') -> stat_biom
  ## messing with biom table format so that the zeroes are represented
  stat_biom %>% 
    filter(ko %in% ko_list) %>% 
    group_by(ko, sample, diet, day_post_inf, purified_diet, high_fat, high_fiber, mouse_id, seq_depth) %>% 
    summarise(taxon_function_abun = sum(taxon_function_abun)) %>% 
    filter(!is.na(day_post_inf)) %>% 
    spread(day_post_inf, taxon_function_abun, fill = 0) %>% 
    gather(-ko, -sample, -diet, -purified_diet, -high_fat, -high_fiber, -mouse_id, -seq_depth,
           key = day_post_inf, value = taxon_function_abun) -> biom_long
  return(biom_long)
}

## file prep 
## butyrate 
but_long <- stat_file_prep(metadata_FP,
                           ko_contrib_FP,
                           but_kos)

## bile acids
bile_long <- stat_file_prep(metadata_FP,
                            ko_contrib_FP,
                            bile_kos)

## butyrate linear model
but_long %>% 
  group_by(ko, day_post_inf) %>% 
  do(tidy(lm(taxon_function_abun ~ high_fat + high_fiber + (purified_diet * seq_depth), 
             data = .))) %>% 
  adjust_pvalue(method = 'BH') %>% 
  filter(p.value <= 0.05) -> buty_lm

## bile acid linear model
bile_long %>% 
  group_by(ko, day_post_inf) %>% 
  do(tidy(lm(taxon_function_abun ~ high_fat + high_fiber + (purified_diet * seq_depth), 
             data = .))) %>% 
  adjust_pvalue(method = 'BH') %>% 
  filter(p.value <= 0.05) -> bile_lm

## saving my outputs as a .tsv
write_tsv(buty_lm,
          './stats/buty_enzyme_lm.tsv')
write_tsv(bile_lm,
          './stats/bile_enzyme_lm.tsv')