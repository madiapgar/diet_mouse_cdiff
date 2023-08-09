## 7-14-23
## this script contains the stats for the wanted microbe family abundances 
## by days relative to infection and diet

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(vegan)
library(viridis)

## input file paths and others
otu_table_FP <- './data/qiime/taxonomy_filtered.qza'
tax_FP <- './data/qiime/taxonomy.qza'
metadata_FP <- './data/misc/processed_metadata.tsv'
wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Tannerellaceae', 'Muribaculaceae', 'Bacteroidaceae', 
                   'Marinifilaceae')

## function
family_abun_file_prep <- function(metadata_fp,
                                  tax_fp,
                                  otu_table_fp,
                                  tax_level,
                                  wanted_tax){
  ## metadata
  metadata <- read_tsv(metadata_FP)
  ## taxonomy
  taxonomy <- read_qza(tax_FP)$data %>% 
    parse_taxonomy() %>% 
    rownames_to_column('asv')
  ## otu table 
  otu_table <- read_qza(otu_table_FP)$data
  otu_table %>% 
    as_tibble(rownames = 'asv') %>% 
    gather(-asv, key = sampleid, value = abun) %>% 
    group_by(sampleid) %>% 
    mutate(rel_abun = abun/sum(abun)) %>% 
    mutate(p_abun = rel_abun + 0.000001) -> otu_table
  ## joining all tables together 
  otu_table %>% 
    left_join(metadata, by = 'sampleid') %>% 
    left_join(taxonomy, by = 'asv') -> abun_table
  abun_table %>% 
    group_by(sampleid, day_post_inf, diet, mouse_id, 
             purified_diet, high_fat, high_fiber, 
             seq_depth, .data[[tax_level]]) %>% 
    summarise(rel_abund = sum(rel_abun)) %>% 
    filter(.data[[tax_level]] %in% wanted_tax) %>% 
    mutate(mouse_fact = as.factor(mouse_id),
           day_fact = as.factor(day_post_inf)) -> abun_filt
  ## creating a list for my outputs
  my_list <- list(Metadata = metadata,
                  Taxonomy = taxonomy,
                  OTUTable = otu_table,
                  AbundanceTable = abun_filt)
  return(my_list)
}

## prepping the file needed to run the linear model
abun_files <- family_abun_file_prep(metadata_FP,
                                    tax_FP,
                                    otu_table_FP,
                                    wanted_level,
                                    wanted_family)

## can still pull out the metadata, taxonomy, and otu table as well
abun_filt <- abun_files$AbundanceTable

## performing the linear model
## looking at each microbe family by days relative to infection and the different diets
## takes out the intercept term and only shows significant p values 
abun_filt %>%
  na.omit() %>% 
  group_by(Family, day_post_inf) %>% 
  do(tidy(lm(rel_abund ~ (purified_diet * seq_depth) + high_fat + high_fiber,
             data =.))) %>% 
  adjust_pvalue(method = 'BH') %>% 
  na.omit() %>% 
  filter(term != '(Intercept)',
         p.value <= 0.05) -> family_abun_lm

## saving my output as a .tsv
write_tsv(family_abun_lm,
          './stats/family_abun_lm.tsv')
