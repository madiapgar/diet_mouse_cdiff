## 5-30-23

## generating the sequencing depth for each sample in my total sum scaled 
## do this by grouping by sample and creating a new column that sums the total abundance per sample 
## do this using my total sum scaled file 

library(ape)
library(ggpubr)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(rstatix)

## reading in pre-total sum scaled file 
biom_fp <- '~/gut_microbiome_metabolomics/CaseyandMadi/euk_filt_mergedDietAim1table_051523-Copy1.qza'
biom <- read_qza(file = biom_fp)
biom <- biom$data 

## filtering lactococcus asvs out of my pre-total sum biom table 
## reading in lactococcus asvs and extracting them 
lacto_asv_fp <- '~/gut_microbiome_metabolomics/CaseyandMadi/lactoOnlydna-sequences.fasta'
lacto_asv <- read.FASTA(lacto_asv_fp)
lacto_asv <- names(lacto_asv)
## actual filtering steps 
biom %>% 
  as_tibble(rownames = 'asv') %>% 
  filter(!(asv %in% lacto_asv)) -> filt_biom

## doing initial total sum scaling steps 
## minus multiplying the abunds by a large number and rounding them 
filt_biom %>% 
  gather(-asv, key = sampleid, value = abund) %>% 
  group_by(sampleid) %>% 
  summarize(seq_depth = sum(abund)) -> total_sum_depth

## writing this out as a tsv 
write_tsv(total_sum_depth, '~/gut_microbiome_metabolomics/total_sum_scaled/updated/tss_seq_depth.tsv')
