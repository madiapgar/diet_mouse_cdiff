## 6-29-23
## calculating the sequencing depth for each sample for later statistical use

## needed libraries 
library(ape)
library(ggpubr)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(rstatix)

## reading in pre-total sum scaled file 
biom_fp <- './data/misc/euk_filt_mergedDietAim1table_051523-Copy1.qza'
biom <- read_qza(file = biom_fp)
biom <- biom$data 

## filtering lactococcus asvs out of my pre-total sum biom table 
## reading in lactococcus asvs and extracting them 
lacto_asv_fp <- './data/misc/lactoOnlydna-sequences.fasta'
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
write_tsv(total_sum_depth, 
          './data/misc/tss_seq_depth.tsv')