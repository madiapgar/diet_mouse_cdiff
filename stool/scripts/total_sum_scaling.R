## 6-29-23
## total sum scaling of 16S rDNA data 
## since we can't rarefy it due to lactococcus contamination

## needed libraries 
library(ape)
library(ggpubr)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(rstatix)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-b",
                    "--biom",
                    dest = "biom_fp",
                    help = "Filepath to biom table in .qza format.")
parser$add_argument("-s",
                    "--sequence",
                    dest = "lacto_asv_fp",
                    help = "Filepath to wanted sequences in .fasta format.")
parser$add_argument("-o",
                    "--output",
                    dest = "output_fp",
                    help = "Filepath to location for output file(s).")

args <- parser$parse_args()

## input file paths
# biom_fp <- './data/misc/euk_filt_mergedDietAim1table_051523-Copy1.qza'
# lacto_asv_fp <- './data/misc/lactoOnlydna-sequences.fasta'

## reading in raw biom table for normalization
biom_pre <- read_qza(args$biom_fp)
biom_pre <- biom_pre$data

## reading in file with the lactococcus asvs
lacto_asv <- read.FASTA(args$lacto_asv_fp)
lacto_asv <- names(lacto_asv)

## filtering lactococcus asvs out of the pre-total sum scaled biom table 
biom_pre %>% 
  as_tibble(rownames = 'asv') %>% 
  filter(!(asv %in% lacto_asv)) -> biom 

## total sum scaling (normalization) steps 
biom %>% 
  gather(-asv, key = sampleid, value = abund) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = abund/sum(abund),
         big_abund = rel_abund*100000,
         big_abund = round(big_abund)) %>% 
  select(asv, sampleid, big_abund) %>% 
  spread(sampleid, big_abund) -> outtab

## writing out total sum scaled table as a .tsv
write_tsv(outtab,
          args$output_fp)
