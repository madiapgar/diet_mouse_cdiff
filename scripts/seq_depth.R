## 6-29-23
## calculating the sequencing depth for each sample for later statistical use

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

## reading in pre-total sum scaled file 
biom <- read_qza(args$biom_fp)
biom <- biom$data 

## filtering lactococcus asvs out of my pre-total sum biom table 
## reading in lactococcus asvs and extracting them 
lacto_asv <- read.FASTA(args$lacto_asv_fp)
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
          args$output_fp)