## 7-19-23
## this script contains a function that filters the ko meta contrib file from picrust2
## to only contain the desired kos (since that file is massive)

## needed libraries
library(ggpubr)
library(magrittr)
library(tidyverse)
library(broom)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-i",
                    "--ko_in",
                    dest = "ko_in_fp",
                    help = "Filepath for metagenome contribution file in .tsv.gz format.")
parser$add_argument("-o",
                    "--ko_out",
                    dest = "ko_out_fp",
                    help = "Filepath for trimmed metagenome contribution output file in .tsv format.")

args <- parser$parse_args()

## input file paths and others
# ko_in <- './data/picrust/out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz'
# ko_out <- './data/picrust/meta_contrib.tsv'
wanted_kos <- c('K00929', 'K01034','K15873', 'K15874')

## function
contrib_red <- function(in_fp, 
                        out_fp,
                        ko_list){
  kometa_contrib <- read_tsv(file = in_fp)
  names(kometa_contrib)[names(kometa_contrib) == 'function'] <- 'ko'
  kometa_contrib %>% 
    filter(ko %in% ko_list) -> min_kometa_contrib
  write_tsv(min_kometa_contrib, out_fp)
}

## using the function
contrib_red(args$ko_in_fp,
            args$ko_out_fp,
            wanted_kos)
