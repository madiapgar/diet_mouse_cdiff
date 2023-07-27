## 7-19-23
## this script contains a function that filters the ko meta contrib file from picrust2
## to only contain the desired kos (since that file is massive)

## needed libraries
packages <- c("ggpubr", 
              "magrittr",
              "tidyverse", 
              "broom")

for(package in packages){
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

## input file paths and others
ko_in <- './data/picrust/out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz'
ko_out <- './data/picrust/tss3_meta_contrib.tsv'
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
contrib_red(ko_in,
            ko_out,
            wanted_kos)
