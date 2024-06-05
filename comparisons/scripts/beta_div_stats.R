## 6-26-23
## Qiime2 core diversity analysis statistical analysis for beta diversity metrics

## needed libraries
library(ggpubr)
library(ggplot2)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(broom)
library(cowplot)
library(rstatix)
library(ape)
library(vegan)
library(ggh4x)
library(viridis)
library(argparse)

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-uu",
                    "--uu_dist",
                    dest = "uu_dist_fp",
                    help = "Filepath to Unweighted UniFrac Distance Matrix in .tsv format.")
parser$add_argument("-wu",
                    "--wu_dist",
                    dest = "wu_dist_fp",
                    help = "Filepath to Weighted UniFrac Distance Matrix in .tsv format")
parser$add_argument("-ua",
                    "--uu_adonis",
                    dest = "uu_adonis_fp",
                    help = "Filepath to Unweighted UniFrac Adonis test results in .tsv format.")
parser$add_argument("-wa",
                    "--wu_adonis",
                    dest = "wu_adonis_fp",
                    help = "Filepath to Weighted UniFrac Adonis test results in .tsv format.")

args <- parser$parse_args()



## functions in order of usage 
##1 
## for distance matrix processing
## for beta diversity statistical analysis 
dist_matrix_prep <- function(metadata_file,
                             dist_matrix_fp){
  ## metadata
  metadata <- metadata_file
  ## distance matrix
  dist <- read_tsv(dist_matrix_fp)
  names(dist)[names(dist) == '...1'] <- 'sampleid'
  dist %>% 
    gather(-sampleid, key = sample_col, value = dist) %>% 
    filter(sampleid %in% metadata$sampleid) %>% 
    filter(sample_col %in% metadata$sampleid) %>% 
    spread(sample_col, dist) -> dist_long
  dist_long %>% 
    select(-sampleid) -> dist_proc
  metadata %>% 
    arrange(sampleid) -> metadata
  metadata %>% 
    filter(sampleid %in% dist_long$sampleid) -> filt_meta
  dist_proc <- as.matrix(dist_proc)
  row.names(dist_proc) <- colnames(dist_proc)
  filt_meta <- filt_meta[order(filt_meta$sampleid),]
  ## list of outputs
  my_list <- list(Metadata = filt_meta,
                  DistanceMatrix = dist_proc)
  return(my_list)
}


## 2
## beta diversity adonis2 testing function
adonis_test <- function(dist_matrix,
                        metadata_file){
  adonis_results <- adonis2(as.dist(dist_matrix) ~ experiment_set + vendor + mouse_sex,
                            data = metadata_file,
                            permutations = 999, 
                            parallel = 4)
  adonis_results <- tidy(adonis_results)
  adonis_results['signif'] <- symnum(adonis_results$p.value,
                                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                     symbols = c("****", "***", "**", "*", "ns"),
                                     abbr.colnames = FALSE,
                                     na = "")
  return(adonis_results)
}


## actually using the functions
## reading in metadata file 
meta <- read_tsv(args$metadata_FP)
names(meta)[names(meta) == '#SampleID'] <- 'sampleid'

## weighted unifrac 
w_dist_files <- dist_matrix_prep(meta,
                                 args$wu_dist_fp)

w_dist <- w_dist_files$DistanceMatrix
stat_meta <- w_dist_files$Metadata

w_adonis <- adonis_test(w_dist,
                        stat_meta)


## unweighted unifrac
uw_dist_files <- dist_matrix_prep(meta,
                                  args$uu_dist_fp)

uw_dist <- uw_dist_files$DistanceMatrix
stat_meta <- uw_dist_files$Metadata

uw_adonis <- adonis_test(uw_dist,
                         stat_meta)


## writing results out as a .tsv file 
write_tsv(w_adonis, args$wu_adonis_fp)
write_tsv(uw_adonis, args$uu_adonis_fp)