## 6-26-23
## Qiime2 core diversity analysis statistical analysis for beta diversity metrics

## needed libraries
packages <- c("ape", 
              "ggpubr", 
              "magrittr", 
              "qiime2R", 
              "tidyverse", 
              "broom", 
              "rstatix",
              "ggh4x",
              "vegan",
              "viridis",
              "cowplot",
              "dunn.test")

chooseCRANmirror(ind = 1)
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)


## input file paths
metadata_FP <- './data/misc/processed_metadata.tsv'
uw_dist_fp <- './data/qiime/core_outputs/uw_dist_matrix.tsv'
w_dist_fp <- './data/qiime/core_outputs/w_dist_matrix.tsv'
unwanted_samples <- c('Mock20220615A', 'Mock_1A', 'Mock_2A',
                      'Mock_3A', 'Mock_4A', 'Mock_5A', 'Mock_6A',
                      'Mock_7A', 'PCR Blank0',
                      'PCR Blank1', 'Mock_7', 'Mock_6',
                      'Mock_5', 'Mock_4', 'PCR blank')

## functions in order of usage 
##1 
## for distance matrix processing
## for beta diversity statistical analysis 
dist_matrix_prep <- function(metadata_file,
                             dist_matrix_fp,
                             sample_filter){ 
  ## metadata filtering
  metadata_file %>% 
    filter(!(sampleid %in% sample_filter)) -> metadata
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


## 4 
## beta diversity adonis2 testing function
adonis_test <- function(dist_matrix,
                        metadata_file){
  adonis_results <- adonis2(as.dist(dist_matrix) ~ purified_diet * seq_depth + high_fat + high_fiber +
                              day_post_inf + study,
                            data = metadata_file,
                            permutations = 999, 
                            parallel = 4)
  adonis_results <- tidy(adonis_results)
  return(adonis_results)
}

## actually using the functions
## reading in metadata file 
meta <- read_tsv(metadata_FP)

## weighted unifrac 
w_dist_files <- dist_matrix_prep(meta,
                                 w_dist_fp,
                                 unwanted_samples)

w_dist <- w_dist_files$DistanceMatrix
stat_meta <- w_dist_files$Metadata

w_adonis <- adonis_test(w_dist,
                        stat_meta)


## unweighted unifrac
uw_dist_files <- dist_matrix_prep(meta,
                                  uw_dist_fp,
                                  unwanted_samples)

uw_dist <- uw_dist_files$DistanceMatrix
stat_meta <- uw_dist_files$Metadata

uw_adonis <- adonis_test(uw_dist,
                         stat_meta)

## writing results out as a .tsv file 
write_tsv(w_adonis, './stats/w_adonis_results.tsv')
write_tsv(uw_adonis, './stats/uw_adonis_results.tsv')