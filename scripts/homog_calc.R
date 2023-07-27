## 6-26-23
## homogeneity calculations for unweighted and weighted UniFrac distance matrices
## outputs a plot and the statistical results

## needed libraries
packages <- c("ggpubr", 
              "magrittr", 
              "qiime2R", 
              "tidyverse", 
              "broom",
              "cowplot")

for(package in packages){
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

## input file paths
metadata_FP <- './data/misc/processed_metadata.tsv'
uu_dist_fp <- './data/qiime/core_outputs/uw_dist_matrix.tsv'
wu_dist_fp <- './data/qiime/core_outputs/w_dist_matrix.tsv'

## lists to redo the diet names on the facet labels of the ggplot created below 
diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')

names(diet_labs) <- c('Chow', 'HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')

## functions in order of usage
## 1 
## this function combines your metadata file and your distance matrix 
## it also pulls out the wanted diet variable from the metadata file 
## requires the input of a distance matrix (dist_mat), metadata file (meta),
## and the list of unique diets (mydiet)
## used in subsequent resil/homog_diet_assembly functions
dist_filter <- function(dist_mat, meta, mydiet){
  meta %>% 
    filter(diet == mydiet) -> min_meta
  
  dist_mat <- data.frame(dist_mat)
  rownames(dist_mat) <- dist_mat[,1]
  dist_mat <- dist_mat[,-1]
  
  sample_idx <- rownames(dist_mat) %in% min_meta$sampleid 
  
  min_dist_mat <- dist_mat[sample_idx, sample_idx]
  min_dist_mat <- as_tibble(min_dist_mat, rownames = 'sampleid')
  colnames(min_dist_mat)[2:ncol(min_dist_mat)] <- min_dist_mat$sampleid
  return(min_dist_mat)
}

## 2
## this function takes the output from dist_filter and the metadata file 
## and pulls all unique distance matrix comparisons out for that particular diet 
## and places them in a new column 
## it also esta the baseline that all sample pairs will be compared to 
## requires the distance matrix output from the dist_filter (min_dist_mat),
## and the metadata file (meta)
homog_dist_comp <- function(min_dist_mat, meta){
  meta %>% 
    select(sampleid, diet, day_post_inf, high_fat, high_fiber, purified_diet,
           seq_depth) -> min_meta
  
  min_dist_mat %>% 
    rename(sampleid_x = sampleid) %>% 
    gather(-sampleid_x, key = sampleid_y, value = dist) %>% 
    merge(min_meta, by.x = 'sampleid_x', by.y = 'sampleid') %>% 
    merge(min_meta, by.x = 'sampleid_y', by.y = 'sampleid') -> long_dist_mat
  
  long_dist_mat %>% 
    mutate(key = paste0(pmin(sampleid_x, sampleid_y), 
                        pmax(sampleid_x, sampleid_y), 
                        sep = "_")) %>% 
    distinct(key, .keep_all = TRUE) %>% 
    filter(sampleid_x != sampleid_y,
           ## for homogenicity testing, do day_post_inf.x == day_post_inf.y 
           ## and don't filter out day -15 
           day_post_inf.x == day_post_inf.y) %>% 
    select(-diet.y) %>% 
    rename(diet = diet.x) -> long_dist_mat
  
  return(long_dist_mat)
}

## 3
## this is a for loop that does the above analysis for each diet variable and places them 
## in a list of dataframes that can be bound together to create one plot for comparison 
## requires the input of a distance matrix (dist_mat), metadata file (meta),
## and the list of unique diets (mydiet)
homog_diet_assembly <- function(dist_mat, 
                                meta,
                                mydiet){
  homog <- list()
  for (i in 1:length(mydiet)){
    tmp_dist_mat <- dist_filter(dist_mat, meta, mydiet[i])
    tmp_homog <- homog_dist_comp(tmp_dist_mat, meta)
    homog[[i]] <- tmp_homog
  }
  homog <- bind_rows(homog)
  return(homog)
}

## file prep for homog plot and stats
meta <- read_tsv(metadata_FP)
uu_dist <- read_tsv(uu_dist_fp)
wu_dist <- read_tsv(wu_dist_fp)

## pulling out unique diets 
diets <- unique(meta$diet)

## homogeneity calculations
## unweighted unifrac
uu_homog <- homog_diet_assembly(uu_dist,
                                meta,
                                diets)

## weighted unifrac
wu_homog <- homog_diet_assembly(wu_dist,
                                meta,
                                diets)

## unweighted UniFrac plot and stats
## plot
uu_homog %>% 
  ggplot(aes(x = day_post_inf.y, y = dist)) +
  geom_boxplot((aes(group = day_post_inf.y)), outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_smooth(se = FALSE, method = 'loess', color = 'blue') + 
  scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) + 
  scale_y_reverse() +
  facet_wrap(~diet, 
             labeller = labeller(diet = diet_labs),
             nrow = 1) +
  theme_bw(base_size = 14) +
  xlab('Days Relative to Infection') +
  ylab('Unweighted UniFrac Distance\n(Pairwise Distances within Day)') +
  ggtitle("Microbiome Homogeneity Over Time") -> uu_homog_plot

## stats
uu_homog %>% 
  group_by(day_post_inf.y) %>% 
  do(tidy(lm(dist ~ (purified_diet.y * seq_depth.y) + high_fat.y + high_fiber.y,
             data = .))) %>% 
  filter(p.value <= 0.05) -> uu_homog_results

## weighted UniFrac plot and stats
## plot
wu_homog %>% 
  ggplot(aes(x = day_post_inf.y, y = dist)) +
  geom_boxplot((aes(group = day_post_inf.y)), outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_smooth(se = FALSE, method = 'loess', color = 'blue') + 
  scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) + 
  scale_y_reverse() +
  facet_wrap(~diet, 
             labeller = labeller(diet = diet_labs),
             nrow = 1) +
  theme_bw(base_size = 14) +
  xlab('Days Relative to Infection') +
  ylab('Weighted UniFrac Distance\n(Pairwise Distances within Day)') +
  ggtitle("Microbiome Homogeneity Over Time") -> wu_homog_plot

## stats
wu_homog %>% 
  group_by(day_post_inf.y) %>% 
  do(tidy(lm(dist ~ (purified_diet.y * seq_depth.y) + high_fat.y + high_fiber.y,
             data = .))) %>% 
  filter(p.value <= 0.05) -> wu_homog_results

## saving my output plots and stats
## plots
ggsave("wu_homogeneity.pdf", 
       plot = wu_homog_plot,
       width = 11, 
       height = 4,
       path = './plots')
ggsave("uu_homogeneity.pdf", 
       plot = uu_homog_plot,
       width = 11, 
       height = 4,
       path = './plots')

## stats 
write_tsv(wu_homog_results,
          './stats/wu_homogeneity.tsv')
write_tsv(uu_homog_results,
          './stats/uu_homogeneity.tsv')

