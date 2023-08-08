## 6-26-23 
## Qiime2 core metrics diversity analysis output plot construction
## unweighted and weighted UniFrac PCoA plots

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
              "cowplot")

chooseCRANmirror(ind = 1)
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)


## input file paths
metadata_FP <- './data/misc/processed_metadata.tsv'
unweighted_FP <- './data/qiime/core_outputs/unweighted_unifrac_pcoa_results.qza'
weighted_FP <- './data/qiime/core_outputs/weighted_unifrac_pcoa_results.qza'

diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')

diet_names_labels <- c('Chow', 
                       'HF/HF', 
                       'HF/LF', 
                       'LF/HF', 
                       'LF/LF')

## functions in order of usage
## 1 
## unweighted/weighted unifrac pcoa result, faith's pd, and shannon entropy file prep 
## going to attempt to return multiple outputs so I can just have one function for file prep
biom_table_prep <- function(unweighted_fp,
                            weighted_fp,
                            faith_fp,
                            shannon_fp,
                            metadata_file){
  ## unweighted pcoa
  unweighted <- read_qza(unweighted_fp)$data
  unweighted_var <- unweighted$ProportionExplained
  unweighted_pcoa <- unweighted$Vectors ##used for pcoa plot
  names(unweighted_pcoa)[names(unweighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## weighted pcoa
  weighted <- read_qza(weighted_fp)$data
  weighted_var <- weighted$ProportionExplained
  weighted_pcoa <- weighted$Vectors
  names(weighted_pcoa)[names(weighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## faith's 
  faith <- read_tsv(faith_fp)
  names(faith)[names(faith) == '#SampleID'] <- 'sampleid'
  ## shannon 
  shannon <- read_tsv(shannon_fp)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  ## unweighted biom 
  unweighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> unweighted_biom
  ## weighted biom
  weighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> weighted_biom
  ## creating a list to return multiple outputs 
  my_list <- list(UnweightedVar = unweighted_var, 
                  WeightedVar = weighted_var,
                  UnweightedBiom = unweighted_biom,
                  WeightedBiom = weighted_biom)
  return(my_list)
}

## 2
## this function will pull out the percent variations from a specified column so you can add it to your pcoa plots 
pcoa_ax_lab <- function(unifrac_var, col_name){
  uni_lab <- as.character(round(unifrac_var[col_name] * 100, 2))
  uni_lab <- paste0(col_name, ' - ', uni_lab, '%')
  return(uni_lab)
}

## 3
## pcoa plot function
## xlab and ylab are outputs from pcoa_ax_lab function
pcoa_plot <- function(biom_file,
                      labels,
                      names_labels,
                      xlab,
                      ylab,
                      title){
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  biom_file %>% 
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = faith_pd), pch = 21, alpha = 0.7) +
    theme_bw(base_size = 14) +
    scale_fill_distiller(palette = 'Spectral', name = "Faith's PD") +
    facet_grid(day_post_inf~diet, 
               labeller = labeller(diet = labs)) +
    theme(legend.text = element_text(size = 8.5),
          strip.text.y = element_text(angle = 0)) +
    ggtitle(title) +
    labs(x = xlab, y = ylab) -> pcoa
  return(pcoa)
}

## core metrics file prep
## metadata file prep
metadata <- read_tsv(metadata_FP)

## preparing core beta diversity files for ggplot
core_files <- biom_table_prep(unweighted_FP,
                              weighted_FP,
                              faith_pd_FP,
                              shannon_FP,
                              metadata)
## extracting core beta diversity files from named list 
uw_var <- core_files$UnweightedVar
w_var <- core_files$WeightedVar
unweighted_biom <- core_files$UnweightedBiom
weighted_biom <- core_files$WeightedBiom

## unweighted pcoa
uw_uni_xlab <- pcoa_ax_lab(uw_var, 'PC1')
uw_uni_ylab <- pcoa_ax_lab(uw_var, 'PC2')

uw_title <- 'Total Sum Scaled Unweighted UniFrac'

unweighted_pcoa <- pcoa_plot(unweighted_biom,
                             diet_labs,
                             diet_names_labels,
                             uw_uni_xlab,
                             uw_uni_ylab,
                             uw_title)
## weighted pcoa 
w_uni_xlab <- pcoa_ax_lab(w_var, 'PC1')
w_uni_ylab <- pcoa_ax_lab(w_var, 'PC2')

w_title <- 'Total Sum Scaled Weighted UniFrac'

weighted_pcoa <- pcoa_plot(weighted_biom,
                           diet_labs,
                           diet_names_labels,
                           w_uni_xlab,
                           w_uni_ylab,
                           w_title)

## saving my plot outputs to the plots folder
ggsave("unweighted_unifrac_pcoa.pdf",
       plot = unweighted_pcoa, 
       width = 12, 
       height = 5, 
       path = './plots')

ggsave("weighted_unifrac_pcoa.pdf",
       plot = weighted_pcoa, 
       width = 12, 
       height = 5, 
       path = './plots')