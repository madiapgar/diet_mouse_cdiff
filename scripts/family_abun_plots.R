## 7-14-23
## putting together my microbe family abundance plots based on the 16S rDNA data
## a different way of visualizing the taxa barplot (might be easier)
## you have two options

## needed libraries
packages <- c("ggpubr", 
              "magrittr", 
              "qiime2R", 
              "tidyverse", 
              "broom",
              "cowplot",
              "vegan",
              "viridis")

install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)

## input file paths and others
otu_table_FP <- './data/qiime/taxonomy_filtered.qza'
tax_FP <- './data/qiime/taxonomy.qza'
metadata_FP <- './data/misc/processed_metadata.tsv'
diet_labs <- 
  c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')
names(diet_labs) <- c('Chow', 
                      'HF/HF', 
                      'HF/LF', 
                      'LF/HF', 
                      'LF/LF')
wanted_level <- 'Family'
wanted_family <- c('Enterobacteriaceae', 'Lactobacillaceae', 'Lachnospiraceae', 'Enterococcaceae',
                   'Staphylococcaceae', 'Tannerellaceae', 'Muribaculaceae', 'Bacteroidaceae', 
                   'Marinifilaceae')

## functions in order of usage 
## 1
family_abun_file_prep <- function(metadata_fp,
                                  tax_fp,
                                  otu_table_fp,
                                  tax_level,
                                  wanted_tax){
  ## metadata
  metadata <- read_tsv(metadata_FP)
  ## taxonomy
  taxonomy <- read_qza(tax_FP)$data %>% 
    parse_taxonomy() %>% 
    rownames_to_column('asv')
  ## otu table 
  otu_table <- read_qza(otu_table_FP)$data
  otu_table %>% 
    as_tibble(rownames = 'asv') %>% 
    gather(-asv, key = sampleid, value = abun) %>% 
    group_by(sampleid) %>% 
    mutate(rel_abun = abun/sum(abun)) %>% 
    mutate(p_abun = rel_abun + 0.000001) -> otu_table
  ## joining all tables together 
  otu_table %>% 
    left_join(metadata, by = 'sampleid') %>% 
    left_join(taxonomy, by = 'asv') -> abun_table
  abun_table %>% 
    group_by(sampleid, day_post_inf, diet, mouse_id, 
             purified_diet, high_fat, high_fiber, 
             seq_depth, .data[[tax_level]]) %>% 
    summarise(rel_abund = sum(rel_abun)) %>% 
    filter(.data[[tax_level]] %in% wanted_tax) %>% 
    mutate(mouse_fact = as.factor(mouse_id),
           day_fact = as.factor(day_post_inf)) -> abun_filt
  ## creating a list for my outputs
  my_list <- list(Metadata = metadata,
                  Taxonomy = taxonomy,
                  OTUTable = otu_table,
                  AbundanceTable = abun_filt)
  return(my_list)
}

## 2 
abun_plots <- function(abundance_table){
  ## first plot
  abundance_table %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = rel_abund)) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
    geom_line(aes(group = mouse_id), alpha = 0.3) +
    geom_smooth(se = FALSE) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    theme_bw(base_size = 16) +
    facet_grid(Family~diet, labeller = labeller(diet = diet_labs)) +
    theme(strip.text.y = element_text(angle = 0)) +
    ggtitle("Microbe Family Relative Abundance") +
    ylab("Relative Abundance") +
    xlab("Days Relative to Infection") -> family_abun1
  ## second plot
  abun_filt %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = mouse_fact, y = day_fact)) +
    geom_tile(aes(fill = rel_abund), alpha = 0.5) +
    scale_fill_viridis(option = "H", name = 'Relative\nAbundance') +
    theme_bw(base_size = 16) +
    facet_grid(Family~diet, scales = 'free',
               labeller = labeller(diet = diet_labs)) +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.x = element_blank()) + 
    xlab("Mouse ID") +
    ylab("Days Relative to Infection") +
    scale_y_discrete(limits = rev) -> family_abun2
  ## creating a list of my two plots
  my_list <- list(FamilyAbundance1 = family_abun1,
                  FamilyAbundance2 = family_abun2)
  return(my_list)
}

## family abundance table prep 
abun_files <- family_abun_file_prep(metadata_FP,
                                    tax_FP,
                                    otu_table_FP,
                                    wanted_level,
                                    wanted_family)

## pulling the abundance table out, you can also take metadata, otu table, and taxonomic info
## out of the list output 
abun_filt <- abun_files$AbundanceTable

## generating the plots
family_abun_plots <- abun_plots(abun_filt)

## scatter plot/ line graph by mouse id and wanted microbe families
abun1 <- family_abun_plots$FamilyAbundance1
## heat map by days relative to infection and mouse id for wanted microbe families 
abun2 <- family_abun_plots$FamilyAbundance2

## saving my plot outputs
## option #1
ggsave("family_abun1.pdf",
       plot = abun1, 
       width = 15, 
       height = 10, 
       path = './plots')
## option #2
ggsave("family_abun2.pdf",
       plot = abun2, 
       width = 25, 
       height = 10, 
       path = './plots')