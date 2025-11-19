## 3-6-24
## this script processes the raw toxin, metabolomics, histopathology, and bile acid files
## can use this outputs to create plots and do comparisons between them

## needed libraries 
library(ggpubr)
library(magrittr)
library(tidyverse)
library(broom)
library(argparse)

#### NOTE!! ####
## study 1 has been filtered out of ALL RESULTS due to lack of CDI!! ##
## madi redid this analysis on 11-18-2025 due to low sample numbers for high fiber diets ##
################

## using argparse for my file paths
## so I can easily edit file paths from my workflow and not have to edit the actual R scripts
parser <- ArgumentParser()
parser$add_argument("-m",
                    "--metadata",
                    dest = "metadata_FP",
                    help = "Filepath to metadata file in .tsv format.")
parser$add_argument("-t",
                    "--toxin",
                    dest = "toxin_FP",
                    help = "Filepath to toxin file in .tsv format.")
parser$add_argument("-h",
                    "--histopathology",
                    dest = "histo_FP",
                    help = "Filepath to histopathology file in .csv format.")
parser$add_argument("-mb",
                    "--metabolomics",
                    dest = "metab_FP",
                    help = "Filepath to metabolomics file in .csv format.")
parser$add_argument("-b",
                    "--bile_acid",
                    dest = "bile_acid_FP",
                    help = "Filepath to corrected bile acid file in .tsv format.")
parser$add_argument("-nto",
                    "--neat_toxin_out",
                    dest = "neatToxin_output_FP",
                    help = "Filepath to location for neat toxin output file in .tsv format.")
parser$add_argument("-dto",
                    "--dil_toxin_out",
                    dest = "dilToxin_output_FP",
                    help = "Filepath to location for diluted toxin output file in .tsv format.")
parser$add_argument("-mo",
                    "--metab_out",
                    dest = "metab_output_FP",
                    help = "Filepath to metabolomics output file in .tsv format.")
parser$add_argument("-ho",
                    "--histo_out",
                    dest = "histo_output_FP",
                    help = "Filepath to histopathology output file in .tsv format.")
parser$add_argument("-bo",
                    "--bile_acid_out",
                    dest = "bileAcid_output_FP",
                    help = "Filepath to bile acid output file in .tsv format.")
parser$add_argument("-br",
                    "--bile_ratio_out",
                    dest = "bileAcid_ratio_FP",
                    help = "Filepath to bile acid ratio output file in .tsv format.")


args <- parser$parse_args()

## input file paths
metadata_FP <- '../data/misc/filt_cecal_processed_metadata.tsv'
mouseid_facil_FP <- '../data/misc/mouseID_facil.tsv'
mouseid_diet_key_FP <- '../data/misc/cecal_key.txt'
toxin_FP <- '../data/misc/toxin_final_data.tsv'
histo_FP <- '../data/misc/histo_data.csv'
metab_FP <- '../data/misc/metabolomics.csv'
bile_acid_FP <- '../data/misc/corrected_bile_acid.tsv'

## lists of wanted and unwanted columns for metabolomics processing
wanted_metabs <- c('Acetic Acid (ug/g)',
                   'Propanoic Acid (ug/g)',
                   'n-Butanoic Acid (ug/g)')
unwanted_columns <- c('2-methyl-propanoic acid (ug/g)',
                      'Isopentanoic Acid (ug/g)',
                      '2-methyl-Butanoic Acid (ug/g)',
                      'Pentanoic Acid (ug/g)',
                      'Notes',
                      'Sample Group',
                      'SCFA Data File',
                      'Acq. Date-Time',
                      'Tube_Label',
                      'Collection Date',
                      'Dil.')
## inhibitors and promoters list for bile acid file processing
inhibitor_list <- c('_a-MCA',
                    '_b-MCA',
                    '_LCA',
                    '_DCA')
promoter_list <- c('_T-CA',
                   '_CA')

## needed functions
## 1 - toxin file prep
toxin_prep <- function(mouseid_facil_fp,
                       toxin_fp){
  ## toxin data 
  toxin <- read_tsv(toxin_FP)
  wanted_ids <- toxin$mouse_id
  ## metadata 
  mouseid_facil <- read_tsv(mouseid_facil_fp)
  
  ## toxin neat concentrations (non_diluted)
  toxin %>% 
    select(!c('Total TcA 1:10', 'Total TcB 1:10')) %>% 
    gather('Total TcA Neat', 'Total TcB Neat', 
           key = neat_toxin, value = neat_conc) -> pre_neat_toxin
  
  pre_neat_toxin$neat_conc[pre_neat_toxin$neat_conc == 'BDL'] <- '0'
  
  pre_neat_toxin %>% 
    mutate(neat_conc = as.numeric(neat_conc),
           diet = case_when(
             Sample_Type == 'Chow' ~ 'Chow',
             Sample_Type == 'Ctrl' ~ 'LF/LF',
             Sample_Type == 'WD' ~ 'HF/LF',
             Sample_Type == 'WD+F' ~ 'HF/HF',
             Sample_Type == 'Ctrl+F' ~ 'LF/HF'
           )) %>% 
    select(!c('Tube_Label',
              'Collection Date',
              'Sample_Type',
              'Extra_Sample')) %>% 
    left_join(mouseid_facil, by = "mouse_id") %>% 
    filter(study != 1) -> neat_toxin
  
  ## toxin diluted concentrations 
  ## chow is not included in this 
  toxin %>% 
    select(!c('Total TcA Neat', 'Total TcB Neat')) %>% 
    gather('Total TcA 1:10', 'Total TcB 1:10',
           key = dil_toxin, value = dil_conc) -> pre_dil_toxin
  
  pre_dil_toxin$dil_conc[pre_dil_toxin$dil_conc == 'BDL'] <- '0'
  pre_dil_toxin$dil_conc[pre_dil_toxin$dil_conc == 'Chow'] <- '0'
  
  pre_dil_toxin %>% 
    mutate(dil_conc = as.numeric(dil_conc),
           diet = case_when(
             Sample_Type == 'Chow' ~ 'Chow',
             Sample_Type == 'Ctrl' ~ 'LF/LF',
             Sample_Type == 'WD' ~ 'HF/LF',
             Sample_Type == 'WD+F' ~ 'HF/HF',
             Sample_Type == 'Ctrl+F' ~ 'LF/HF'
           )) %>% 
    filter(diet != 'Chow') %>%
    select(!c('Tube_Label',
              'Collection Date',
              'Sample_Type',
              'Extra_Sample')) %>% 
    left_join(mouseid_facil, by = "mouse_id") %>% 
    filter(study != 1) -> dil_toxin
  ## creating a list of my outputs
  my_list <- list(NeatToxin = neat_toxin,
                  DilToxin = dil_toxin)
  return(my_list)
}

## 2 - metabolomics file prep
metab_prep <- function(mouseid_facil_fp,
                       metab_fp,
                       metab_col_filter,
                       metab_filter){
  ## metabolomics file prep 
  metab <- read_csv(metab_FP)
  ## metadta file prep 
  mouseid_facil <- read_tsv(mouseid_facil_fp) 
  
  metab %>% 
    select(-(all_of(metab_col_filter))) %>% 
    gather(all_of(metab_filter), key = metabolite, value = concentration) -> pre_metab

  ## changes all 'ND' values in the concentration column to 0 
  pre_metab$concentration[pre_metab$concentration == 'ND'] <- 0
  
  pre_metab %>% 
    mutate(concentration = ifelse(is.na(concentration), "not_enough_sample", concentration),
           diet = case_when(
             Sample_Type == 'Chow' ~ 'Chow',
             Sample_Type == 'Ctrl' ~ 'LF/LF',
             Sample_Type == 'WD' ~ 'HF/LF',
             Sample_Type == 'WD+F' ~ 'HF/HF',
             Sample_Type == 'Ctrl+F' ~ 'LF/HF'
           )) %>% 
    select(-Sample_Type) %>% 
    left_join(mouseid_facil, by = "mouse_id") %>% 
    filter(study != 1) -> big_metab
  
  return(big_metab)
}

## 3 - histopathology file prep
histo_prep <- function(mouseid_facil_fp,
                       mouseid_diet_fp,
                       histo_fp){
  ## reading in metadata file(s)
  mouseid_facil <- read_tsv(mouseid_facil_fp)
  mouseid_diet_key <- read_tsv(mouseid_diet_fp)
  ## reading in histo file
  histo <- read_csv(histo_fp) %>% 
    filter(!is.na(mouse_id)) %>% 
    gather(cecum, colon, key = tissue, value = score) %>% 
    left_join(mouseid_facil, by = "mouse_id") %>% 
    left_join(mouseid_diet_key, by = "mouse_id") %>% 
    select(!c('date', 'tube_numb')) %>% 
    filter(study != 1) -> big_histo 

  return(big_histo)
}

## 4 - bile acid file prep
bile_acid_prep <- function(bile_acid_fp,
                           inhibitors,
                           promoters){
  ## reading in bile acid file
  bile_acid <- read_tsv(bile_acid_fp)
  ## putting together inhibitor table
  bile_acid %>% 
    select(diet, tube_numb, tube_label, date, 
           mouse_id, facility, study, contains(unlist(inhibitors))) %>% 
    gather(contains('acid'), key = bile_acid, value = concentration) %>% 
    mutate(c_diff_effect = paste('inhibitor')) %>% 
    filter(concentration != '#VALUE!') -> inhibit_bile
  ## putting together promoter table
  bile_acid %>% 
    select(diet, tube_numb, tube_label, date, 
           mouse_id, facility, study, contains(unlist(promoters))) %>% 
    gather(contains('acid'), key = bile_acid, value = concentration) %>% 
    mutate(c_diff_effect = paste('promoter')) -> promote_bile
  ## putting tables back together
  big_bile_acid <- rbind(inhibit_bile,
                         promote_bile)
  ## changing all undetectable values to 0 
  big_bile_acid$concentration[big_bile_acid$concentration == '<LOD' | big_bile_acid$concentration == '<LOQ'] <- 0
  ## adding 2 to all values so that all points will be on the plot when transformed to log10
  big_bile_acid %>%
    mutate(concentration = ifelse(is.na(concentration), "not_enough_sample", concentration),
           conc_normalized = ifelse(concentration != "not_enough_sample", as.numeric(concentration) + 2, concentration)) %>% 
    filter(study != 1) -> big_bile_acid
  ## creating list of outputs
  my_list <- list(NonProcBile = bile_acid,
                  ProcBile = big_bile_acid)
  return(my_list)
}

## 5 - bile acid ratio file prep
bile_ratio_prep <- function(proc_bile_table){
  ## summing the normalized concentration of promoters to inhibitors by mouse id for each diet 
  proc_bile_table %>% 
    select(-bile_acid) %>% 
    group_by(diet, c_diff_effect, mouse_id) %>% 
    filter(conc_normalized != "not_enough_sample") %>% 
    summarise(sum = sum(as.numeric(conc_normalized))) %>% 
    ungroup() -> bile_sum
  ## them dividing the normalized concentration of promoters/inhibitors for each mouse id by diet
  bile_sum %>% 
    spread(c_diff_effect, sum) %>% 
    group_by(diet, mouse_id) %>% 
    mutate(ratio = (promoter/inhibitor)) %>% 
    ungroup() %>% 
    mutate(ratio_label = paste('C. difficile Promoter:Inhibitor')) -> bile_ratio
  return(bile_ratio)
}


## processing my files so I can use them for plot construction and statistical analysis
## toxin
toxin_files <- toxin_prep(mouseid_facil_FP,
                          toxin_FP)
neat_toxin <- toxin_files$NeatToxin
dil_toxin <- toxin_files$DilToxin

## metabolomics 
metab <- metab_prep(mouseid_facil_FP,
                    metab_FP,
                    unwanted_columns,
                    wanted_metabs)
## histopathology
histo <- histo_prep(mouseid_facil_FP,
                    mouseid_diet_key_FP,
                    histo_FP)
## bile acid 
bile_files <- bile_acid_prep(bile_acid_FP,
                             inhibitor_list,
                             promoter_list)
proc_bile_acid <- bile_files$ProcBile
noProc_bile_acid <- bile_files$NonProcBile

## bile acid ratio (dependent on proc_bile_acid file generated above)
bile_ratio <- bile_ratio_prep(proc_bile_acid)

## saving my processed files as a .tsv so I can continue to use them 
neatToxin_output_FP <- '../data/misc/v2_proc_neat_toxin.tsv'
dilToxin_output_FP <- '../data/misc/v2_proc_dil_toxin.tsv'
metab_output_FP <- '../data/misc/v2_proc_metab.tsv'
histo_output_FP <- '../data/misc/v2_proc_histo.tsv'
bileAcid_output_FP <- '../data/misc/v2_proc_bile_acid.tsv'
bileAcid_ratio_FP <- '../data/misc/v2_proc_ratio_bileAcid.tsv'

write_tsv(neat_toxin,
          neatToxin_output_FP)
write_tsv(dil_toxin,
          dilToxin_output_FP)
write_tsv(metab,
          metab_output_FP)
write_tsv(histo,
          histo_output_FP)
write_tsv(proc_bile_acid,
          bileAcid_output_FP)
write_tsv(bile_ratio,
          bileAcid_ratio_FP)
