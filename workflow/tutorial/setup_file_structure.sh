#!/bin/bash

## creating a function to read in lists of file paths to make my life easier
read_in_files () {
    for link in ${*} ;
        do
            wget ${link}
        done
}

 workflow_list=("https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/snakefile" 
                "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/run_snakemake.sh"
                )

rules_list=("https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/01_demux_dada2.smk" 
            "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/02_phylogeny.smk" 
            "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/03_tss.smk" 
            "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/04_core_metrics.smk" 
            "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/05_microbiome_r.smk" 
            "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/rules/06_madis_cecalAnalysis.smk"
            )

envs_list=("https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/envs/r_env.yml" 
           "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/envs/install_envs_linux.sh" 
           "https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/envs/install_envs_macos.sh"
           )

config_list=("https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/templates/config_template.yml")


## creating a directory for the analysis
mkdir madis_16s_workflow

## navigate into said directory
cd ./madis_16s_workflow

## creating my_data and workflow subdirectories
mkdir workflow

## navigating into workflow subdirectory
cd ./workflow
read_in_files ${workflow_list[*]}

## creating subdirectories in workflow 
mkdir config_files
mkdir envs
mkdir rules

## download files into rules folder
cd ./rules
read_in_files ${rules_list[*]}

## download files into envs folder
cd ../envs
read_in_files ${envs_list[*]}

## download files into config_files foler
cd ../config_files
read_in_files ${config_list[*]}




