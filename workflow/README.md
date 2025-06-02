# 16S rRNA Data Analysis Pipeline for this Study 

This directory contains the snakemake workflow written for the upper 16S rRNA sequencing data analysis. The workflow heavily uses QIIME2 with some small deviations and while it does have the capacity to perform some basic downstream analysis in R, it's currently not set up to produce the plots/stats featured in the paper. 

> [!IMPORTANT]
> If you have pulled the sequencing data off of QIITA and want to run it through this pipeline, I would recommend using the barcodes files located [here](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/arizona_experiments/stool/data/misc) for SEQ040, SEQ044, and SEQ045 due to sample id discrepancies! All other sequencing runs should be fine. 

## Directory contents are as follows:

- [**config_files:**](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/config_files) Config files containing inputs needed to run the pipeline for the associated data analysis
- [**envs:**](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/envs) Conda environment .yaml files and bash scripts for easy install of the version of QIIME2 used in this analysis
- [**rules:**](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/rules) Sub-workflows written for pipeline modularization that go in order of the number in the file name (you shouldn't have to touch these but you can look at them)
- [**templates:**](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/templates) Includes templates for basic downstream R analysis (see tutorial) and a blank config file template with additional information
- [**tutorial:**](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/tutorial) The written tutorial for pipeline install/usage and associated scripts

## Other scripts in this directory:

- [**run_snakemake.sh:**](https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/run_snakemake.sh) A bash script that has the command to run snakemake with the necessary parameters to run the pipeline 
- [**snakefile:**](https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/snakefile) The master script of the pipeline that pulls together inputs from the config file and sub-workflows in the `rules` directory

