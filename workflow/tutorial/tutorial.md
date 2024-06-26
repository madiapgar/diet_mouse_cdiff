# **Running Madi's 16S Sequencing Microbiome Profiling Data Analysis Workflow!**

## **Introduction**

I initially wrote this workflow to help myself out with 16S data analysis for a gut microbiome study I worked on for my graduate thesis.
Since then, it has proven helpful in subsequent 16S data analysis projects for myself and others, so I decided to put together a quick
tutorial on how to use it! 

This workflow can take raw 16S sequencing FASTA files and run them through QIIME2 mirobiome profiling software to return alpha and beta diversity measures and taxa barcharts along with generating basic plots and statistical results in R. Coming soon, this workflow will have an option to run PICRUSt2 analysis as well. 

*Disclaimer: While I am actively working on improving my workflow and making it more user-friendly, I reccommend that those who want to use it in it's current state have some amount of bioinformatics/software experience and be familiar with typical microbiome profiling analysis.* 

## **Workflow Options**

There are a few different places that you can start and stop my workflow depending on which parts of the 16S data analysis you need done:

1. You can run the analysis from start to finish, raw FASTA sequencing files to R plots and statistical results. 
2. You can come into the analysis with BIOM table and representative sequences .qza files. 
3. You can run QIIME2 analysis only without generating any plots or statistical results in R. 
4. For the R analysis, you have a few options depending on whether your samples were taken at multiple time points (longitudinal) or one time point (not longitudinal).

*Disclaimer: Parts of this workflow are extremely tailored to my current data analysis - I would not reccommend choosing the not longitudinal R option as the longitudinal R workflow can be tailored to a greater variety of analyses, even if they're not taken at multiple time points.* 

## **Tutorial**

If you've made it this far, congratulations, and I'm so sorry. In this tutorial, I will try my best to help you run my workflow on your 16S microbiome analysis!

### **Some brief organization that will make your life easier**

Based on how my workflow is designed and best practices, I'm offering some file system organization that will make your life easier. You can choose not to do this but let's just say, I warned you. I was also super nice and wrapped the file system organization in a `bash` script. If you run the following commands, it will setup your file system organization for you and download all needed files to run my workflow!

```bash
## getting said bash script off of my github
wget 
```

```bash
## creating a directory for the analysis
mkdir practice_workflow

## navigate into said directory
cd ./practice_workflow

## creating my_data and workflow subdirectories
mkdir my_data
mkdir workflow

## navigating into workflow subdirectory
cd ./workflow

## creating subdirectories in workflow 
mkdir config_files
mkdir envs
mkdir rules
```

So, after all those commands you should have a file system setup that looks something like this:

```bash
$ tree practice_workflow
practice_workflow
├── my_data
└── workflow
    ├── config_files
    │   └── config_template.yml
    ├── envs
    │   ├── install_envs_linux.sh
    │   ├── install_envs_macos.sh
    │   ├── picrust2.yml
    │   └── r_env.yml
    ├── rules
    │   ├── 01_demux_dada2.smk
    │   ├── 02_phylogeny.smk
    │   ├── 03_tss.smk
    │   ├── 04_core_metrics.smk
    │   ├── 05_longitudinal.smk
    │   └── 06_notLongitudinal.smk
    ├── run_snakemake.sh
    └── snakefile
```

**Next steps:**

1. Put the data that you want analyzed in the `my_data` directory. 


### **Installing snakemake and needed conda environments**

This workflow was written via Snakemake so you will first need to install it into a conda environment. I would reccommend creating a new conda environment for this installation so Snakemake is not installed into your base conda environment. 

```bash
## activating base conda environment if not activated already
conda activate base

## creating a conda environment named "snakemake_env" and directly installing snakemake into it
conda create -c conda-forge -c bioconda --name snakemake_env snakemake
```

If installing Snakemake via conda-forge and bioconda doesn't work out for you, you can install it via pip as well. Be warned, the Snakemake installation via pip doesn't have full capabilites (but will still work for running this workflow). 

```bash
## activating base conda environment if not activated already
conda activate base

## creating a new conda environment named "snakemake_env"
conda create --name snakemake_env

## activating the "snakemake_env" you just created
conda activate snakemake_env

## installing snakemake into this environment via pip
pip install snakemake
```

Since QIIME2 and R are used in this workflow, you will need to install separate conda environments for them prior to running the workflow. Luckily for you, I have written `.yaml` files and a `bash` script for the conda environment installation which are under `workflow/envs`. QIIME2 has different installation instructions based on the OS of your local computer; Linux users will run the `install_envs_linux.sh` script and MacOS users will run the `install_envs_macos.sh` script. 

```bash
## linux users
sh practice_workflow/workflow/envs/install_envs_linux.sh

## macos (apple silicon/arm64) users
sh practice_workflow/workflow/envs/install_envs_macos.sh
```

I currently do not have a `bash` script put together for Windows or non-Apple Silicon MacOS users but QIIME2 installation instructions for those operating systems can be found [here](https://docs.qiime2.org/2024.5/install/native/#install-qiime-2-within-a-conda-environment). 

Hopefully you have now successfully installed all needed conda environments to run the workflow! To double check, run the following code:

```bash
conda env list
```

Where you should see `snakemake_env`, `qiime2-2023.5`, and `r_env` listed among your other conda environments.




