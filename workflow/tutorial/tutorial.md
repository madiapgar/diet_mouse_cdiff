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
wget https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/tutorial/setup_file_structure.sh

## running said bash script
sh setup_file_structure.sh
```

So, after running `setup_file_structure.sh` you should have a file system setup that looks something like this:

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

### **Setting up your config file**

If you navigate to your `practice_workflow/workflow/config_files` directory, you'll notice that I already put a `config_template.yml` file there. Let's open `config_template.yml` and take a look at it. 

```yaml
## config file template
## if you need additional examples for how to handle your config file, check out the workflow/config_files directory!

## the directory the data you're analyzing lives in (this should NOT be in the workflow directory) - in this case, it's your my_data directory
## you only have to name this once so all your subdirectories for raw sequencing files and metadata can be written as if you're already in my_data
dataset_dir: "my_data/"

## which QIIME2 environment did you install? 
qiime_env: "qiime2-2023.5"

## options: yes and no
raw_sequences: "are you starting the analysis with raw 16S fasta files?"

## 01: demux and dada2
raw_seq_dir: "the subdirectory(s) that your raw 16S fasta files live in"
raw_seqs: "a list of file and/or directory names/prefixes of your raw sequences (if you have more than one) - these are like wildcards 
    ex: ["seq1_run/seq1_", "seq2_run/seq2_", "seq3_run/seq3_"]"
barcodes: "this will most likely be the same as your raw_seqs entry, I just had my barcode files named differently"
## where do you want to trim and truncate your sequences during DADA2?
## pro tip: make sure these values are the same for every set sequencing experiments that you want to compare!
dada2_trim_left_for: 0
dada2_trim_left_rev: 0
dada2_trunc_len_for: 0
dada2_trunc_len_rev: 0

## 02/03: phylogeny and total sum scaling
r_script_dir: "the subdirectory that holds the r scripts you want snakemake to reference"
## if you ran raw sequences, your biom_table and rep_seqs will be under the file path/names below
## if not, you need to tell snakemake the file path to your biom_table and rep_seqs
biom_table: "data/qiime/merged_table.qza"
rep_seqs: "data/qiime/merged_rep_seqs.qza"

## 04: core metrics analysis
metadata: "what is your QIIME2-approved metadata file called?"
## what sampling depth do you want to use for your core metrics analysis? if you're not sure, I'd consult data/qiime/taxonomy_filtered.qzv
core_metrics_sampling_depth: 0

## 05: longitudinal dataset plots/stats in R 
## options: yes, no, or no downstream
longitudinal_dataset: "was data collected at multiple time points and do you have a column for that in your metadata? do you want R analysis run on your samples?"
## if you're running the cecal sample/non-longitudinal analysis, it will generate a processed metadata file for you so you just need to say what you want it to be called
processed_metadata: "what is/what do you want your processed metadata file called? - metadata doesn't need to be processed, but I use processed metadata files"

## 06: non-longitudinal dataset plots/stats in R 
## the non-longitudinal dataset delim is mainly for my cecal sample dataset and is extremely tailored to it rn
## I had to do some funky metadata wrangling and that's what these files are for
sampleID_key: NA
mouseID_facil_key: NA
## raw data files for desired quantified substances/inflammatory measures
## these get processed all nicely for you as well
bile_acid: NA
toxin: NA
histo: NA
metab: NA
hypoxia: NA
```
Once you have your config file set up the way you want it, you only have one more major thing to set up before you can run my workflow. 

### **Prepping your R scripts**

Unfortunately, data visualizations and statistical analyses are typically incredibly unique to whichever study you're doing. Because of this, I have provided some R script templates [here](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/templates) that you can take and make your own. I would only suggest making sure that your inputs/outputs match up with the `argparse` arguments at the top of the script or else my workflow may get mad at you. 

I would also suggest putting all of your R scripts in their own subdirectory under `my_data` (I like to just call it `scripts`) and update your config file accordingly. 

### **Running the actual workflow (finally!!)**

So, after all of your hard work, it's time to attempt to run my workflow. Be warned, you may have to do some debugging but once you get it working, it runs beautifully (much like GitHub). 

I've also included a handy `bash` script under `practice_workflow/workflow/` named `run_snakemake.sh`. This basically allows you to freely edit the snakemake command, which could mean switching out your config file, altering the amount of cores on your computer snakemake uses, and adding flags for the workflow as you see fit. So, let's take a look at `run_snakemake.sh`. 

```bash
snakemake \
    -s workflow/snakefile \ ## points to where the snakefile is 
    -c 7 \ ## the amount of cores snakemake will use
    --use-conda \ ## tells snakemake to use your conda environments 
    --keep-going \ ## tells snakemake to keep going if rules fail
    --configfile workflow/config_files/config.yml ## points to where your config file is
    
## if you would like to dry run your workflow, add --dry-run to the above command
```
Once you're satisified with your snakemake command in `run_snakemake.sh`, navigate into `practice_workflow` and activate your `snakemake_env` conda environment. From there you can run:

```bash
sh workflow/run_snakemake.sh
```

Happy analyzing! :) 




