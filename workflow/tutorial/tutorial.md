# **Running Madi's 16S Sequencing Microbiome Profiling Data Analysis Workflow!**

## **Introduction**

I initially wrote this workflow to help myself out with 16S data analysis for a gut microbiome study I worked on for my graduate thesis. Since then, it has proven helpful in subsequent 16S data analysis projects for myself and others, so I decided to put together a quick tutorial on how to use it!   

This workflow can take raw 16S sequencing FASTQ files and run them through QIIME2 microbiome profiling software to return alpha and beta diversity measures and taxa barcharts along with generating basic plots and statistical results in R (see disclaimer under **Prepping your R scripts**). While users can just directly use QIIME2 for their 16S analysis, I've found this workflow incredibly helpful when analyzing multiple sequencing runs and in the case of human error, which I'm extremely prone to, rerunning the analysis. 

*Disclaimer: While I am actively working on improving my workflow and making it more user-friendly, I reccommend that those who want to use it in it's current state have some amount of bioinformatics/software experience and be familiar with typical microbiome profiling analysis.* 

## **Workflow Options**

There are a few different places that you can start and stop my workflow depending on which parts of the 16S data analysis you need done:

1. You can run the analysis from start to finish, raw FASTQ sequencing files through alpha/beta diversity.
    - via all global options in the config file being set to **"yes"**
2. You can run Demux and DADA2 separately so you can take the BIOM table and representative sequences and insert them elsewhere.
    - via `raw_sequences:`**"yes"** in the config file
    - this is really helpful when you have multiple sequencing runs since this step automatically combines the BIOM tables and representative sequences for each run into one file
3. You can come into the analysis with BIOM table and representative sequences .qza files and run taxonomic classification/core metrics analysis from there.
    - via `raw_sequences`:**"no"** in the config file
    - via `tax_class` and `core_metrics`:**"yes"** in the config file
    - make sure that you include the file paths to your BIOM table and representative sequences under `biom_table`/`rep_seqs` in the config file
4. You can replicate what was done in the project featured in this github repository, where total sum scaling was used to "rarefy" the BIOM table due to an inert signal of *Lactococcus* contamination.
    - via `total_sum_scaling`:**"yes"** in the config file
    - make sure that you specify the directory containing `total_sum_scaling.R` under `r_script_dir` in the config file
5. You can stop the workflow after taxonomic classification to determine your sampling depth for core metrics analysis. 
    - via `tax_class`:**"yes"** and `core_metrics`:**"no"** in the config file 
    - once you know your sampling depth and have updated `core_metrics_sampling_depth` in your config file, you can change `core_metrics:` **"yes"** and rerun the workflow 
6. You can generate basic plots/stats in R.
    - via `microbiome_r_outputs`:**"yes"** in the config file
    - make sure that you specify the directory with you R scripts based on the templates provided (see disclaimer under **Prepping your R scripts**) under `r_script_dir` in the config file

If this seems all a little complicated, don't worry, we're going over the config file in more depth later. 

## **Tutorial**

In this tutorial, I will try my best to help you run my workflow on your 16S microbiome analysis (or the sequences associated with this project that you pulled off of QIITA)!

### **Some brief organization that will make your life easier**

Based on how my workflow is designed and best practices, I'm offering some file system organization that will make your life easier. I wrote a handy `bash` [script](https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/tutorial/setup_file_structure.sh) that will set up your file system organization for you and download all needed files to run my workflow!

```bash
## getting said bash script off of my github
wget https://github.com/madiapgar/diet_mouse_cdiff/blob/master/workflow/tutorial/setup_file_structure.sh

## running said bash script
bash setup_file_structure.sh
```

So, after running `setup_file_structure.sh` you should have a file system setup that looks something like this:

```bash
$ tree madis_16s_workflow
madis_16s_workflow
└── workflow
    ├── config_files
    │   └── config_template.yml
    ├── envs
    │   ├── install_envs_linux.sh
    │   ├── install_envs_macos.sh
    │   └── r_env.yml
    ├── rules
    │   ├── 01_demux_dada2.smk
    │   ├── 02_phylogeny.smk
    │   ├── 03_tss.smk
    │   ├── 04_core_metrics.smk
    │   ├── 05_microbiome_r.smk
    │   └── 06_madis_cecalAnalysis.smk
    ├── run_snakemake.sh
    └── snakefile
```

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

Since QIIME2 and R are used in this workflow, you will need to install separate conda environments for them prior to running the workflow. Luckily for you, I have written `.yaml` files and a `bash` script for the conda environment installation which are under `workflow/envs`. QIIME2 has different installation instructions based on the OS of your local computer; Linux users will run the `install_envs_linux.sh` script and MacOS users will run the `install_envs_macos.sh` script. **If you already have QIIME2 installed on your computer, you can skip this step.**

```bash
## linux users
bash madis_16s_workflow/workflow/envs/install_envs_linux.sh

## macos (apple silicon/arm64) users
bash madis_16s_workflow/workflow/envs/install_envs_macos.sh
```

I currently do not have a `bash` script put together for Windows or non-Apple Silicon MacOS users but QIIME2 installation instructions for those operating systems can be found [here](https://docs.qiime2.org/2024.5/install/native/#install-qiime-2-within-a-conda-environment). 

Hopefully you have now successfully installed all needed conda environments to run the workflow! To double check, run the following:

```bash
conda env list
```

Where you should see `snakemake_env`, `qiime2-2023.5`, and `r_env` listed among your other conda environments.

### **Setting up your config file**

If you navigate to your `madis_16s_workflow/workflow/config_files` directory, you'll notice that I already put a `config_template.yml` file there. Let's open `config_template.yml` and take a look at it. 

*A note: Like any other software tool, my workflow is something that may need to be run multiple times with differing parameters in order to get to the end result you desire, it just depends on how complicated your data is. The config file is incredibly flexible and can be easily edited between runs to reflect new desired parameters or files so don't be afraid to switch things up!* 

```yaml
## config file template
## if you need additional examples for how to handle your config file, check out the workflow/config_files directory!

dataset_dir: "the directory the data you're analyzing lives in (this should NOT be in the workflow directory)"
qiime_env: "which QIIME2 environment did you install?"

## --GLOBAL OPTIONS--
## are you starting with raw 16S sequences?
raw_sequences: "no" ## options: yes and no 

## do you want to run taxonomic classification? chances are yes since you're running this workflow
tax_class: "yes" ## options: yes and no, default is yes 

## does your data need to go through total sum scaling? chances are that it doesn't
total_sum_scaling: "yes" ## options: yes and no

## do you want to run core metrics analysis on your data to get alpha/beta diversity? you HAVE to run this if you want the microbiome R outputs
core_metrics: "yes" ## options: yes and no

## would you like basic R plots/stats for your microbiome data? if so, this requires you to edit your R scripts directory for the exact visualizations you want
## (see more in the tutorial)
microbiome_r_outputs: "yes" ## options: yes and no
## -------------------

## --NEEDED FILE PATHS--
## include if raw_sequences = "yes"
raw_seq_dir: "the subdirectory(s) that your raw 16S fasta files live in"
raw_seqs: "a list of file and/or directory names/prefixes of your raw sequences (if you have more than one) - these are wildcards"
"ex: ["seq1_run/seq1_", "seq2_run/seq2_", "seq3_run/seq3_"]" ## assumes that everything after this is "_paired_end_seqs.qza"
barcodes: "this will most likely be the same as your raw_seqs entry, I just had my barcode files named differently" ## assumes that everything after this is "_barcodes.txt"
## where do you want to trim and truncate your sequences during DADA2?
## pro tip: make sure these values are the same for every set sequencing experiments that you want to compare!
dada2_trim_left_for: 0
dada2_trim_left_rev: 0
dada2_trunc_len_for: 0
dada2_trunc_len_rev: 0

## include if tax_class = "yes"
## if you ran raw sequences, your biom_table and rep_seqs will be under the file path/names below
## if not, you need to tell snakemake the file path to your biom_table and rep_seqs
biom_table: "data/qiime/merged_table.qza"
rep_seqs: "data/qiime/merged_rep_seqs.qza"
metadata: "what is your QIIME2-approved metadata file called?"

## include if total_sum_scaling = "yes"
r_script_dir: "the subdirectory that holds the r scripts you want snakemake to reference"

## include if core_metrics = "yes"
## what sampling depth do you want to use for your core metrics analysis? 
## if you're not sure, I'd consult data/qiime/taxonomy_filtered.qzv
core_metrics_sampling_depth: 0

## include if microbiome_r_outputs = "yes"
## r script directory if you haven't already (r_script_dir: "scripts/")
processed_metadata: "what is/what do you want your processed metadata file called? - metadata doesn't need to be processed, but I use 
processed metadata files"
## -------------------

## madis cecal dataset plots/stats in R (measurements at only ONE time point) - you will NOT use this 
## include if running that analysis 
madis_cecal_analysis: "no"
sampleID_key: NA
mouseID_facil_key: NA
bile_acid: NA
toxin: NA
histo: NA
metab: NA
hypoxia: NA
```
Once you have your config file set up the way you want it, you only have one more major thing to set up before you can run my workflow. 

### **Prepping your R scripts**

Unfortunately, data visualizations and statistical analyses are typically incredibly unique to whichever study you're doing. Because of this, I have provided some R script templates [here](https://github.com/madiapgar/diet_mouse_cdiff/tree/master/workflow/templates/r_script_templates/) that you can take and make your own. I would only suggest making sure that your inputs/outputs match up with the `argparse` arguments at the top of the script or else my workflow may get mad at you. 

I would also suggest putting all of your R scripts in their own subdirectory under `madis_16s_workflow` (I like to just call it `scripts` or `workflow_src`) and update your config file accordingly. 

### **Running the actual workflow (finally!!)**

So, after all of your hard work, it's time to attempt to run my workflow. Be warned, you may have to do some debugging but once you get it working, it runs beautifully (much like GitHub).

The most important part of this analysis are the raw 16S sequences so let's place the final piece of the puzzle. Let's create a subdirectory under `madis_16s_workflow/` to put the raw sequences in.

```bash
mkdir my_raw_16s_data
```

Copy your raw 16S sequences into your new `my_raw_16s_data` subdirectory. When you're done, you should have two subdirectories, `my_raw_16s_data` and `workflow` (which you received when you ran the `setup_file_structure.sh` script). Since the workflow doesn't take the completely raw 16S sequences due to them being either single or paired end, the 'raw' data the workflow is expecting a QIIME2 imported object and its associated `barcodes.txt` file. 

> [!IMPORTANT]
> My workflow currently **only supports** paired-end analysis!

```bash
## single end qiime import command
qiime tools import \
    --type EMPSingleEndSequences \
    --input-path raw_seqs \ ## this is a directory with the forward, reverse, and barcodes .fastq.gz files
    --output-path single_end_raw_seqs.qza

## paired end qiime import command
qiime tools import \
    --type EMPPairedEndSequences \
    --input-path raw_seqs \ ## this is a directory with the forward, reverse, and barcodes .fastq.gz files
    --output-path paired_end_raw_seqs.qza
``` 

I've also included a handy `bash` script under `madis_16s_workflow/workflow/` named `run_snakemake.sh`. This basically allows you to freely edit the snakemake command, which could mean switching out your config file, altering the amount of cores on your computer snakemake uses, and adding flags for the workflow as you see fit. So, let's take a look at `run_snakemake.sh`. 

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




