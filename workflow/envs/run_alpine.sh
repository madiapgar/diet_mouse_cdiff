#!/bin/bash

set -e
set -u
set -o

echo "--------creating qiime environment"
## conda update conda
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-linux-conda.yml
rm qiime2-2023.5-py38-linux-conda.yml
conda activate qiime2-2023.5
pip install snakemake
conda deactivate

echo "--------creating R envrionment"
mamba env create -f r_env.yml

echo "--------creating picrust environment"
mamba env create -f picrust2.yml
