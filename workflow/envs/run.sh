#!/bin/bash

set -e
set -u
set -o

echo "--------creating qiime environment"
conda update conda
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
CONDA_SUBDIR=osx-64 conda env create -n qiime2_testing --file qiime2-2023.5-py38-osx-conda.yml
conda config --env --set subdir osx-64
rm qiime2-2023.5-py38-osx-conda.yml
pip install snakemake

echo "--------creating R envrionment"
conda env create -f workflow/envs/r_env.yml

echo "--------creating picrust environment"
conda env create -f workflow/envs/picrust2.yml