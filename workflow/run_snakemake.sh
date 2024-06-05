#!/bin/bash

## my command to run my snakemake workflow bc I'm lazy
snakemake -s workflow/snakefile -c 7 --use-conda --keep-going --configfile workflow/newExp_comp_d15_config.yml