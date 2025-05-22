#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=newExp_qiime
#SBATCH --output=myscratch/diet_cdiff_newExp/new_experiments/logs/ne_qiime_%J.log
#SBATCH --error=myscratch/diet_cdiff_newExp/new_experiments/logs/ne_qiime_%J.err
#SBATCH --mail-user=madison.apgar@cuanschutz.edu
#SBATCH --mail-type=FAIL,END

cd myscratch/diet_cdiff_newExp/

module load mambaforge

snakemake -s workflow/snakefile -c 20 --use-conda