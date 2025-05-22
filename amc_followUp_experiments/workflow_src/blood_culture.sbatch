#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=blood_culture_qiime
#SBATCH --output=myscratch/diet_cdiff_bloodCulture/logs/bc_qiime_%J.log
#SBATCH --error=myscratch/diet_cdiff_bloodCulture/logs/bc_qiime_%J.err
#SBATCH --mail-user=madison.apgar@cuanschutz.edu
#SBATCH --mail-type=FAIL,END

cd myscratch/diet_cdiff_bloodCulture/

module load mambaforge

snakemake -s workflow/blood_culture_snakefile.smk -c 20 --use-conda