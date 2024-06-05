#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=comp_r
#SBATCH --output=myscratch/diet_mouse_cdiff_comparisons/logs/comp_r_%J.log
#SBATCH --error=myscratch/diet_mouse_cdiff_comparisons/logs/comp_r_%J.err
#SBATCH --mail-user=madison.apgar@cuanschutz.edu
#SBATCH --mail-type=FAIL,END

cd myscratch/diet_mouse_cdiff_comparisons/

module load mambaforge

snakemake -c 20 --use-conda
