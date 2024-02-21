#!/bin/bash
#SBATCH -J LPS_cape
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=8G # memory pool for all cores
#SBATCH -t 0-48:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-user=anna.tyler@jax.org
#SBATCH --mail-type=END
#for now set arguments in Rmd file
# example use: sbatch -q long LPS_cape

cd $SLURM_SUBMIT_DIR

module load singularity

singularity exec ../../../Containers/R.sif R -e 'rmarkdown::render(here::here("Documents", "LPS_Sensitivity.Rmd"))'
