#!/bin/bash

#SBATCH --job-name=sims_adj_latest_dt
#SBATCH --ntasks-per-node=100
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=jwbowers@illinois.edu
#SBATCH --export=PATH,R_LIBS,CLUSTER,CORES,LD_LIBRARY_PATH
#SBATCH --partition=IllinoisComputes
#
#module load gnu/gnu-9.3.0 gnu/openmpi-3.1.6-gnu-9.3.0
#source ~/miniforge3/bin/activate
conda activate mambaR

## set InputDir=/data/keeling/a/jwbowers/Documents/PROJECTS/manytests-paper
## set RunDir=/data/keeling/a/jwbowers/Documents/PROJECTS/manytests-paper

## cd /data/keeling/a/jwbowers/Documents/PROJECTS/manytests-paper

## These next for pulling into R if needd
export R_DATATABLE_NUM_THREADS=1
# export CLUSTER="cc-socket"
export CORES=100

echo "Started on: $(date --rfc-3339=seconds)"
echo "Hostname: $(hostname)"
echo "Working directory: $PWD"

#Rscript r-multihost.R
Rscript -e "remotes::install_github('bowers-illinois-edu/TreeTestSim')"
Rscript --verbose Simple_Analysis/local_adj_dt_sims.R
# R --file=Analysis/siusims.R


echo "Finished on: $(date --rfc-3339=seconds)"

