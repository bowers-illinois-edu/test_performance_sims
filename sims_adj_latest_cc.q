#!/bin/bash

#SBATCH --job-name=sims_adj_latest_dt
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=jwbowers@illinois.edu
#SBATCH --export=PATH,R_LIBS,CLUSTER,CORES,LD_LIBRARY_PATH,MACHINE
#SBATCH --partition=IllinoisComputes
#
conda activate mambaR


## These next for pulling into R if needd
export R_DATATABLE_NUM_THREADS=1
# export CLUSTER="cc-socket"
export CORES=128
export MACHINE="CampusCluster"

echo "Started on: $(date --rfc-3339=seconds)"
echo "Hostname: $(hostname)"
echo "Working directory: $PWD"

#Rscript r-multihost.R
Rscript -e "remotes::install_github('bowers-illinois-edu/TreeTestSim')"
Rscript --verbose Simple_Analysis/local_adj_dt_sims.R
# R --file=Analysis/siusims.R


echo "Finished on: $(date --rfc-3339=seconds)"

