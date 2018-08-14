#!/bin/bash
#SBATCH -A snic2017-1-573
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5-00:00:00
 
module load R/3.4.3

Rscript R/sim_f2.R
Rscript R/sim_f2_diverse.R
Rscript R/sim_f8.R
