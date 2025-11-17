#!/bin/sh
#SBATCH --job-name=DMR #Job name
#SBATCH --mail-type=NONE
#SBATCH --mail-user=gmiao@usf.edu # Where to send mail
#SBATCH --error=batch_%j.err
#SBATCH --output=batch_%j.out
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=5Gb
#SBATCH --account=group
#SBATCH --qos=group
#SBATCH --time=100:00:00       
#SBATCH --array=1
date;hostname;pwd

module load combined-pvalues
#https://github.com/brentp/combined-pvalues/tree/master/examples

comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 -p DMR/DMR dt_DMR.bed
