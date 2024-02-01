#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29

module load R/4.2.2
Rscript ~/gh_quant_traits/scripts/prepare_gh_files.R ${SGE_TASK_ID}
echo "finished"
