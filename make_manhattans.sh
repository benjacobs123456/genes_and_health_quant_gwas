#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 5

module load R/4.2.2
Rscript /data/home/hmy117/gh_quant_traits/scripts/make_manhattans.R ${SGE_TASK_ID}

