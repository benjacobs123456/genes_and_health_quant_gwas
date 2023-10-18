#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 29

# get trait ID
module load R/4.2.2
Rscript ~/gh_quant_traits/scripts/process_mr_mega_res.R ${SGE_TASK_ID}

