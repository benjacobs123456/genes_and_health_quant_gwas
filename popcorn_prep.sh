#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=128G
#$ -l h_rt=240:0:0
#$ -l highmem
#$ -cwd
#$ -o /data/scratch/hmy117

# after setup run once, run this 
source ~/popcorn/bin/activate

# compute scores for SAS vs EUR 
popcorn compute -v 3 --bfile1 /data/scratch/hmy117/gwas_raw_results/kg_ref/g1000_sas_hg38_chrpos \
--bfile2 /data/scratch/hmy117/gwas_raw_results/kg_ref/g1000_eur_hg38_chrpos \
~/gh_quant_traits/outputs/sas_vs_eur_scores.txt \
--maf 0.05 \
--SNPs_to_store 300000
