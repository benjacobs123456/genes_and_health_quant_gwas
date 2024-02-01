#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 29

source ~/popcorn/bin/activate

# get trait ID
trait=$(awk -F"," -v x=${SGE_TASK_ID} 'NR==x+1{print $2}' ~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv)

# define locus for LD ref

# unzip gwas
cd /data/scratch/hmy117/

# calculate cross-ancestry LD 
awk 'BEGIN {OFS="\t"}; NR==1{print "SNP","A1","A2","N","beta","SE"}; NR>1{print $12,$4,$3,$6,$7,$8}' /data/scratch/hmy117/gwas_raw_results/susiex_inputs/$trait\.tsv > /data/scratch/hmy117/gwas_raw_results/susiex_inputs/popcorn_$trait\.tsv
awk 'BEGIN {OFS="\t"}; NR==1{print "SNP","A1","A2","N","beta","SE"}; NR>1{print $6,$3,$4,$9,$7,$8}' /data/scratch/hmy117/gwas_raw_results/susiex_inputs/UKB_$trait\_EUR.tsv > /data/scratch/hmy117/gwas_raw_results/susiex_inputs/popcorn_UKB_$trait\_EUR.tsv

popcorn fit -v 1 \
--maf 0.05 \
--cfile ~/gh_quant_traits/outputs/sas_vs_eur_scores.txt \
--sfile1 /data/scratch/hmy117/gwas_raw_results/susiex_inputs/popcorn_$trait\.tsv \
--sfile2 /data/scratch/hmy117/gwas_raw_results/susiex_inputs/popcorn_UKB_$trait\_EUR.tsv \
~/gh_quant_traits/outputs/cross_ancestry_correlation_$trait\.txt
