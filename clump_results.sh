#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 29

# get trait ID
trait=$(awk -F"," -v x=${SGE_TASK_ID} 'NR==x+1{print $2}' ~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv)

echo doing $trait
~/plink --bfile /data/scratch/hmy117/gwas_raw_results/kg_ref/g1000_eur_hg38_qc \
--clump /data/scratch/hmy117/gwas_raw_results/mr_mega/$trait\_mr_mega_res.tsv \
--clump-p1 1.72e-9 --clump-r2 0.001 \
--out /data/scratch/hmy117/gwas_raw_results/clump_results/$trait \
--clump-snp-field MarkerName \
--clump-field P.value_association \
--clump-kb 1000 

# clean 
awk '{print $3}' /data/scratch/hmy117/gwas_raw_results/clump_results/$trait\.clumped > /data/scratch/hmy117/gwas_raw_results/sig_lead_snps_$trait
