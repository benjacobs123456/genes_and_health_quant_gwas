#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 1

# get trait ID
trait=$(awk -F"," -v x=${SGE_TASK_ID} 'NR==x+1{print $2}' ~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv)

# define locus for LD ref

# unzip gwas
cd /data/scratch/hmy117/gwas_raw_results/susiex_inputs/

# harmonise gwas sum stats
module load R/4.2.2
Rscript /data/home/hmy117/gh_quant_traits/scripts/susie_harmonisation.R $trait

# light filtering of ref data
for ancestry in eur afr eas;
  do
    ANC=${ancestry^^}
    ~/plink2 --bfile ../kg_ref/g1000_$ancestry\_hg38_chrpos \
    --extract UKB_$trait\_$ANC\.tsv \
    --make-bed \
    --out g1000_$ancestry\_hg38_chrpos_filtered_$trait
  done 

~/plink2 --bfile ../kg_ref/g1000_sas\_hg38_chrpos \
--extract $trait\.tsv \
--make-bed \
--out g1000_sas_hg38_chrpos_filtered_$trait
  

# Harmonise with ref data
Rscript /data/home/hmy117/gh_quant_traits/scripts/remove_dups.R g1000_eur_hg38_chrpos_filtered_$trait\.bim UKB_Albumin_EUR.tsv EUR 
Rscript /data/home/hmy117/gh_quant_traits/scripts/remove_dups.R g1000_eas_hg38_chrpos_filtered_$trait\.bim UKB_Albumin_EAS.tsv EAS 
Rscript /data/home/hmy117/gh_quant_traits/scripts/remove_dups.R g1000_afr_hg38_chrpos_filtered_$trait\.bim UKB_Albumin_AFR.tsv AFR 
Rscript /data/home/hmy117/gh_quant_traits/scripts/remove_dups.R g1000_sas_hg38_chrpos_filtered_$trait\.bim Albumin.tsv SAS 


# more filtering to harmonised SNPs
for ancestry in eas afr eur sas;
  do
    ANC=$(echo ${ancestry~~})
    ~/plink2 --bfile g1000_$ancestry\_hg38_chrpos_filtered_$trait \
    --extract snps_to_include_for_susiex_$ANC\.tsv \
    --make-bed --out g1000_$ancestry\_hg38_chrpos_filtered2_$trait

     rm g1000_$ancestry\_hg38_chrpos_filtered_$trait\*
  done

# run SUSIE

# grab sig hits for this trait 

awk -v x=$trait -F "," 'NR==1{print $1,$2,$3,$9,$10,$15};NR>1{if($9<5e-8 && $10 <5e-8 && $15==x) print $1,$2,$3,$9,$10,$15}' "/data/home/hmy117/gh_quant_traits/outputs/all_meta_sig_hits_pre_vep.tsv" > susiex_loci_$trait
# get nrows 
