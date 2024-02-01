#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:29
#$ -tc 29

# get trait ID
trait=$(awk -F"," -v x=${SGE_TASK_ID} 'NR==x+1{print $2}' ~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv)

# define locus for LD ref
# MUST HAVE RUN SUSIE_PREP.SH BEFOREHAND

# unzip gwas
cd /data/scratch/hmy117/gwas_raw_results/susiex_inputs/
# run SUSIE

# get nrows
nrow=$(awk 'END{print NR}' susiex_loci_$trait)

# run susiex
module load plink/1.9-170906

for (( i=2; i<=$nrow; i++ ))
  do
    # clean up
    rm $trait\_afr_ld_*
    rm $trait\_eas_ld_*
    rm $trait\_eur_ld_*
    rm $trait\_sas_ld_*

    locus_name=$(awk -v i=$i 'NR==i{print $6}' susiex_loci_$trait)
    locus_start=$(awk -v i=$i 'NR==i{print $3-5e5}' susiex_loci_$trait)
    locus_end=$(awk -v i=$i 'NR==i{print $3+5e5}' susiex_loci_$trait)
    locus_chr=$(awk -v i=$i 'NR==i{print $2}' susiex_loci_$trait)
    n1=$(awk 'NR==2{print $6}' $trait\.tsv)
    n2=$(awk 'NR==2{print $9}' UKB_$trait\_EUR.tsv)
    n3=$(awk 'NR==2{print $9}' UKB_$trait\_AFR.tsv)
    n4=$(awk 'NR==2{print $9}' UKB_$trait\_EAS.tsv)
    n5=$(awk 'NR==2{print $9}' UKB_$trait\_CSA.tsv)

    locus_name=$(echo $locus_name\_$i)
    echo "doing $locus_name"

    # run susiex
  /data/home/hmy117/SuSiEx/bin_static/SuSiEx \
  --sst_file=$trait\.tsv,UKB_$trait\_EUR.tsv,UKB_$trait\_AFR.tsv,UKB_$trait\_EAS.tsv,UKB_$trait\_CSA.tsv \
  --n_gwas=$n1,$n2,$n3,$n4,$n5 \
  --ref_file=g1000_sas_hg38_chrpos_filtered2_$trait\,g1000_eur_hg38_chrpos_filtered2_$trait\,g1000_afr_hg38_chrpos_filtered2_$trait\,g1000_eas_hg38_chrpos_filtered2_$trait\,g1000_sas_hg38_chrpos_filtered2_$trait \
  --ld_file=$trait\_sas_ld_filtered2,$trait\_eur_ld_filtered2,$trait\_afr_ld_filtered2,$trait\_eas_ld_filtered2,$trait\_sas_ld_filtered2 \
  --out_dir=/data/scratch/hmy117/gwas_raw_results/susie_res \
  --out_name=$locus_name \
  --chr=$locus_chr \
  --bp=$locus_start,$locus_end \
  --chr_col=1,1,1,1,1 \
  --snp_col=12,6,6,6,6 \
  --bp_col=2,2,2,2,2 \
  --a1_col=4,3,3,3,3 \
  --a2_col=3,4,4,4,4 \
  --eff_col=7,7,7,7,7 \
  --se_col=8,8,8,8,8 \
  --pval_col=9,10,10,10,10 \
  --plink=plink \
  --mult-step=T \
  --threads=$NSLOTS \
  --n_sig=10 \
  --maf 0.001 \
  --pval_thresh=1.72e-9 \
  --max_iter=100
done
