library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/")

args = commandArgs(trailingOnly = T)

# read in codex 
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()

# read in UKB N 
ukb_n = read_csv("~/gh_quant_traits/inputs/ukb_gwas_numbers.csv") %>%
  dplyr::select(1:7) %>%
  na.omit()

# combine
filelist = filelist %>% 
  left_join(ukb_n,by="ukb_phenotype")

i = as.numeric(args[1]) 
# get path
this_file = filelist[i,]
# get path 
path_to_gwas = str_remove_all(paste0("/data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/",this_file$ukb_path),".bgz")
message("Doing trait: ",this_file$ukb_phenotype)

# read in 
ukb_gwas = read_tsv(path_to_gwas,col_types = cols_only(
chr = col_character(),
pos = col_double(),
ref = col_character(),
alt = col_character(),
af_meta = col_double(),
beta_meta = col_double(),
se_meta = col_double(),
neglog10_pval_meta = col_double(),
af_AFR = col_double(),
af_AMR = col_double(),
af_CSA = col_double(),
af_EAS = col_double(),
af_EUR = col_double(),
af_MID = col_double(),
beta_AFR = col_double(),
beta_AMR = col_double(),
beta_CSA = col_double(),
beta_EAS = col_double(),
beta_EUR = col_double(),
beta_MID = col_double(),
se_AFR = col_double(),
se_AMR = col_double(),
se_CSA = col_double(),
se_EAS = col_double(),
se_EUR = col_double(),
se_MID = col_double(),
neglog10_pval_AFR = col_double(),
neglog10_pval_AMR = col_double(),
neglog10_pval_CSA = col_double(),
neglog10_pval_EAS = col_double(),
neglog10_pval_EUR = col_double(),
neglog10_pval_MID = col_double(),
low_confidence_AFR = col_logical(),
low_confidence_AMR = col_logical(),
low_confidence_CSA = col_logical(),
low_confidence_EAS = col_logical(),
low_confidence_EUR = col_logical(),
low_confidence_MID = col_logical()
))

# filter to maf > 0.01 
# filter to high confidence in all cohorts 
ukb_gwas = ukb_gwas %>%
  filter(af_meta >= 0.01 & af_meta <= 0.99) %>% 
  filter_at(vars(contains("low_confidence_")),all_vars(. == F))

# liftover to hg38  
hg38_coords = read_table("/data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg38.bed",
col_names = F,
col_types = cols(.default = "c"))

ukb_gwas = ukb_gwas %>%
  mutate(snp = paste0(chr,":",pos)) %>%
  filter(snp %in% hg38_coords$X4) 

hg38_coords = hg38_coords %>% 
dplyr::rename("snp"= X4,"pos"=X3) %>% 
  filter(snp %in% ukb_gwas$snp) %>%
  dplyr::select(snp,pos)

# unique positions 
hg38_coords = hg38_coords %>% distinct(snp,.keep_all=T)

# join 
ukb_gwas = ukb_gwas %>%
  dplyr::select(-pos) %>%
  left_join(hg38_coords,by="snp")

# split by ancestry
ancestries = c("AFR","AMR","CSA","EAS","EUR","MID")

for(ancestry in ancestries){
  # format for MR-MEGA
  this_ancestry_n = this_file[[ancestry]]
  this_ancestry = ukb_gwas %>% 
    dplyr::select(1,pos,2,3,4,contains(ancestry))
  colnames(this_ancestry) = str_remove_all(colnames(this_ancestry),paste0("_",ancestry))

  this_ancestry = this_ancestry %>%
    mutate(MARKERNAME = paste0(chr,":",pos), N = this_ancestry_n) %>%
    dplyr::rename(
      "POSITION" = pos,
      "EA" = alt,
      "EAF" = af,
      "NEA" = ref,
      "CHROMOSOME" = chr,
      "BETA" = beta,
      "SE" = se)

  this_ancestry = this_ancestry %>% 
    dplyr::select(CHROMOSOME,POSITION,EA,NEA,EAF,MARKERNAME,BETA,SE,N,neglog10_pval) %>%
    mutate(P = 10^(-1 * neglog10_pval) ) %>%
    dplyr::select(-neglog10_pval)
  
  this_ancestry = this_ancestry %>% mutate(cohort = "UKB")
  path_out = paste0("/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/UKB_",this_file$gh_phenotype,"_",ancestry,".tsv")
  write_tsv(this_ancestry,path_out)
  }

