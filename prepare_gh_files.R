library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/raw")

# read in codex
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()

args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

# loop through each file and process for MR-MEGA
# get path
this_file = filelist[i,]
message("Doing trait ",this_file$ukb_phenotype)
# get path
path_to_gwas = this_file$gh_path
message("Doing trait: ",this_file$gh_phenotype)
# read in
dat = read_table(path_to_gwas,
               col_types = cols_only(
  CHROM = col_double(),
  GENPOS = col_double(),
  ID = col_character(),
  ALLELE0 = col_character(),
  ALLELE1 = col_character(),
  A1FREQ = col_double(),
  N = col_double(),
  BETA = col_double(),
  SE = col_double(),
  LOG10P = col_double(),
)) %>%
  mutate(P = 10^-LOG10P,
         trait = this_file$gh_phenotype) %>%
  filter(A1FREQ >= 0.01 & A1FREQ <= 0.99) %>%
  dplyr::select(-LOG10P) %>%
  mutate(cohort = "GH") %>%
  mutate(MARKERNAME = paste0(CHROM,":",GENPOS)) %>%
  dplyr::rename(
    "POSITION" = GENPOS,
    "EA" = ALLELE1,
    "EAF" = A1FREQ,
    "NEA" = ALLELE0,
    "CHROMOSOME" = CHROM
    ) %>%
    dplyr::select(-ID)

path_out = paste0("/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/",this_file$gh_phenotype,".tsv")
write_tsv(dat,path_out)
