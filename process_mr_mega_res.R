library(tidyverse)
args = commandArgs(trailingOnly=T)

setwd("/data/scratch/hmy117/gwas_raw_results/mr_mega/")

# read in codex
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()


i=as.numeric(args[1])
  trait = filelist$gh_phenotype[i]

  # read in results
  res = read_tsv(paste0(trait,"_mr_mega_res.result"),col_types = cols_only(
  MarkerName = col_character(),
  Chromosome = col_double(),
  Position = col_double(),
  EA = col_character(),
  NEA = col_character(),
  EAF = col_double(),
  Nsample = col_double(),
  Ncohort = col_double(),
  beta_0 = col_double(),
  se_0 = col_double(),
  beta_1 = col_double(),
  se_1 = col_double(),
  beta_2 = col_double(),
  se_2 = col_double(),
  beta_3 = col_double(),
  se_3 = col_double(),
  beta_4 = col_double(),
  se_4 = col_double(),
  chisq_association = col_double(),
  ndf_association = col_double(),
  `P-value_association` = col_double(),
  chisq_ancestry_het = col_double(),
  ndf_ancestry_het = col_double(),
  `P-value_ancestry_het` = col_double(),
  chisq_residual_het = col_double(),
  ndf_residual_het = col_double(),
  `P-value_residual_het` = col_double(),
  lnBF = col_double(),
  Comments = col_character()
  ))

  # fix P values
  res$`P.value_association` = pchisq(res$chisq_association, res$ndf_association, lower.tail=FALSE)
  res$`P.value_ancestry_het` = pchisq(res$chisq_ancestry_het, res$ndf_ancestry_het, lower.tail=FALSE)
  res$`P.value_residual_het` = pchisq(res$chisq_residual_het, res$ndf_residual_het, lower.tail=FALSE)

  message("fixed p values")

  # filter to SNPs tested in all cohorts

res = res %>%
  filter(is.na(Comments)) %>%
  mutate(p_het = `P.value_ancestry_het`) %>%
  filter(Ncohort == 7) %>%
  filter(EAF >= 0.01 & EAF <= 0.99) %>%
  filter(`P.value_residual_het`>0.05) %>%
  filter(Chromosome != "23")



# save
outfile = paste0(trait,"_mr_mega_res.tsv")
write_tsv(res,outfile)

