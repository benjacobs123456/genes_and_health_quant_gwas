library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/mr_mega/")

args = commandArgs(trailingOnly = T)
# read in codex
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()

# loop through each trait
i = as.numeric(args[1])

trait = filelist$gh_phenotype[i]

# find inputs
files_to_read = list.files("../mr_mega_inputs",
           full.names = T,
           pattern = trait)
files_to_read = files_to_read[files_to_read!=paste0("../mr_mega_inputs/UKB_",trait,".tsv")]
# write input
input_path = paste0(trait,"_mr_mega_input")
output_path = paste0(trait,"_mr_mega_res")

write_tsv(data.frame(files_to_read),input_path,col_names = F)

# run MR MEGA
cmd = paste0("~/MR-MEGA -i ",input_path," --qt -o ",output_path," --pc 4 --gc")
system(cmd)
