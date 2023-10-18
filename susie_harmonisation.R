library(tidyverse)

args = commandArgs(trailingOnly=T)
trait = args[1]

# gwas files
gwas_files = list.files("../mr_mega_inputs/",pattern=trait,full.names=T)

# read in
gwas_res = purrr::map(gwas_files,function(x){
read_table(x,col_types = cols_only(
  MARKERNAME = col_character(),
  EA = col_character(),
  NEA = col_character())) %>%
  mutate(cohort = x)
})

# combine
gwas_res_long = do.call("bind_rows",gwas_res)

# filter to SNPs in all cohorts
filtered_snps = gwas_res_long %>%
  mutate(full_id = paste0(MARKERNAME,"_",NEA,"_",EA)) %>%
  dplyr::count(full_id) %>%
  filter(n==7)

# filter to these SNPs
gwas_res = purrr::map(gwas_files,function(x){
dat = read_table(x,col_types = cols(.default = "c")) %>%
  mutate(full_id = paste0(MARKERNAME,"_",NEA,"_",EA)) %>%
  filter(full_id %in% filtered_snps$full_id) %>%
  dplyr::select(-full_id)
  y = str_remove_all(x,"..\\/mr_mega_inputs/")
  y = str_remove_all(y,"\\/")
  message("writing to ",y)
  write_tsv(dat,y)
})
