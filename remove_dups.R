library(tidyverse)

# get args
args = commandArgs(trailingOnly = T)
bim_path = args[1]
gwas_path = args[2]

# read in 
bim_file = read_table2(bim_path,col_names=F,col_types=cols(.default="c"))
gwas_file = read_table2(gwas_path,col_types = cols(.default="c"))

# combine 
bim_file = bim_file %>% 
  dplyr::select(2,5,6) %>%
  dplyr::rename("MARKERNAME" = X2)

gwas_file = gwas_file %>%
  filter(MARKERNAME %in% bim_file$MARKERNAME)

snps_to_keep = gwas_file %>%
  dplyr::select(MARKERNAME,EA,NEA) %>%
  left_join(bim_file,by="MARKERNAME") %>%
  filter(
    (EA == X5 & NEA == X6 ) | 
    (EA == X6 & NEA == X5 ) 
    ) %>%
  dplyr::select(MARKERNAME)

# get rid of dup positions 
snps_to_keep = snps_to_keep %>%
  dplyr::count(MARKERNAME) %>%
  filter(n==1) %>%
  dplyr::select(MARKERNAME)
  
outfile = paste0("snps_to_include_for_susiex_",args[3],".tsv")
write_tsv(snps_to_keep,outfile,col_names=F)  

