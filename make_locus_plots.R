library(tidyverse)

args = commandArgs(trailingOnly = T)
trait = args[1]

# read in GWAS files
files = list.files("/data/scratch/hmy117/gwas_raw_results/susiex_inputs/",pattern=trait,full.names=T)
files = files[grepl("tsv",files)]

# read in susiex loci
susiex_loci = read_table(
  paste0(
    "/data/scratch/hmy117/gwas_raw_results/susiex_inputs/susiex_loci_",
    trait
    ), col_types = "cddddc"
  )

# define locus
for(i in c(1:nrow(susiex_loci))){
  locus_chr = susiex_loci$Chromosome[i]
  locus_bp = susiex_loci$Position[i]


  # read in and filter
  res_sas_file = paste0("/data/scratch/hmy117/gwas_raw_results/susiex_inputs/",trait,".tsv")
  res_sas = read_tsv(res_sas_file, col_types = cols_only(
      CHROMOSOME = col_double(),
      POSITION = col_double(),
      MARKERNAME = col_character(),
      P = col_double(),
      cohort = col_character())) %>%
      filter(
        CHROMOSOME == locus_chr &
        POSITION > locus_bp - 5e5 &
        POSITION < locus_bp + 5e5 )

    res_eur_file = paste0("/data/scratch/hmy117/gwas_raw_results/susiex_inputs/UKB_",trait,"_EUR.tsv")
  res_eur = read_tsv(res_eur_file, col_types = cols_only(
      CHROMOSOME = col_double(),
      POSITION = col_double(),
      MARKERNAME = col_character(),
      P = col_double(),
      cohort = col_character())) %>%
      filter(
        CHROMOSOME == locus_chr &
        POSITION > locus_bp - 5e5 &
        POSITION < locus_bp + 5e5 )


combo_df = bind_rows(
  res_eur %>% mutate(ancestry="EUR"),
  res_sas %>% mutate(ancestry="SAS")
  )



  # locus plot

snp_file = paste0(
  "/data/scratch/hmy117/gwas_raw_results/susie_res/",
  trait,"_",i+1,".snp"
  )

# read in fine mapping results
res = read_tsv(snp_file,col_types = cols(.default="c")) %>%
dplyr::select(1,2,contains("PIP")) %>%
pivot_longer(cols = contains("PIP")) %>%
mutate(pos = as.numeric(BP),pip = as.numeric(value))

# filter GWAS to fine mapped SNPs
combo_df = combo_df %>% filter(MARKERNAME %in% res$SNP)

library(ggsci)
library(Hmisc)
############
# sas
############

# get lead SNP in SAS
snps_for_ld = combo_df %>% distinct(MARKERNAME)
write_tsv(snps_for_ld,"/data/scratch/hmy117/ld_snps_locus_tmp.tsv",col_names=F)

sas_bim = read_table("/data/scratch/hmy117/gwas_raw_results/kg_ref/g1000_sas_hg38_chrpos.bim",col_names=F,col_types = cols_only(
  X2 = col_character()
  ))

sas_lead = combo_df %>%
  filter(ancestry=="SAS") %>%
  filter(MARKERNAME %in% sas_bim$X2) %>%
  slice_min(P,with_ties=F)

# SAS LD
cmd = paste0(
  "~/plink --bfile ",
  "/data/scratch/hmy117/gwas_raw_results/susiex_inputs/g1000_sas_hg38_chrpos_filtered2_",trait," ",
  "--ld-snp ",sas_lead$MARKERNAME," --r2 --ld-window-r2 0.000001 --ld-window 99999 ",
  "--out /data/home/hmy117/ld_snps_for_plot_sas"
  )
  system(cmd)

sas_ld = read_table("/data/home/hmy117/ld_snps_for_plot_sas.ld",col_types = cols(.default="c"))

sas_dat_for_plot = combo_df %>%
  filter(ancestry=="SAS") %>%
  left_join(sas_ld %>%
      dplyr::rename("MARKERNAME" = SNP_B),by="MARKERNAME") %>%
      mutate(r2 = as.numeric(R2))

sas_dat_for_plot$r2_bin = Hmisc::cut2(sas_dat_for_plot$r2,cuts = c(0,0.2,0.4,0.6,0.8,1),drop=F,oneval=F,onlycuts=T)
sas_dat_for_plot$r2_bin = fct_rev(sas_dat_for_plot$r2_bin)

sas_plot = ggplot(sas_dat_for_plot %>% filter(r2 > 0.1),
  aes(POSITION,-log10(P),fill=r2_bin))+
  geom_point(shape=21,color="black")+
  theme_bw()+
  ggtitle("SAS")+
  scale_fill_locuszoom(drop=F)+
  theme(legend.position="none")+
  geom_point(data = sas_dat_for_plot %>% filter(r2 <= 0.1),shape=21,color="black",fill="grey",alpha=0.1)



# EUR
eur_bim = read_table("/data/scratch/hmy117/gwas_raw_results/kg_ref/g1000_eur_hg38_chrpos.bim",col_names=F,col_types = cols_only(
  X2 = col_character()
  ))

eur_lead = combo_df %>%
  filter(ancestry=="EUR") %>%
  filter(MARKERNAME %in% eur_bim$X2) %>%
  slice_min(P,with_ties=F)

# eur LD
cmd = paste0(
  "~/plink --bfile ",
  "/data/scratch/hmy117/gwas_raw_results/susiex_inputs/g1000_eur_hg38_chrpos_filtered2_",trait," ",
  "--ld-snp ",eur_lead$MARKERNAME," --r2 --ld-window-r2 0.000001 --ld-window 99999 ",
  "--out /data/home/hmy117/ld_snps_for_plot_eur"
  )
  system(cmd)

eur_ld = read_table("/data/home/hmy117/ld_snps_for_plot_eur.ld",col_types = cols(.default="c"))

eur_dat_for_plot = combo_df %>%
  filter(ancestry=="EUR") %>%
  left_join(eur_ld %>%
      dplyr::rename("MARKERNAME" = SNP_B),by="MARKERNAME") %>%
      mutate(r2 = as.numeric(R2))

eur_dat_for_plot$r2_bin = Hmisc::cut2(eur_dat_for_plot$r2,cuts = c(0,0.2,0.4,0.6,0.8,1),drop=F,oneval=F,onlycuts=T)
eur_dat_for_plot$r2_bin = fct_rev(eur_dat_for_plot$r2_bin)

eur_plot = ggplot(eur_dat_for_plot %>% filter(r2 > 0.1),
  aes(POSITION,-log10(P),fill=r2_bin))+
  geom_point(shape=21,color="black")+
  theme_bw()+
  ggtitle("EUR")+
  scale_fill_locuszoom(drop=F)+
  theme(legend.position="none")+
  geom_point(data = eur_dat_for_plot %>% filter(r2 <= 0.1),shape=21,color="black",fill="grey",alpha=0.1)


# plot
finemap = ggplot(res %>% filter(pip > 0.1),aes(pos,pip,fill=name))+
geom_point(shape=21,color="black",size=3)+
theme_bw()+
scale_fill_brewer(palette="Paired")+
theme(legend.position="none")+
labs(y="PIP",x="POSITION")+
ggtitle("Multi-ancestry fine mapping")+
geom_point(data = res %>% filter(pip <= 0.1),shape=21,color="black",fill="grey",alpha=0.1)


outfile = paste0("/data/home/hmy117/gh_quant_traits/outputs/locus_plots/",trait,"_locus_chr_",locus_chr,"_locus_bp_",locus_bp,".png")

png(outfile,res=900,units="in",height=8,width=6)
print(
  cowplot::plot_grid(eur_plot,sas_plot,finemap,align = "v",ncol=1,rel_heights = c(1/3,1/3,1/3))
  )
dev.off()
}
