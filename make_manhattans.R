library(tidyverse)
library(ggrepel)
setwd("/data/home/hmy117/gh_quant_traits/")

# read in codex
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()


# loop through each trait
args = commandArgs(trailingOnly=T)
i = as.numeric(args[1])

trait = filelist$gh_phenotype[i]

# read in results
res = read_tsv(paste0("/data/scratch/hmy117/gwas_raw_results/mr_mega/",trait,"_mr_mega_res.tsv"), col_types = cols_only(
MarkerName = col_character(),
Chromosome = col_double(),
Position = col_double(),
EA = col_character(),
NEA = col_character(),
EAF = col_double(),
Nsample = col_double(),
Ncohort = col_double(),
P.value_association = col_double(),
P.value_ancestry_het = col_double(),
P.value_residual_het = col_double())
)

# read in ukb data
ukb_eur = read_tsv(paste0(
  "/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/UKB_",trait,"_EUR.tsv")
)

ukb_afr = read_tsv(paste0(
  "/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/UKB_",trait,"_AFR.tsv")
)

ukb_eas = read_tsv(paste0(
  "/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/UKB_",trait,"_EAS.tsv")
)

ukb_sas = read_tsv(paste0(
  "/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/UKB_",trait,"_CSA.tsv")
)

# read in GH SAS
gh_sas = read_tsv(paste0(
  "/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/",trait,".tsv"
), col_types = "ddccdddddccc")

# read in GH clump res
sig_snps_gh = read_tsv(paste0("/data/scratch/hmy117/gwas_raw_results/gh_sas_sig_lead_snps_",trait),col_types = "c")

# read in meta clump res
sig_snps = read_tsv(paste0("/data/scratch/hmy117/gwas_raw_results/sig_lead_snps_",trait),col_types = "c")

# make manhattan

## get chrom midpoints
chr_coords = res %>%
group_by(Chromosome) %>%
summarise(min_bp = min(Position),
        max_bp = max(Position),
        midpoint = mean(Position)) %>%
ungroup %>%
mutate(cumbp = cumsum(max_bp)) %>%
mutate(Chromosome = Chromosome + 1)

res = res %>%
  left_join(chr_coords,by="Chromosome") %>%
  mutate(cum_pos = ifelse(is.na(cumbp),Position,Position + cumbp))

midpoints = chr_coords %>%
  mutate(midpoint_cum = cumbp - (max_bp - min_bp)/2 )

sig_hits = res %>%
filter(MarkerName %in% sig_snps$SNP)



# loop through each GH hit and see if there are known hits nearby
p_study_sig = 5e-8 / 29
sig_hits_gh = gh_sas %>% filter(P<p_study_sig)
specific_snps = data.frame()

# combine UKB GWAS
ukb_combined = bind_rows(ukb_eur,ukb_eas,ukb_sas,ukb_afr) %>%
filter(P < 1e-5) %>%
dplyr::select(CHROMOSOME,POSITION)

for(i in c(1:nrow(sig_hits_gh))){
message(i, "of ",nrow(sig_hits_gh))
this_snp = sig_hits_gh[i,]

  matches = ukb_combined %>%
    filter(CHROMOSOME == this_snp$CHROMOSOME &
             POSITION > this_snp$POSITION - 1e6 &
             POSITION < this_snp$POSITION + 1e6)

  if(nrow(matches)==0){
    specific_snps <<- bind_rows(specific_snps,this_snp)
  }

}

# filter to individual loci
## take lowest P value SNP within 1MB window
specific_loci = data.frame()

if(nrow(specific_snps)>0){
for(i in c(1:nrow(specific_snps))){
  message(i)
  this_snp = specific_snps[i,]
  more_sig_snps = specific_snps %>%
    filter(CHROMOSOME == this_snp$CHROMOSOME &
             POSITION > this_snp$POSITION - 1e6 &
             POSITION < this_snp$POSITION + 1e6 &
             P < this_snp$P)
  if(nrow(more_sig_snps)==0){
    specific_loci <<- specific_loci %>%
      bind_rows(this_snp)
  }
}
}

# forest plot

# combine UKB_EUR and GH
combined_gwas = bind_rows(
  gh_sas %>% mutate(CHROMOSOME = as.character(CHROMOSOME),ancestry="SAS"),
  ukb_eur %>% mutate(ancestry = "EUR"),
  ukb_sas %>% mutate(ancestry = "SAS"),
  ukb_eas %>% mutate(ancestry = "EAS"),
  ukb_afr %>% mutate(ancestry = "AFR"))

if(nrow(specific_snps)!=0){
  for(i in c(1:nrow(specific_snps))){
  this_snp = specific_snps[i,]

  # variant
  # forest plot
  plot_dat = combined_gwas %>%
    filter(MARKERNAME == this_snp$MARKERNAME)

  # align to effect allele
  ea1 = this_snp$EA
  plot_dat = plot_dat %>%
    mutate(BETA = ifelse(EA==ea1,BETA,BETA*-1)) %>%
    mutate(EAF = ifelse(EA==ea1,EAF,1-EAF))

  # forest
  plot_dat = plot_dat %>%
    arrange(BETA) %>%
    mutate(pop = paste0(cohort,"_",ancestry))
  plot_dat$cohort = factor(plot_dat$pop,levels = plot_dat$pop,ordered=T)



  p=ggplot(plot_dat,aes(BETA,cohort,size=EAF,color=ancestry))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.5)+
    geom_errorbarh(mapping= aes(y=cohort,xmin = BETA - 1.96*SE,xmax = BETA+1.96*SE),color="black",height=0.1,size=0.1)+
  geom_point()+
  theme_bw()+
  scale_color_brewer(palette="Paired")+
  labs(size = "EAF",color="Ancestry",
       x = paste0("Effect size (beta)\n of ",this_snp$MARKERNAME,"-",ea1," on ",trait))
  varname = str_replace_all(this_snp$MARKERNAME,":","_")
  png(paste0("/data/home/hmy117/gh_quant_traits/outputs/forest_plots/",trait,"_",varname,".png"),res=900,units="in",width=6,height=4)
  print(p)
  dev.off()
  }
}

# make plots
res = res %>%
mutate(color = case_when(
  `P.value_association` < p_study_sig & MarkerName %in% specific_snps$MARKERNAME ~ "gwas_sig_specific",
  `P.value_association` < p_study_sig ~ "gwas_sig",
  Chromosome %% 2 == 0 ~ "even",
  Chromosome %% 2 != 0 ~ "odd",
))

pal = c("orange1","red","lavenderblush1","lavenderblush2")
names(pal) = c("gwas_sig_specific","gwas_sig","even","odd")

# manhattan
n_loci = nrow(sig_snps)
n_spec = nrow(specific_loci)

p1=ggplot(res %>%
       filter(`P.value_association`< 0.05),
     aes(cum_pos,-log10(`P.value_association`),col = color))+
geom_point(size = 0.5)+
scale_color_manual(values = pal) +
theme_bw()+
scale_x_continuous(breaks = midpoints$midpoint_cum,
                   labels = c(1:22))+
theme(panel.grid = element_blank(),legend.position = "none",
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank())+
labs(y="-log10(P)")+
ggtitle(paste0("GWAS Meta-analysis of ",trait,"\nN loci: ",n_loci,"\nN specific loci: ",n_spec))


# het manhattan
sig_het_hits = res %>%
filter(`P.value_ancestry_het`< p_study_sig) %>%
mutate(rounded_cumpos = round(cum_pos/1e6,0)) %>%
group_by(rounded_cumpos) %>%
slice_min(`P.value_ancestry_het`) %>%
dplyr::select(1,2,3,`P.value_ancestry_het`) %>%
ungroup() %>%
dplyr::select(-rounded_cumpos)

# recolour
res = res %>%
mutate(color = case_when(
  `P.value_ancestry_het` < p_study_sig & `P.value_association`< p_study_sig  ~ "gwas_sig_het",
  Chromosome %% 2 == 0 ~ "even",
  Chromosome %% 2 != 0 ~ "odd",
))

pal = c("orange1","lavenderblush1","lavenderblush2")
names(pal) = c("gwas_sig_het","even","odd")


p2=ggplot(res %>%
       filter(`P.value_ancestry_het`< 0.05),
     aes(cum_pos,-log10(`P.value_ancestry_het`),col = color))+
geom_point(size = 0.5)+
scale_color_manual(values = pal) +
theme_bw()+
scale_x_continuous(breaks = midpoints$midpoint_cum,
                   labels = c(1:22))+
theme(panel.grid = element_blank(),legend.position = "none",
      axis.text.x = element_text(hjust=0,angle=90))+
labs(x="Chromosome",y="-log10(P-Het)") +
ggtitle("Ancestral heterogeneity")

# get results for writing
results_for_file = res %>%
  mutate(sig_snp = ifelse(MarkerName %in% sig_snps$SNP,"lead_sig_snp"," ")) %>%
  mutate(gh_specific_locus = ifelse(MarkerName %in% specific_snps$MARKERNAME,"specific_locus"," ")) %>%
  filter(  `P.value_association` < p_study_sig ) %>%
  dplyr::select(-min_bp,-max_bp,-midpoint,-cum_pos,-cumbp) %>%
  mutate(trait = trait)
outfile = paste0("/data/home/hmy117/gh_quant_traits/outputs/sig_hits_",trait,".csv")
write_csv(results_for_file,outfile)


# GH results for file
results_for_file_gh = gh_sas %>%
  mutate(sig_snp = ifelse(MARKERNAME %in% sig_snps_gh$SNP,"lead_sig_snp"," ")) %>%
  mutate(gh_specific_locus = ifelse(MARKERNAME %in% specific_snps$MARKERNAME,"specific_locus"," ")) %>%
  filter( P < p_study_sig )
outfile = paste0("/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_",trait,".csv")
write_csv(results_for_file_gh,outfile)

# combine with UKB & GH
res = res %>%
  left_join(ukb_eur %>%
    dplyr::select(MARKERNAME,P) %>%
    dplyr::rename("MarkerName" = MARKERNAME,
                  "UKB_P" = P),
                  by="MarkerName") %>%
  left_join(gh_sas %>%
    dplyr::select(MARKERNAME,P) %>%
    dplyr::rename("MarkerName" = MARKERNAME,
                  "GH_P" = P),
                  by="MarkerName")

# Make GH manhattan
# make plots
res = res %>%
mutate(color = case_when(
  GH_P < p_study_sig & MarkerName %in% specific_snps$MARKERNAME ~ "gwas_sig_specific",
  GH_P < p_study_sig ~ "gwas_sig",
  Chromosome %% 2 == 0 ~ "even",
  Chromosome %% 2 != 0 ~ "odd",
))

pal = c("orange1","red","lavenderblush1","lavenderblush2")
names(pal) = c("gwas_sig_specific","gwas_sig","even","odd")

# manhattan

p3=ggplot(res %>%
       filter(GH_P< 0.05),
     aes(cum_pos,-log10(GH_P),col = color))+
geom_point(size = 0.5)+
scale_color_manual(values = pal) +
theme_bw()+
scale_x_continuous(breaks = midpoints$midpoint_cum,
                   labels = c(1:22))+
theme(panel.grid = element_blank(),legend.position = "none",
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank())+
labs(y="-log10(P)")+
ggtitle("Genes & Health SAS GWAS")

# make UKB manhattan

# make plots
res = res %>%
mutate(color = case_when(
UKB_P < p_study_sig & MarkerName %in% specific_snps$MARKERNAME ~ "gwas_sig_specific",
UKB_P < p_study_sig ~ "gwas_sig",
Chromosome %% 2 == 0 ~ "even",
Chromosome %% 2 != 0 ~ "odd",
))

pal = c("orange1","red","lavenderblush1","lavenderblush2")
names(pal) = c("gwas_sig_specific","gwas_sig","even","odd")

# manhattan

p4=ggplot(res %>%
   filter(UKB_P< 0.05),
 aes(cum_pos,-log10(UKB_P),col = color))+
geom_point(size = 0.5)+
scale_color_manual(values = pal) +
theme_bw()+
scale_x_continuous(breaks = midpoints$midpoint_cum,
               labels = c(1:22))+
theme(panel.grid = element_blank(),legend.position = "none",
  axis.line.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank())+
labs(y="-log10(P)")+
ggtitle("UKB EUR GWAS")


png(paste0("./outputs/plots/manhattans_",trait,".png"),
  res=900,units="in",height=8,width=8)
print(cowplot::plot_grid(p1, p3,p4,p2, align = "v", nrow = 4, rel_heights = c(0.25,0.25,0.25,0.25)))
dev.off()
