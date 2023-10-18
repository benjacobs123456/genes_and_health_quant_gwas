# Preamble
This repo contains code used to produce the results of "Genetic architecture of routinely acquired blood tests in a UK cohort of South Asian ancestry reveals ancestry-specific causal variants"

Contact: b.jacobs@qmul.ac.uk
18-10-23

# Code

## Download GH GWAS data
```unix
cd /data/scratch/hmy117/gwas_raw_results/raw/
~/google-cloud-sdk/bin/gsutil cp gs://genesandhealth_publicdatasets/results_GWAS_public/2023_07_version007phenotypes/2023_09_quant_trait_gwas_paper/* ./

# unzip
wait
for file in *.gz;
  do
    gunzip $file &
  done
```

## Download pan-UKB GWAS
```unix
cd /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30610-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30620-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30650-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30670-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30680-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30690-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30700-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30710-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30730-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30740-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30750-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30760-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30780-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30810-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30840-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30860-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30870-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30880-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30890-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30150-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30030-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30020-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30120-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30050-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30040-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30130-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30140-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30080-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30010-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30070-both_sexes-irnt.tsv.bgz &
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30000-both_sexes-irnt.tsv.bgz &

wait
for file in *.bgz;
  do
    bgzip -d $file &
  done
```

## Process GH GWAS sum stats
```r
library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/raw")

# read in codex
filelist = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv") %>%
  dplyr::select(1:4) %>%
  na.omit()

# loop through each file and process for MR-MEGA
for(i in c(1:nrow(filelist))){

  # get path
  this_file = filelist[i,]
  message("Doing trait ",this_file$ukb_phenotype)
  # get path
  path_to_gwas = paste0("/data/scratch/hmy117/gwas_raw_results/raw/",this_file$gh_path)
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
    INFO = col_double(),
    N = col_double(),
    BETA = col_double(),
    SE = col_double(),
    LOG10P = col_double(),
  )) %>%
    mutate(P = 10^-LOG10P,
           trait = this_file$gh_phenotype) %>%
    filter(A1FREQ >= 0.01 & A1FREQ <= 0.99 & INFO > 0.7) %>%
    dplyr::select(-INFO,-LOG10P) %>%
    mutate(cohort = "GH") %>%
    dplyr::rename(
      "MARKERNAME" = ID,
      "POSITION" = GENPOS,
      "EA" = ALLELE1,
      "EAF" = A1FREQ,
      "NEA" = ALLELE0,
      "CHROMOSOME" = CHROM
      )

  path_out = paste0("/data/scratch/hmy117/gwas_raw_results/mr_mega_inputs/",this_file$gh_phenotype,".tsv")
  write_tsv(dat,path_out)
}

```

## Liftover UKB GWAS
```unix

# liftover ukb gwas
awk 'NR>1{print $1,$2-1,$2,$1":"$2}' /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/biomarkers-30600-both_sexes-irnt.tsv > /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg19

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg19 \
/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
ukb_gwas_hg38.bed unlifted.bed

```


## Format UKB files
```unix
qsub ~/gh_quant_traits/scripts/prepare_ukb_files.sh
```

## Run MR-MEGA
```unix
qsub ~/gh_quant_traits/scripts/mr_mega.sh
```


## Plot ancestral PCs
```r
library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/")

# read in log to check PCs
mr_mega_log = read_table("./mr_mega/Albumin_mr_mega_res.log",skip=50,n_max = 7)
mr_mega_log$PCs = str_remove_all(mr_mega_log$PCs,"../mr_mega_inputs/")
mr_mega_log$PCs = str_remove_all(mr_mega_log$PCs,".tsv")

a = ggplot(mr_mega_log,aes(PC0,PC1,label=PCs))+geom_point()+ggrepel::geom_text_repel(show.legend = F)+
  theme_minimal()+
    scale_color_brewer(palette="Paired")
b = ggplot(mr_mega_log,aes(PC1,PC2,label=PCs))+geom_point()+ggrepel::geom_text_repel(show.legend = F)+
  theme_minimal()+
  scale_color_brewer(palette="Paired")
png("pc_plots.png",res=900,units="in",width=8,height=4)
gridExtra::grid.arrange(a,b,ncol=2)
dev.off()
```

## Prepare 1kg files for clumping
```unix
cd /data/scratch/hmy117/gwas_raw_results/kg_ref

# download ref files
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip &
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_sas.zip &
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_amr.zip &
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_afr.zip &
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip &

wait
# unzip
unzip -o g1000_eur.zip &
unzip -o g1000_sas.zip &
unzip -o g1000_amr.zip &
unzip  -o g1000_eas.zip &
unzip  -o g1000_afr.zip &

# liftover
awk '{print "chr"$1,$4-1,$4,$2}' g1000_eur.bim > g1000_eur_hg19_bed
awk '{print "chr"$1,$4-1,$4,$2}' g1000_sas.bim > g1000_sas_hg19_bed
awk '{print "chr"$1,$4-1,$4,$2}' g1000_afr.bim > g1000_afr_hg19_bed
awk '{print "chr"$1,$4-1,$4,$2}' g1000_eas.bim > g1000_eas_hg19_bed
awk '{print "chr"$1,$4-1,$4,$2}' g1000_amr.bim > g1000_amr_hg19_bed

# liftOver
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_eur_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_eur_hg38_bed unlifted.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_sas_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_sas_hg38_bed unlifted2.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_eas_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_eas_hg38_bed unlifted3.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_afr_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_afr_hg38_bed unlifted4.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_amr_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_amr_hg38_bed unlifted5.bed &

# update SNP IDs
awk '{print $4,$3}' g1000_eur_hg38_bed > eur_vars
awk '{print $4,$3}' g1000_sas_hg38_bed > sas_vars
awk '{print $4,$3}' g1000_eas_hg38_bed > eas_vars
awk '{print $4,$3}' g1000_afr_hg38_bed > afr_vars
awk '{print $4,$3}' g1000_amr_hg38_bed > amr_vars

# update files to hg38
for ancestry in eur sas eas amr afr;
  do
    ~/plink2 --bfile g1000_$ancestry \
    --update-map $ancestry\_vars \
    --make-pgen \
    --out g1000_$ancestry\_hg38 \
    --sort-vars &
  done

for ancestry in eur sas eas amr afr;
  do
    ~/plink2 --pfile g1000_$ancestry\_hg38 \
    --set-all-var-ids @:# \
    --rm-dup exclude-all \
    --out g1000_$ancestry\_hg38_chrpos \
    --make-bed
  done


# light QC
for ancestry in eur sas eas amr afr;
  do
    ~/plink2 --bfile g1000_$ancestry\_hg38_chrpos \
    --maf 0.01 \
    --hwe 1e-5 \
    --geno 0.1 \
    --make-bed \
    --out g1000_$ancestry\_hg38_qc &
  done

```

## Process MR-MEGA results
```unix
qsub ~/gh_quant_traits/scripts/process_mr_mega_res.sh
```

## Clump results
```unix
qsub ~/gh_quant_traits/scripts/clump_results.sh
qsub ~/gh_quant_traits/scripts/clump_gh_results.sh
```

## Plots
```unix
qsub ~/gh_quant_traits/scripts/make_manhattans.sh
```

## Combine results for summary tables
```unix
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = read_csv(files,col_types = "ddcccdddddcccc")

vep_input = dat %>%
mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
dplyr::select(1,start,end,alleles,strand)
write_tsv(vep_input,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv",col_names=F)

# repeat for meta
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = read_csv(files,col_types = "ddcccdddddcccc")

vep_input = dat %>%
mutate(start = Position,end = Position,alleles = paste0(EA,"/",NEA),strand="+") %>%
dplyr::select(1,start,end,alleles,strand)
write_tsv(vep_input,"/data/home/hmy117/gh_quant_traits/outputs/all_meta_sig_hits_for_vep.tsv",col_names=F)

```

## VEP
```unix
module load perl
module load perl5lib
cd /data/home/hmy117/gh_quant_traits/outputs/
~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv \
--cache /data/scratch/hmy117/.vep \
--force_overwrite \
--most_severe \
--tab --fields "Uploaded_variation,Consequence"

cd /data/home/hmy117/gh_quant_traits/outputs/
~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv \
--cache /data/scratch/hmy117/.vep \
--force_overwrite \
--everything \
-o gh_sig_hits_all_annotations.tsv

cd /data/home/hmy117/gh_quant_traits/outputs/
~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
-o nearest_gene.tsv \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Feature,Feature_type,Consequence"

 ~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_meta_sig_hits_for_vep.tsv --cache --dir_cache /data/scratch/hmy117/.vep --force_overwrite --nearest symbol -o meta_nearest_gene.tsv --tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Feature,Feature_type,Consequence"
```


## Combine results for summary tables (with VEP results)
```r
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annot",files)]

# read VEP results
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/nearest_gene.tsv",skip=32,col_types = cols(.default="c"))
vep_res2 = read_table("/data/home/hmy117/gh_quant_traits/outputs/variant_effect_output.txt",skip=26,col_types = cols(.default="c"))
vep_res3 = read_table("/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_all_annotations.tsv",skip=105,col_types = cols(.default="c"))

# take one gene per SNP
vep_res = vep_res %>% distinct(`#Uploaded_variation`,NEAREST)
vep_res2 = vep_res2 %>% distinct(`#Uploaded_variation`,.keep_all=T)


# read in gh gwas
dat = read_csv(files,col_types = "ddcccdddddcccc") %>%
  mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(CHROMOSOME,"_",POSITION,"_",EA,"/",NEA)) %>%
  left_join(vep_res %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  dplyr::select(1,2,full_snp_id,NEA,EA,EAF,N,BETA,SE,P,trait,sig_snp,gh_specific_locus,NEAREST) %>%
  left_join(vep_res2 %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id")

write_csv(dat,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_mapped.csv")

# read in gh gwas
dat2 = read_csv(files,col_types = "ddcccdddddcccc") %>%
  mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(CHROMOSOME,"_",POSITION,"_",EA,"/",NEA)) %>%
  left_join(vep_res3 %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  mutate(MARKERNAME = full_snp_id)

write_csv(dat2,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_mapped_full_anno.csv")

# sum plots
n_loci = dat %>% group_by(trait) %>%
dplyr::count(sig_snp) %>%
na.omit() %>%
  left_join(
    dat %>%
      group_by(trait) %>%
      summarise(n = max(N)),
      by="trait")

      library(ggrepel)
p=ggplot(n_loci,aes(n.y,n.x,col=trait,label=trait))+
  geom_point()+
  theme_bw()+
  geom_text_repel(color="black",min.segment.length=0)+
  labs(x="N GWAS", y="N loci")+
  theme(legend.position="none")


png("/data/home/hmy117/gh_quant_traits/outputs/plots/n_loci_vs_n_gwas.png",
  res=900,units="in",height=5,width=6)
print(p)
dev.off()

# repeat for meta

files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = read_csv(files,col_types = "ddcccdddddcccc")

# read VEP results
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/meta_nearest_gene.tsv",skip=32,col_types = cols(.default="c"))

# take one gene per SNP
vep_res = vep_res %>% distinct(`#Uploaded_variation`,NEAREST)


# read in gh gwas
dat = dat %>%
  mutate(start = Position,end = Position,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(Chromosome,"_",Position,"_",EA,"/",NEA)) %>%
  left_join(vep_res %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  mutate(MarkerName = full_snp_id)

write_tsv(dat,"/data/home/hmy117/gh_quant_traits/outputs/all_meta_anno.tsv",col_names=F)
```
## Run SuSiEx
```unix
qsub ./scripts/susie_prep.sh
qsub ./scripts/susie.sh
```

## Read in fine mapping results
```r
library(tidyverse)
# get files
files = list.files("/data/scratch/hmy117/gwas_raw_results/susie_res",full.names = T,pattern="summary")

# read in all credible sets
finemap_res = purrr::map(files,function(x){
  df = read_table(x,skip = 1,col_types = cols(.default = "c"))
  locus_name =
    colnames(
      read_table(x,n_max =1,col_types = cols(.default = "c"))
    )[2]
  df$locus_name = locus_name
  trait = str_split_fixed(str_split_fixed(x,pattern = "res\\/",n=2)[1,2],pattern="_",n=2)[1,1]
  df$trait = trait
  df
})
finemap_res = do.call("bind_rows",finemap_res)

# explore data
counts = finemap_res %>%
  distinct(trait,locus_name) %>%
  dplyr::count(trait)

ggplot(counts,aes(n,trait))+
  geom_col()

# get median credible set
finemap_res %>%
  filter(CS_LENGTH==1)


# SAS-specific
sas_gh_specific = finemap_res %>%
  filter(`POST-HOC_PROB_POP1` > 0.7 & `POST-HOC_PROB_POP5` > 0.7) %>%
  filter_at(vars(
             c(`POST-HOC_PROB_POP2`,
               `POST-HOC_PROB_POP3`,
               `POST-HOC_PROB_POP4`
               )),
             all_vars( . < 0.7 )

             )

# save to file
write_csv(sas_gh_specific,"~/gh_quant_traits/outputs/sas_specific_finemap.csv")

# filter further to those that replicate
write_csv(finemap_res,"~/gh_quant_traits/outputs/all_finemap.csv")

# locus plot

locus = "16:88649138"
snp_files = files = list.files("/data/scratch/hmy117/gwas_raw_results/susie_res",full.names = T,pattern=".snp")

res = read_tsv("/data/scratch/hmy117/gwas_raw_results/susie_res/Vitamin_D_2.snp",col_types = cols(.default="c")) %>%
dplyr::select(1,2,contains("PIP")) %>%
pivot_longer(cols = contains("PIP")) %>%
mutate(pos = as.numeric(BP),pip = as.numeric(value))


ggplot(res,aes(pos,pip,col=name))+
geom_point()+
theme_bw()

```

## Locus plots
```unix
qsub ./make_locus_plots.sh
```

## Pleiotropy plot GH
```r
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = read_csv(files,col_types = "ddcccdddddcccc")

# pleio plot

# variant
snp = "16:88716656"

# forest plot

plot_dat = dat %>%
  filter(MARKERNAME == snp)

# align to effect allele
ea1 = plot_dat$EA[1]
plot_dat = plot_dat %>%
  mutate(BETA = ifelse(EA==ea1,BETA,BETA*-1)) %>%
  mutate(EAF = ifelse(EA==ea1,EAF,1-EAF))

# forest
plot_dat = plot_dat %>%
  arrange(BETA)
plot_dat$trait = factor(plot_dat$trait,levels = plot_dat$trait,ordered=T)



  p=ggplot(plot_dat,aes(BETA,trait,size=-log10(P),color=trait))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.5)+
    geom_errorbarh(mapping= aes(y=trait,xmin = BETA - 1.96*SE,xmax = BETA+1.96*SE),color="black",height=0.1,linewidth=0.1)+
  geom_point()+
  theme_bw()+
  scale_color_brewer(palette="Paired")+
  labs(size = "-log10(P)",color="Trait",
       x = paste0("Effect size (beta)\n of ",snp,"-",ea1," on trait"))
  varname = str_replace_all(snp,":","_")
  png(paste0("/data/home/hmy117/gh_quant_traits/outputs/pleio_plots/GH_",varname,".png"),res=900,units="in",width=6,height=4)
  print(p)
  dev.off()

  # gwas catalogue
  # get GWAS catalogue
gwas_cat = read_tsv("/data/scratch/hmy117/full") %>%
  mutate(snp_id = paste0(CHR_ID,":",CHR_POS)) %>%
  filter(snp_id == snp)

# loop through
overall_res = list()
for(i in c(1:nrow(meta_res))){
  message(i)
  this_snp = meta_res[i,]

  # filter gwas_cat
  gwas_cat_assocs = gwas_cat %>%
    filter(CHR_ID == this_snp$Chromosome) %>%
    filter(as.numeric(CHR_POS) > (as.numeric(this_snp$Position) - 5e5)) %>%
    filter(as.numeric(CHR_POS) < (as.numeric(this_snp$Position) + 5e5)) %>%
    dplyr::count(`DISEASE/TRAIT` )

  this_snp$associations = list(gwas_cat_assocs$`DISEASE/TRAIT`)
  overall_res[[i]] = this_snp
}
overall_res = do.call("bind_rows",overall_res)


```

## Phewas
```r
library(tidyverse)
library(ggrepel)
setwd("/data/home/hmy117/gh_quant_traits/")

trait = commandArgs(trailingOnly = TRUE)[1]


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
), col_types = "ddcccdddddc")

# combine UKB_EUR and GH
combined_gwas = bind_rows(
  gh_sas %>% mutate(CHROMOSOME = as.character(CHROMOSOME),ancestry="SAS"),
  ukb_eur %>% mutate(ancestry = "EUR"),
  ukb_sas %>% mutate(ancestry = "SAS"),
  ukb_eas %>% mutate(ancestry = "EAS"),
  ukb_afr %>% mutate(ancestry = "AFR"))

# variant
snp = "16:88716656"
# forest plot

plot_dat = combined_gwas %>%
  filter(MARKERNAME == snp)

# align to effect allele
ea1 = plot_dat$EA[1]
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
    geom_errorbarh(mapping= aes(y=cohort,xmin = BETA - 1.96*SE,xmax = BETA+1.96*SE),color="black",height=0.1,linewidth=0.1)+
  geom_point()+
  theme_bw()+
  scale_color_brewer(palette="Paired")+
  labs(size = "EAF",color="Ancestry",
       x = paste0("Effect size (beta)\n of ",snp,"-",ea1," on ",trait))
  varname = str_replace_all(snp,":","_")
  png(paste0("/data/home/hmy117/gh_quant_traits/outputs/forest_plots/",trait,"_",varname,".png"),res=900,units="in",width=6,height=4)
  print(p)
  dev.off()
```
