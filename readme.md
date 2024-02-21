# Preamble
This repo contains code used to produce the results of "Genetic architecture of routinely acquired blood tests in a UK cohort of South Asian ancestry reveals ancestry-specific causal variants"

Link to preprint: [click here](https://www.researchsquare.com/article/rs-3438851/v1)

Contact: b.jacobs@qmul.ac.uk
12-02-23

# Code
This repo contains the following steps/sections:
1. [Step 0: GWAS and analysis within the Genes & Health TRE](#phenotype-preparation-and-gwas)
2. [Step 1: munge GWAS summary statistics](#process-gwas-data)
3. [Step 2: run & analyse MR-MEGA meta-analysis](#meta-analysis)
4. [Step 3: fine-mapping](#fine-mapping)
5. [Step 4: miscellaneous (pleiotropy, phewas, popcorn)](#misc)

## Phenotype preparation and GWAS
Note that this section was run inside the Genes & Health Google Cloud TRE.

### Phenotype prep
```R
library(tidyverse)
setwd("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/")

## Read in phenotype data
# get short names for phenos
short_names = list.files("/genesandhealth/red/quantitative_traits_v2/outputs/",
                         pattern = "stats",
                         full.names = F)

phenos = stringr::str_split_fixed(pattern = "_per", n=2, short_names)

# read in shortlist
shortlist = read_csv("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/inputs/trait_shortlist.csv")

# list files
file_to_read = paste0("/genesandhealth/red/quantitative_traits_v2/outputs/",shortlist$trait,
                      "_per_individual_stats_qc_2023-12-15.csv")

# read in files
results = readr::read_csv(file_to_read, show_col_types = FALSE)
input_file = read_csv("/genesandhealth/red/quantitative_traits_v2/inputs/trait_input_file.csv")

# read in unique timepoint(raw readings)
file_to_read = paste0("/genesandhealth/red/quantitative_traits_v2/outputs/",shortlist$trait,
                      "_readings_at_unique_timepoints_readings_2023-12-15.csv")

# read in files
raw_results = readr::read_csv(file_to_read, show_col_types = FALSE)

counts_per_source = raw_results %>% distinct(pseudo_nhs_number,trait,source) %>% dplyr::count(trait,source)

counts_all = raw_results %>% distinct(pseudo_nhs_number,trait) %>% dplyr::count(trait)  %>% mutate(source="both")
counts_all = counts_all %>% bind_rows(counts_per_source)
prop_primary_care = counts_all %>% pivot_wider(id_cols = trait,names_from = source,values_from = n) %>% mutate(prop = primary_care / both * 100) %>% arrange(prop)
prop_primary_care$trait = factor(prop_primary_care$trait,levels = prop_primary_care$trait,ordered=T)
ggplot(prop_primary_care,aes(prop,trait,fill=trait))+
  geom_col(color="black",show.legend = F)+
  theme_minimal()+
  labs(x="Proportion of individuals with at least 1 primary care reading")


counts = results %>%
  dplyr::count(trait) %>%
  arrange(desc(n))

counts2 = results %>%
  filter(n_values >= 2) %>%
  dplyr::count(trait) %>%
  arrange(desc(n))

n_readings = results %>%
  group_by(trait) %>%
  dplyr::summarise(n = sum(n_values)) %>%
  dplyr::rename("total_n_readings" = n)

counts_tbl = counts %>%
  left_join(counts2,by="trait") %>%
  left_join(n_readings,by="trait") %>%
  dplyr::rename("N" = n.x,
                "N with >=2 readings" = n.y)
counts_tbl = counts_tbl %>%
  left_join(input_file,by="trait")

full_trait_names = shortlist %>% dplyr::select(1,2)
counts_tbl = counts_tbl %>%
  left_join(full_trait_names,by="trait")


write_csv(counts_tbl %>%
            dplyr::select(1,9,2,3,4,5),
          "/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/quant_traits_counts.csv")

# sum stats
counts_tbl %>%
  summarise(median(N), range(N), median(total_n_readings), range(total_n_readings))

# make histogram
counts_tbl = counts_tbl %>%
  arrange(desc(N))

counts_tbl$trait = factor(counts_tbl$`Trait (full name)`,levels=counts_tbl$`Trait (full name)`,ordered = T)
p=ggplot(counts_tbl,aes(N,trait))+
  geom_point(size=3)+
  theme_minimal()+
  labs(x = "N individuals per trait",y="Trait")
p

plot_fx(p,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/counts")


# repeat for primary care only
pc_counts = raw_results %>%
  filter(source=="primary_care") %>%
  group_by(pseudo_nhs_number,trait,unit) %>%
  dplyr::count()

pc_counts = pc_counts %>%
  ungroup() %>%
  left_join(full_trait_names,by="trait")

# see what % of tests come from primary care

source_prop = raw_results %>%
  group_by(trait) %>%
  dplyr::count(source) %>%
  mutate(pct = 100* n/sum(n)) %>%
  arrange(desc(pct)) %>%
  ungroup()
source_prop = source_prop %>% left_join(full_trait_names,by="trait") %>% mutate(trait =`Trait (full name)`)
lvls = source_prop %>% filter(source=="primary_care")
source_prop$trait = factor(source_prop$trait,levels = unique(lvls$trait),ordered=T)
p = ggplot(source_prop,aes(trait,pct,fill=source,label=paste0(round(pct,1),"%")))+
  geom_col(position = position_stack(),col="black")+
  geom_text(data = source_prop %>% filter(source=="primary_care"),mapping = aes(trait,90))+
  scale_fill_brewer(palette="Set1",labels = c("Primary care","Secondary care"))+
  labs(x="Trait",y="% of readings from primary/secondary care",fill="Source of reading")+
  coord_flip()+
  theme_minimal()
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/source_of_reading.png",res=900,units="in",width=8,height=8)
p
dev.off()


# distribution of no tests
plot_dat = results %>%
  arrange(desc(n_values))

ggplot(plot_dat,aes(n_values,trait))+
  geom_boxplot(outlier.shape = NA)+
  labs(x="N readings per participant",y="Trait")+
  theme_minimal()+
  scale_x_continuous(limit = c(0,50))

plot_dat$n_values_bin = Hmisc::cut2(plot_dat$n_values,cuts=c(0,1,2,5,10,20,50,100,5000))
levels(plot_dat$n_values_bin) = c("0","1","2-5","5-10","10-20","20-50","50-100",">100")

plot_dat = plot_dat %>%
  group_by(trait) %>%
  dplyr::count(n_values_bin)

# get nice trait names
plot_dat = plot_dat %>% left_join(shortlist,by="trait")
p=ggplot(plot_dat,aes(n_values_bin,n,fill=`Trait (full name)`))+
  geom_col(color="black",alpha=0.8,show.legend = F)+
  facet_wrap(~`Trait (full name)`,ncol=4)+
  theme_minimal()+
  labs(x="N readings per participant",y="N")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/no_of_readings",res=900,units="in",width=12,height=12)
p
dev.off()

# relationship between multiple readings and values

results = results %>% left_join(shortlist,by="trait")
plots = list()
for(this_trait in unique(results$`Trait (full name)`)){
  plot_dat = results %>% filter(`Trait (full name)` == this_trait)

  # models

  spearman = cor.test(plot_dat$n_values,plot_dat$mean,method="spearman")
  graph_label = paste0(this_trait,"\nRho = ",round(spearman$estimate,2))
  plots[[length(plots)+1]] = ggplot(plot_dat,aes(n_values,mean))+
    geom_point()+
    ggtitle(graph_label)+
    theme_minimal()+
    labs(x="N values",y="Mean value")

}

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/no_of_readings_vs_value.png",res=900,units="in",width=16,height=12)
do.call("grid.arrange",plots)
dev.off()


# relationship between record length and mean value
plots = list()
for(this_trait in unique(results$`Trait (full name)`)){
  plot_dat = results %>% filter(`Trait (full name)` == this_trait)

  # models

  spearman = cor.test(plot_dat$time_from_earliest_to_latest,plot_dat$mean,method="spearman")
  graph_label = paste0(this_trait,"\nRho = ",round(spearman$estimate,2))
  plots[[length(plots)+1]] = ggplot(plot_dat,aes(time_from_earliest_to_latest,mean))+
    geom_point()+
    ggtitle(graph_label)+
    theme_minimal()+
    labs(x="N values",y="Mean value")

}

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/length_of_record_vs_value.png",res=900,units="in",width=24,height=12)
do.call("grid.arrange",plots)
dev.off()

# influence of N readings and length of record
res_df = data.frame()
for(this_trait in unique(results$`Trait (full name)`)){
  plot_dat = results %>% filter(`Trait (full name)` == this_trait)

  # normalise mean
  plot_dat$mean = RNOmni::RankNorm(plot_dat$mean)

    # model
  model = summary(lm(data = plot_dat, mean ~ n_values))

  out = data.frame(
    trait = this_trait,
    r2 = model$r.squared,
    beta = model$coefficients[2,1],
    se = model$coefficients[2,2]
  )
  res_df <<- bind_rows(res_df,out)
}

# plot
res_df = res_df %>%
  tibble() %>%
  arrange(desc(beta))
res_df$trait = factor(res_df$trait,levels = res_df$trait,ordered=T)
p = ggplot(res_df,aes(beta,trait))+
  geom_vline(xintercept = 0)+
  geom_errorbarh(mapping = aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se, y = trait),height=0.1)+
  geom_point(size=3)+
  theme_bw()+
  labs(x="Beta (relationship between N readings\nand mean value)",y="Trait")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/n_readings_vs_value_betas",res=900,units="in",width=6,height=6)
p
dev.off()

# influence of length of record
res_df = data.frame()
for(this_trait in unique(results$`Trait (full name)`)){
  plot_dat = results %>% filter(`Trait (full name)` == this_trait)

  # normalise mean
  plot_dat$mean = RNOmni::RankNorm(plot_dat$mean)

    # model
  model = summary(lm(data = plot_dat, mean ~ time_from_earliest_to_latest))

  out = data.frame(
    trait = this_trait,
    r2 = model$r.squared,
    beta = model$coefficients[2,1],
    se = model$coefficients[2,2]
  )
  res_df <<- bind_rows(res_df,out)
}

# plot
res_df = res_df %>%
  tibble() %>%
  arrange(desc(beta))
res_df$trait = factor(res_df$trait,levels = res_df$trait,ordered=T)
p = ggplot(res_df,aes(beta,trait))+
  geom_vline(xintercept = 0)+
  geom_errorbarh(mapping = aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se, y = trait),height=0.1)+
  geom_point(size=3)+
  theme_bw()+
  labs(x="Beta (relationship between time from earliest\n to latest reading and mean value)",y="Trait")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/span_of_readings_vs_value_betas",res=900,units="in",width=6,height=6)
p
dev.off()


# consistency within individuals

res_df = data.frame()
for(this_trait in unique(results$`Trait (full name)`)){
  this_trait_dat = raw_results %>% left_join(full_trait_names,by="trait") %>%
    filter(`Trait (full name)` == this_trait)

  # find people with >1 reading
  ppl_with_multiple_readings = this_trait_dat %>% dplyr::count(pseudo_nhs_number) %>% filter(n>=2)
  n_pre = nrow(this_trait_dat)
  this_trait_dat = this_trait_dat %>%
    filter(pseudo_nhs_number %in% ppl_with_multiple_readings$pseudo_nhs_number)
  n_post = nrow(this_trait_dat)
  message("excluded ",n_pre-n_post," records")

  # compare each person's variance with total variance
  variances = this_trait_dat %>%
#    mutate(value = RNOmni::RankNorm(value)) %>%
    group_by(pseudo_nhs_number) %>%
    mutate(var_person = var(value)) %>%
    ungroup() %>%
    mutate(overall_var = var(value)) %>%
    mutate(more_var_than_pop = ifelse(var_person > overall_var,"more_variable_than_pop","less_variable_than_pop")) %>%
    distinct(pseudo_nhs_number,more_var_than_pop,trait)
  out = variances %>% dplyr::count(more_var_than_pop) %>%
    mutate(trait = this_trait, total = sum(n), prop = n/sum(n))
  res_df <<- bind_rows(res_df,out)
}

plot_dat = res_df %>%
  tibble() %>%
  filter(more_var_than_pop == "less_variable_than_pop") %>%
  arrange(desc(prop))
plot_dat$trait = factor(plot_dat$trait,levels = plot_dat$trait,ordered=T)
plot_dat$prop %>% summary()

p=ggplot(plot_dat,aes(prop*100,trait))+
  geom_point(size=3)+
  theme_bw()+
  labs(x="% of people with >=2 readings\nwhose intra-individual variance < population variance",y="Trait")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/individual_variance",res=900,units="in",width=6,height=6)
p
dev.off()

# repeat for primary care only

res_df = data.frame()
for(this_trait in unique(results$`Trait (full name)`)){
  this_trait_dat = raw_results %>% filter(source == "primary_care") %>%
    left_join(full_trait_names,by="trait") %>%
    filter(`Trait (full name)` == this_trait)

  # find people with >1 reading
  ppl_with_multiple_readings = this_trait_dat %>% dplyr::count(pseudo_nhs_number) %>% filter(n>=2)
  n_pre = nrow(this_trait_dat)
  this_trait_dat = this_trait_dat %>%
    filter(pseudo_nhs_number %in% ppl_with_multiple_readings$pseudo_nhs_number)
  n_post = nrow(this_trait_dat)
  message("excluded ",n_pre-n_post," records")

  # compare each person's variance with total variance
  variances = this_trait_dat %>%
#    mutate(value = RNOmni::RankNorm(value)) %>%
    group_by(pseudo_nhs_number) %>%
    mutate(var_person = var(value)) %>%
    ungroup() %>%
    mutate(overall_var = var(value)) %>%
    mutate(more_var_than_pop = ifelse(var_person > overall_var,"more_variable_than_pop","less_variable_than_pop")) %>%
    distinct(pseudo_nhs_number,more_var_than_pop,trait)
  out = variances %>% dplyr::count(more_var_than_pop) %>%
    mutate(trait = this_trait, total = sum(n), prop = n/sum(n))
  res_df <<- bind_rows(res_df,out)
}

plot_dat = res_df %>%
  tibble() %>%
  filter(more_var_than_pop == "less_variable_than_pop") %>%
  arrange(desc(prop))
plot_dat$trait = factor(plot_dat$trait,levels = plot_dat$trait,ordered=T)
plot_dat$prop %>% summary()

p=ggplot(plot_dat,aes(prop*100,trait))+
  geom_point(size=3)+
  theme_bw()+
  labs(x="% of people with >=2 readings\nwhose intra-individual variance < population variance",y="Trait")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/individual_variance_primary_care",res=900,units="in",width=6,height=6)
p
dev.off()

# overall variance
variances = raw_results %>%
  group_by(trait,source) %>%
  mutate(variance = var(value)) %>%
  distinct(trait, source, variance) %>%
  ungroup() %>%
  pivot_wider(id_cols = trait,values_from = variance,names_from = source) %>%
  mutate(ratio = secondary_care / primary_care)
variances = variances %>%
  arrange(desc(ratio))
plot_dat = variances %>% left_join(full_trait_names,by="trait")
plot_dat$trait = factor(plot_dat$`Trait (full name)`,levels = plot_dat$`Trait (full name)`,ordered=T)

p=ggplot(plot_dat,aes(ratio,trait))+
  geom_point(size=3)+
  scale_x_log10()+
  theme_bw()+
  labs(x="Ratio of secondary care variance\n to primary care variance",y="Trait")+
  geom_vline(xintercept=1,linetype="dashed")

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/overall_variance_primary_care_vs_secondary",res=900,units="in",width=6,height=6)
p
dev.off()





# relationship between primary and secondary care

plots = list()
indiv_means = raw_results %>% group_by(pseudo_nhs_number,source, trait) %>%
  summarise(mean_val = mean(value))
indiv_means = indiv_means %>%
  ungroup %>%
  pivot_wider(id_cols = c(pseudo_nhs_number,trait),
              values_from = mean_val,
              names_from = source)
indiv_means = indiv_means %>%
  left_join(full_trait_names,by="trait")
indiv_means = indiv_means %>% filter(trait != "ESR")

pop_means = indiv_means %>% group_by(`Trait (full name)`) %>%
  summarise(mean_primary_care = mean(primary_care,na.rm=T), mean_secondary_care = mean(secondary_care,na.rm=T))

# plot
p = ggplot(pop_means,aes(mean_primary_care,mean_secondary_care,fill=`Trait (full name)`,label = `Trait (full name)`)) +
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  ggrepel::geom_text_repel(min.segment.length = 0,force_pull = 0,force=10)+
  geom_point(size=3,color="black",shape=21) +
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position="none")+
  labs(x="Mean in primary care data",y="Mean in secondary care data")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/primary_vs_secondary_care_popmeans.png",res=900,units="in",width=8,height=8)
p
dev.off()

for(this_trait in unique(indiv_means$`Trait (full name)`)){
  message(this_trait)
  plot_dat = indiv_means %>% filter(`Trait (full name)` == this_trait)
  colnames(plot_dat)[c(3,4)] = paste0("mean_",colnames(plot_dat)[c(3,4)])
  spearman = cor.test(plot_dat$mean_primary_care,plot_dat$mean_secondary_care,method="spearman")
  t_test = t.test(plot_dat$mean_primary_care,plot_dat$mean_secondary_care)

  lab = paste0(ifelse(sign(t_test$statistic)[[1]]==-1,
         "Higher in secondary care",
         "Lower in secondary care"))

  graph_label = paste0(this_trait,"\nRho = ",round(spearman$estimate,2),"\n",lab)
  plots[[length(plots)+1]] = ggplot(plot_dat,aes(mean_primary_care,mean_secondary_care))+
    geom_point()+
    ggtitle(graph_label)+
    theme_minimal()+
    labs(x="Primary care mean",y="Secondary care mean")+
    geom_abline(intercept=0,slope=1)

}

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/primary_vs_secondary_care.png",res=900,units="in",width=16,height=16)
do.call("grid.arrange",plots)
dev.off()

# slope
delta_values = results %>%
  filter(n_values>1) %>%
  mutate(delta = latest_value - earliest_value,
         delta_time = as.numeric(latest_date - earliest_date) / 365.25) %>%
  mutate(delta_time > 1) %>%
  dplyr::select(1,2,3,delta,delta_time) %>%
  mutate(delta_per_year = delta / delta_time)

hba1c_delta = delta_values %>%
  filter(trait == "HbA1c")

summary(hba1c_delta$delta_per_year)

# Z-score normalise
results = results %>%
  group_by(trait) %>%
  mutate(z_score = (mean - mean(mean)) / sd(mean))

# correlation
results_wide = results %>%
  pivot_wider(id_cols = pseudo_nhs_number,
              names_from = `Trait (full name)` ,
              values_from = mean)

# correlations
cor = cor(results_wide[,-1], use = "complete.obs")
cor_pval = corrplot::cor.mtest(results_wide[,-1], use = "complete.obs")
n = ncol(cor)
sig_levels = 0.05 / ( ( n - 1 ) / 2 * n )

# nice function to print correlation & P
print_corr = function(trait){
  rho = cor[rownames(cor)==trait,]
  p = cor_pval$p[rownames(cor_pval$p)==trait,]
  df = data.frame(rho,p)
  df$target = rownames(df)
  df %>%
    dplyr::arrange(desc(abs(rho))) %>%
    tibble() %>%
    dplyr::select(3,1,2)
}

print_corr("HbA1c")
print_corr("Urea")
print_corr("Serum_Na")
png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/cor_plot.png",res=600,units="in",width=10,height=10)
corrplot::corrplot(cor,
                   method = "circle",
                   order = "alphabet",
                   p.mat = cor_pval$p,
                   sig.level = sig_levels,
                   insig="blank",
                   diag = F,
                   type="lower",
                   tl.col = "black",tl.srt = 45,
                   tl.cex = 0.7
)    
dev.off()

multiple_readings = results %>%
  filter(n_values>1)
cor.test(multiple_readings$min,multiple_readings$max)
cor.test(multiple_readings$earliest_value,multiple_readings$latest_value)


traits = unique(results$trait)
# confounders
gender_res = list()
results$gender = relevel(factor(results$gender),ref="M")
for(trait_to_test in traits){
  # filter
  this_trait = results %>%
    filter(trait == trait_to_test)

  # normalise
  this_trait$norm_trait = RNOmni::RankNorm(this_trait$mean)

  # regression on age
  null_model = lm(
    norm_trait ~ gender,
    data = this_trait
  )
  summ = summary(null_model)
  r2 = summ$r.squared
  res = summary(null_model)$coefficients %>%
    data.frame() %>%
    mutate(trait = trait_to_test, r2 = r2)
  res = res[!grepl("Intercept",rownames(res)),]
  res$var = rownames(res)
  gender_res[[length(gender_res)+1]] = res
}

gender_res = do.call("bind_rows",gender_res)
colnames(gender_res)= c("beta","se","t","p","trait","r2","variable")
gender_res = gender_res %>%
  tibble %>%
  mutate(r2_pct = r2*100) %>%
  arrange(desc(r2))
gender_res %>%
  print(n=100)

summary(gender_res$r2_pct)

# all confounders

# list files
file_to_read = paste0(
  "/genesandhealth/red/quantitative_traits_v2/outputs/",
  unique(results$trait),
  "_readings_at_unique_timepoints_readings_2023-12-15.csv"
)

# read in files
results_timepoints = readr::read_csv(file_to_read, show_col_types = FALSE)

#
get_r2_from_model = function(model){
  summ = summary(model)
  r2 = summ$r.squared
  r2
}  

confounder_res = list()
for(trait_to_test in traits){
  message(trait_to_test)
  # filter
  this_trait = results_timepoints %>%
    filter(trait == trait_to_test)

  # normalise
  this_trait$norm_trait = RNOmni::RankNorm(this_trait$value)

  # get year
  this_trait$year_of_test = as.numeric(format(this_trait$date,format="%Y"))

  # regression on age
  gender_model = lm(
    norm_trait ~ gender,
    data = this_trait
  )
  age_model = lm(
    norm_trait ~ age_at_test,
    data = this_trait
  )
  source_model = lm(
    norm_trait ~ source,
    data = this_trait
  )
  year_model = lm(
    norm_trait ~ year_of_test,
    data = this_trait
  )
  combined_model = lm(
    norm_trait ~ gender + age_at_test + source + year_of_test,
    data = this_trait
  )

  all_models = list(gender_model,age_model,source_model,year_model,combined_model)
  # get r2

  r2_vals = lapply(all_models,get_r2_from_model)%>% unlist()
  res = data.frame(
    trait_to_test,
    r2_vals,
    variable = c("Gender","Age","Source","Year","Combined")
  )
  confounder_res[[length(confounder_res)+1]] = res
}

confounder_res = do.call("bind_rows",confounder_res)

summary_confounders = confounder_res %>%
  mutate(r2_pct = r2_vals * 100) %>%
  group_by(variable) %>%
  summarise(median_r2_pct = median(r2_pct),
            range_r2_pct = range(r2_pct)
            )

confounder_res_wide = confounder_res %>%
  mutate(r2_pct = r2_vals * 100) %>%
  pivot_wider(id_cols = trait_to_test,values_from = r2_pct,names_from = variable)

# plot
confounder_res$variable = factor(confounder_res$variable,levels = c("Age","Gender","Source","Year","Combined"),ordered=T)

p=ggplot(confounder_res %>% mutate("trait" = trait_to_test) %>%
           left_join(full_trait_names,by="trait"),
       aes(`Trait (full name)`,r2_vals*100,fill=variable))+
         geom_col(color="black")+
  theme_minimal()+
  scale_fill_brewer(palette = "Paired")+
  coord_flip()+
  facet_wrap(~variable,nrow=1)+
  labs(y="% variance explained",x="Trait",fill="Model")

confounder_res_wide = confounder_res_wide %>% arrange(desc(Source))
confounder_res_wide$trait_to_test = factor(confounder_res_wide$trait_to_test,levels = confounder_res_wide$trait_to_test,ordered = T)
ggplot(confounder_res_wide,
       aes(trait_to_test,Source))+
         geom_col(color="black")+
  theme_minimal()+
  scale_fill_brewer(palette = "Paired")+
  coord_flip()+
  labs(y="% variance explained",x="Trait")

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/var_explained.png",res=900,units="in",width=10,height=5)
p
dev.off()

confounder_res_wide %>% arrange(desc(Combined))

plots = list()
# distributions
for(trait_to_test in traits){
  # filter
  this_trait = results %>%
    filter(trait == trait_to_test)



  plots[[length(plots)+1]] = ggplot(this_trait,
         aes(mean)
  )+geom_histogram(color="black")+
    labs(x=trait_to_test)+
    theme_minimal()
}

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/all_histograms.png",res=900,units="in",width=12,height=12)
grid.arrange(grobs=plots)
dev.off()


plots2 = list()
# distributions
for(trait_to_test in traits){
  # filter
  this_trait = results %>%
    filter(trait == trait_to_test)

  # log10
  this_trait$log10_trait = log10(this_trait$mean)

  plots2[[length(plots2)+1]] = ggplot(this_trait,
         aes(log10_trait)
  )+geom_histogram(color="black")+
    labs(x=this_trait$`Trait (full name)`[1])+
    theme_minimal()
}

png("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/all_log10_histograms.png",res=900,units="in",width=12,height=12)
grid.arrange(grobs=plots2,ncol=4)
dev.off()

```

### Genotype QC
```unix
cd /home/ivm/genetics

## Write list of good variants to keep
rm good_vars_maf_0.01_info_0.3
# filter by INFO > 0.3 and MAF > 0.01
for i in {1..22}
do
	awk '{if($7>0.01 && $10 > 0.3) print $1,$7,$9,$10}' /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Oct2023_51k_TOPMED-r2_Imputation_b38/unfiltered_merged_vcfs/chr$i\_variant_info.tab >> good_vars_maf_0.01_info_0.3
echo "finishied copying list of good vars for chrom $i\ "
done

## Copy plink file to IVM
cd /home/ivm/genetics
cp /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Oct2023_51k_TOPMED-r2_Imputation_b38/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2_mac10_51176samples* ./

cd /home/ivm/genetics

## Basic filtering in PLINK
plink2 --pfile chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2_mac10_51176samples \
--extract good_vars_maf_0.01_info_0.3 \
--out chr_all_filtered \
--make-pgen \
--snps-only just-acgt


## More QC in PLINK
plink2 --pfile chr_all_filtered \
--snps-only just-acgt \
--maf 0.01 \
--hwe 1e-15 \
--geno 0.1 dosage \
--mind 0.1 dosage \
--make-pgen \
--out chr_all_filtered2
```

#### Update sample IDs
```R
# modify psam file
geno = read_tsv("/home/ivm/genetics/chr_all_filtered2.psam")
geno = geno %>%
mutate(old_iid = IID) %>%
mutate(old_fid = 1) %>%
separate(IID,sep="_",into=c("IID","FID","other")) %>%
select(old_fid,old_iid,FID,IID)
write_tsv(geno,"/home/ivm/genetics/new_iids")
```
#### Update sample IDs and SNP IDs as CPRA format
```unix
# rename IIDs
# rename SNPs as chr:pos:r:a
cd /home/ivm/genetics
plink2 --pfile chr_all_filtered2 \
--update-ids new_iids \
--out chr_all_filtered3 \
--make-just-psam

plink2 --pfile chr_all_filtered2 \
--set-all-var-ids @:#\$r:\$a \
--out chr_all_filtered3 \
--make-just-pvar

# rename pgen
mv chr_all_filtered2.pgen chr_all_filtered3.pgen
```
#### Prep traits for GWAS & covar file
```R
library(tidyverse)
# read link file
link_file = read_csv("/genesandhealth/library-red/genesandhealth/2023_11_08_pseudoNHSuniq_55176GSAOct2023__44028exomesJuly2023.csv")

# modify covar file
covars = read_tsv("/genesandhealth/library-red/genesandhealth/GSAv3EAMD/Oct2023_51k_TOPMED-r2_Imputation_b38/GNH.51170.noEthnicOutliers.covariates.20PCs.tab")

covars = covars %>% tidyr::separate(IID,sep="_",into = c("IID","FID","other"))
covars = covars %>%
  dplyr::select(FID,IID,4,5,10,c(11:30))
write_tsv(covars,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/covars.tsv")


# restrict to traits of interest
gwas_files = paste0("/genesandhealth/red/quantitative_traits_v2/outputs//",
  unique(results$trait),
  "_regenie_51kNov2023_GSA_pheno2023-12-15.tsv")

gwas_res = purrr::map(gwas_files,function(x){
  read_tsv(x)  %>%
    mutate(IID = str_remove_all(IID,"GNH-")) %>%
    tidyr::separate(IID,"_",into=c("IID","other")) %>%
    dplyr::select(-other)

})

# combine all files
gwas_all = do.call("bind_rows",gwas_res )

# consolidate to one person per row
gwas_all = gwas_all %>%
  group_by(FID,IID) %>%
  summarise_each(funs(mean(.,na.rm=T)))

# get rid of unneeded cols
gwas_all = gwas_all %>%
  ungroup() %>%
  dplyr::select(-pseudo_nhs_number,-trait,-unit)

psam = read_table("/home/ivm/genetics/chr_all_filtered3.psam", col_types = cols(.default ="c")) %>%
  dplyr::rename("FID" = "#FID")
pheno = gwas_all %>%
  ungroup %>%
  dplyr::select(-FID) %>%
  left_join(psam,by="IID")

pheno = pheno %>%
  dplyr::select(FID,IID,c(1:ncol(pheno)),-SEX)
write_tsv(pheno,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/pheno_modified.tsv")
```

#### Define high-quality SNPs for REGENIE step 1 (INFO > 0.98)
```R
# get high-qual snps for REGENIE step 1
pvar = read_tsv("/home/ivm/genetics/chr_all_filtered3.pvar",skip=47)

# split info field
pvar = pvar %>% separate(INFO,sep=";",into=c("AF","MAF","R2","IMPUTED","INFO","MISSING","AN","AC"))

# make numeric
pvar = pvar %>%
  mutate(INFO = as.numeric(str_remove_all(INFO,"INFO=")),
         MAF = as.numeric(str_remove_all(MAF,"MAF=")))

# filter on INFO and AF
filtered_snps = pvar %>% filter(INFO >= 0.98 & MAF > 0.01)

write_tsv(filtered_snps %>% dplyr::select(ID),
          "/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/snps_for_step1_regenie",
          col_names=F)
```

#### Basic descriptives
```R
covar = read_tsv("../outputs/covars.tsv")
pheno = read_tsv("../outputs/pheno_modified.tsv")

pheno = pheno %>% dplyr::select(FID,IID,contains("mean"))
write_tsv(pheno,"../outputs/pheno_just_means.tsv")

# summarise demogs
covar = covar %>%
    filter(IID %in% pheno$IID)
covar %>%
    group_by(pop.GSA.51k) %>%
  dplyr::count(S1QST_Gender) %>%
  mutate(pct = 100 * n/sum(n))

covar %>%
  dplyr::count(pop.GSA.51k) %>%
  mutate(pct = 100 * n/sum(n))

z_score = function(x){
  ( x - mean(x) ) / sd(x)
}
covar %>%
  mutate(PC1 = z_score(PC1),
         PC2 = z_score(PC2) ) %>%
  group_by(pop.GSA.51k) %>%
  dplyr::summarise_at(.vars = vars(AgeAtRecruitment,PC1,PC2),
                    .funs = c(mean,sd))

# summary values
summary_values = pheno %>%
  dplyr::select(contains("mean_primary_care")) %>%
  summarise_all(.funs = c("mean" = mean,"sd" = sd,"min" = min,"max" = max),
               na.rm=T) %>%
  pivot_longer(cols = everything()) %>%
  separate(name,sep = ".mean_primary_care_",into=c("variable","parameter")) %>%
  pivot_wider(id_cols = variable,values_from = value,names_from = parameter)
counts = pheno %>%
  dplyr::select(contains("mean_primary_care")) %>%
  pivot_longer(cols=everything()) %>%
  group_by(name) %>%
  dplyr::count(vals = !is.na(value)) %>%
  filter(vals == T) %>%
  mutate(variable = str_remove_all(name,"\\.mean_primary_care")) %>%
  ungroup() %>%
  dplyr::select(variable,n)
shortlist = read_csv("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/inputs/trait_shortlist.csv")

summary_values = summary_values %>%
  left_join(counts,by="variable") %>%
  dplyr::rename("trait" = variable) %>%
  left_join(shortlist,by="trait") %>%
  dplyr::select(7,1,8,2,3,4,5,6)
write_csv(summary_values,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/table2.csv")
write_csv(summary_values,"/home/ivm/table2.csv")

```
#### Run REGENIE GWAS
```unix
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 1 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--extract snps_for_step1_regenie \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--bsize 1000 \
--lowmem \
--qt \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out /home/ivm/regenie_temp/regenie_step1

### step 2
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 2 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--bsize 1000 \
--qt \
--apply-rint \
--pThresh 0.01 \
--pred /home/ivm/regenie_temp/regenie_step1_pred.list \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out quant_trait_gwas \
--lowmem

```

#### Conditional PIEZO1 GWAS for HbA1c
```unix

cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

echo "16:88716656G:T" > snps_to_condition_hba1c
regenie \
--step 1 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--extract snps_for_step1_regenie \
--phenoFile pheno_just_means.tsv \
--phenoCol HbA1c.mean_primary_care \
--condition-list snps_to_condition_hba1c \
--covarFile covars.tsv \
--bsize 1000 \
--lowmem \
--qt \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out /home/ivm/regenie_temp/regenie_step1

### step 2
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 2 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--phenoCol HbA1c.mean_primary_care \
--condition-list snps_to_condition_hba1c \
--bsize 1000 \
--qt \
--apply-rint \
--pThresh 0.01 \
--pred /home/ivm/regenie_temp/regenie_step1_pred.list \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out quant_trait_gwas_conditioned_hba1c \
--lowmem

```

#### HbA1c GWAS stratified by DM status
##### Prep phenotype
```R
library(tidyverse)

all_icd_codes = read_tsv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/3digitICD10/regenie_report_3dICD10_GSAIDs.tab")

codes = c("E08","E09","E10","E11","E13")
dm_cases = all_icd_codes %>%
  dplyr::select(FID,IID,all_of(codes)) %>%
  filter_at(.vars = c(3:7),.vars_predicate = any_vars(. == 1)) %>%
  dplyr::select(IID) %>%
  separate(IID,sep="_",into=c("IID","FID","other")) %>%
  dplyr::select(FID,IID)
write_tsv(dm_cases,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv")

# prep for binary trait gwas
counts = all_icd_codes %>%
  dplyr::select(-c(1:2)) %>%
  colSums()
high_prev_diseases = counts[counts>1000]
all_cases = all_icd_codes %>%
  dplyr::select(FID,IID,all_of(names(high_prev_diseases))) %>%
  separate(IID,sep="_",into=c("IID","FID","other")) %>%
  dplyr::select(-other,FID,IID,c(1:165))
write_tsv(all_cases,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/all_icd_pheno.tsv")

```
##### Run GWAS
```unix
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 1 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--extract snps_for_step1_regenie \
--phenoFile pheno_just_means.tsv \
--phenoCol HbA1c.mean_primary_care \
--remove /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv \
--covarFile covars.tsv \
--bsize 1000 \
--lowmem \
--qt \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out /home/ivm/regenie_temp/no_diabetics_hba1c_regenie_step1

### step 2
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 2 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--phenoCol HbA1c.mean_primary_care \
--remove /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv \
--bsize 1000 \
--qt \
--apply-rint \
--pThresh 0.01 \
--pred /home/ivm/regenie_temp/no_diabetics_hba1c_regenie_step1_pred.list \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out quant_trait_gwas_no_diabetics_hba1c \
--lowmem

cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 1 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--extract snps_for_step1_regenie \
--phenoFile pheno_just_means.tsv \
--phenoCol HbA1c.mean_primary_care \
--keep /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv \
--covarFile covars.tsv \
--bsize 1000 \
--lowmem \
--qt \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out /home/ivm/regenie_temp/just_diabetics_hba1c_regenie_step1

### step 2
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 2 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--phenoCol HbA1c.mean_primary_care \
--keep /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv \
--bsize 1000 \
--qt \
--apply-rint \
--pThresh 0.01 \
--pred /home/ivm/regenie_temp/just_diabetics_hba1c_regenie_step1_pred.list \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out quant_trait_gwas_just_diabetics_hba1c \
--lowmem


```

#### Statin-stratified LDL GWAS
##### Prep phenotype
```R
library(tidyverse)

all_meds = read_csv("../../../Medications/outputs/cleaned_prescription_data_16_11_23.csv")

# find statins
statins = all_meds %>%
  filter(grepl("statin",drugname)) %>%
  filter(!grepl("nystatin",drugname))

# get full FIDs
ldl = read_tsv("../../outputs/LDL-C_regenie_51kNov2023_GSA_pheno2023-12-15.tsv")

statins = statins %>%
  distinct(pseudo_nhs_number) %>%
  left_join(ldl,by="pseudo_nhs_number") %>%
  dplyr::select(IID) %>%
  separate(IID,sep="_",into=c("IID","FID","Other"))
statins = statins %>% dplyr::select(FID,IID) %>%
  filter(!is.na(IID))
write_tsv(statins,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/statin_pts_to_exclude.tsv")

```
##### Run GWAS
```unix
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 1 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--extract snps_for_step1_regenie \
--phenoFile pheno_just_means.tsv \
--phenoCol LDL-C.mean_primary_care \
--remove /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/statin_pts_to_exclude.tsv \
--covarFile covars.tsv \
--bsize 1000 \
--lowmem \
--qt \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out /home/ivm/regenie_temp/no_statins_ldl_regenie_step1

### step 2
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/

regenie \
--step 2 \
--pgen /home/ivm/genetics/chr_all_filtered3 \
--phenoFile pheno_just_means.tsv \
--covarFile covars.tsv \
--phenoCol LDL-C.mean_primary_care \
--remove /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/statin_pts_to_exclude.tsv \
--bsize 1000 \
--qt \
--apply-rint \
--pThresh 0.01 \
--pred /home/ivm/regenie_temp/no_statins_ldl_regenie_step1_pred.list \
--catCovarList S1QST_Gender,pop.GSA.51k \
--covarColList AgeAtRecruitment,PC{1:20} \
--apply-rint \
--out quant_trait_gwas_no_statins_ldl \
--lowmem


```

#### Analyze DM-stratified GWAS
```R
library(tidyverse)
hba1c_nondiab = read_table("../outputs/quant_trait_gwas_no_diabetics_hba1c_HbA1c.mean_primary_care.regenie")
hba1c_diab = read_table("../outputs/quant_trait_gwas_just_diabetics_hba1c_HbA1c.mean_primary_care.regenie")
hba1c_original = read_table("../outputs/quant_trait_gwas_HbA1c.mean_primary_care.regenie")

# check if any new signals
nondiab_sig = hba1c_nondiab %>%
  filter(LOG10P > -log10(1.72e-9))

# loop through each
n_sig_hits = list()
for(i in c(1:nrow(nondiab_sig))){
  this_snp = nondiab_sig[i,]

  # filter original gwas
  original_gwas_this_locus = hba1c_original %>%
    filter(CHROM == this_snp$CHROM & GENPOS > this_snp$GENPOS -5e5 & GENPOS < this_snp$GENPOS +5e5)

  n_sig_hits[[i]] = original_gwas_this_locus %>%
    filter(LOG10P > -log10(1.72e-9)) %>%
    nrow()

}

nondiab_sig = nondiab_sig %>%
  mutate(n_sig_hits_original_gwas = unlist(n_sig_hits))


# miami plot


# join
combo_gwas = hba1c_diab %>%
  mutate(cohort = "Diabetics") %>%
  bind_rows(
    hba1c_nondiab %>%
      mutate(cohort = "Non-diabetics")
  ) %>%
  mutate(P = 10^ - LOG10P) %>%
  dplyr::rename("CHR" = CHROM,"BP" = GENPOS)


# miami plot

chr_coords = combo_gwas %>%
  group_by(CHR) %>%
  summarise(min_bp = min(BP), max_bp = max(BP)) %>%
  mutate(cum_bp = cumsum(max_bp)) %>%
  mutate(CHR = CHR + 1)

combo_gwas = combo_gwas %>%
  left_join(chr_coords,by="CHR") %>%
  mutate(genpos = BP + cum_bp) %>%
  mutate(genpos = ifelse(CHR==1,BP,genpos))


  combo_gwas = combo_gwas %>%
    mutate(color_chrom = case_when(
      P < 1.72e-9 ~ "sig",
      CHR%%2==0 ~ "even",
      CHR%%2!=0 ~ "odd"
))


midpoints = combo_gwas %>%
  group_by(CHR) %>%
  summarise(midpoint = median(genpos))
pal = c("lavenderblush1","lavenderblush2","red")
pal = setNames(pal,c("even","odd","sig"))

# negate P vals for diabetics
combo_gwas = combo_gwas %>%
  mutate(log10p = ifelse(cohort == "Diabetics",-LOG10P,LOG10P))

# truncate to 30
combo_gwas = combo_gwas %>%
  mutate(log10p = ifelse(log10p>30,30,log10p)) %>%
  mutate(log10p = ifelse(log10p<(-30),-30,log10p))

p=ggplot(combo_gwas,
aes(genpos,log10p,color=color_chrom))+
geom_point()+
scale_color_manual(values = pal)+
theme_bw()+
theme(panel.grid = element_blank())+
geom_hline(yintercept=-log10(1.72e-9),linetype="dashed")+
geom_hline(yintercept=log10(1.72e-9),linetype="dashed")+
theme(legend.position = "none")+
scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHR)+
labs(x="Chromosome")+
labs(y="-log10(P) for non-diabetics (upper)\nlog10(P) for diabetics (lower)")

png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/miami_plots/miami_hba1c_diab_vs_non.png"),res=900,units="in",height=6,width=12)
print(p)
dev.off()


# find sig nondiab snps
sig_nondiab_snps = hba1c_nondiab %>%
  filter(LOG10P > -log10(5e-8)) %>%
  left_join(hba1c_original,by=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1"))
ggplot(sig_nondiab_snps,aes(BETA.x,BETA.y))+
  geom_point()+
  theme_bw()+
  labs(x="Non-diabetic beta",y="Whole cohort beta")+
  geom_abline(intercept=0,slope=1)
ggplot(sig_nondiab_snps,aes(LOG10P.x,LOG10P.y))+
  geom_point()+
  theme_bw()+
  labs(x="Non-diabetic log10p",y="Whole cohort log10p")+
  geom_abline(intercept=0,slope=1)+
  geom_vline(xintercept = -log10(5e-8))+
  geom_hline(yintercept = -log10(5e-8))

sig_nondiab_snps %>%
  filter(LOG10P.y < -log10(1e-5)) %>%
  print(n=300)


# look at piezo1
chr = 16
pos = 88716656

# locus plot
hba1c_nondiab = hba1c_nondiab %>% mutate(GWAS = "Non-diabetics") %>% filter(CHROM == chr & GENPOS > pos - 5e+5  & GENPOS < pos + 5e+5)
hba1c_diab = hba1c_diab %>% mutate(GWAS = "Diabetics") %>% filter(CHROM == chr & GENPOS > pos - 5e+5  & GENPOS < pos + 5e+5)

p = ggplot(hba1c_nondiab,aes(GENPOS,LOG10P))+geom_point()+
  theme_bw()+
  ggtitle("HbA1c GWAS: chr16:88250000 - 89250000\nin nondiabetics")+
  geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  labs(x="Position on chromosome 16",y="-log10(P)")+
  ggrepel::geom_text_repel(data = hba1c_nondiab %>% filter(GENPOS == pos),aes(label="rs563555492"),min.segment.length = 0,nudge_y=1)+
  ggrepel::geom_text_repel(data = hba1c_nondiab %>% filter(GENPOS == 88649138),aes(label="rs72549200"),min.segment.length = 0,nudge_y = 1)

p2 = ggplot(hba1c_diab,aes(GENPOS,LOG10P))+geom_point()+
  theme_bw()+
  ggtitle("HbA1c GWAS: chr16:88250000 - 89250000\nin diabetics")+
  geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  labs(x="Position on chromosome 16",y="-log10(P)")+
  ggrepel::geom_text_repel(data = hba1c_diab %>% filter(GENPOS == pos),aes(label="rs563555492"),min.segment.length = 0,nudge_y=1)+
  ggrepel::geom_text_repel(data = hba1c_diab %>% filter(GENPOS == 88649138),aes(label="rs72549200"),min.segment.length = 0,nudge_y = 1)



png("../outputs/piezo1_diab_vs_non.png",res=900,height=6,width=8,units="in")
cowplot::plot_grid(p2,p,align="v",ncol=1)
dev.off()
```

#### Analyze LDL statin GWAS
```R
library(tidyverse)
ldl_original = read_table("../outputs/quant_trait_gwas_LDL-C.mean_primary_care.regenie")
ldl_no_statins = read_table("../outputs/quant_trait_gwas_no_statins_ldl_LDL-C.mean_primary_care.regenie")

# check if any new signals
ldl_no_statins_sig = ldl_no_statins %>%
  filter(LOG10P > -log10(1.72e-9))

# loop through each
n_sig_hits = list()
for(i in c(1:nrow(ldl_no_statins_sig))){
  this_snp = ldl_no_statins[i,]

  # filter original gwas
  original_gwas_this_locus = ldl_original %>%
    filter(CHROM == this_snp$CHROM & GENPOS > this_snp$GENPOS -5e5 & GENPOS < this_snp$GENPOS +5e5)

  n_sig_hits[[i]] = original_gwas_this_locus %>%
    filter(LOG10P > -log10(1.72e-9)) %>%
    nrow()

}

ldl_no_statins_sig = ldl_no_statins_sig %>%
  mutate(n_sig_hits_original_gwas = unlist(n_sig_hits))


# miami plot


# join
combo_gwas = ldl_no_statins %>%
  mutate(cohort = "No statins") %>%
  bind_rows(
    ldl_original %>%
      mutate(cohort = "Whole cohort")
  ) %>%
  mutate(P = 10^ - LOG10P) %>%
  dplyr::rename("CHR" = CHROM,"BP" = GENPOS)


# miami plot

chr_coords = combo_gwas %>%
  group_by(CHR) %>%
  summarise(min_bp = min(BP), max_bp = max(BP)) %>%
  mutate(cum_bp = cumsum(max_bp)) %>%
  mutate(CHR = CHR + 1)

combo_gwas = combo_gwas %>%
  left_join(chr_coords,by="CHR") %>%
  mutate(genpos = BP + cum_bp) %>%
  mutate(genpos = ifelse(CHR==1,BP,genpos))


  combo_gwas = combo_gwas %>%
    mutate(color_chrom = case_when(
      P < 1.72e-9 ~ "sig",
      CHR%%2==0 ~ "even",
      CHR%%2!=0 ~ "odd"
))


midpoints = combo_gwas %>%
  group_by(CHR) %>%
  summarise(midpoint = median(genpos))
pal = c("lavenderblush1","lavenderblush2","red")
pal = setNames(pal,c("even","odd","sig"))

# negate P vals for diabetics
combo_gwas = combo_gwas %>%
  mutate(log10p = ifelse(cohort == "Whole cohort",-LOG10P,LOG10P))

# truncate to 30
combo_gwas = combo_gwas %>%
  mutate(log10p = ifelse(log10p>30,30,log10p)) %>%
  mutate(log10p = ifelse(log10p<(-30),-30,log10p))

p=ggplot(combo_gwas,
aes(genpos,log10p,color=color_chrom))+
geom_point()+
scale_color_manual(values = pal)+
theme_bw()+
theme(panel.grid = element_blank())+
geom_hline(yintercept=-log10(1.72e-9),linetype="dashed")+
geom_hline(yintercept=log10(1.72e-9),linetype="dashed")+
theme(legend.position = "none")+
scale_x_continuous(breaks = midpoints$midpoint,labels = midpoints$CHR)+
labs(x="Chromosome")+
labs(y="-log10(P) in people not taking statins (upper)\nlog10(P) for whole cohort (lower)")

  png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/miami_plots/miami_ldl_statin_vs_no.png"),res=900,units="in",height=6,width=12)
  print(p)
  dev.off()

```
#### GWAS power calcs
```R
# power
maf = 0.03
beta = -0.4

find_power = function(beta){
pvals = list()
for(i in c(1:1000)){

# simulate genos
df = data.frame(genos = rbinom(n = 6396, size = 2, prob = maf))
rare_hom = df %>% filter(genos == 2)
het = df %>% filter(genos == 1)
hom = df %>% filter(genos == 0)

# simulate phenos
rare_hom$pheno = rnorm(n = nrow(rare_hom),mean = 0 + 2*beta,sd=1)
het$pheno = rnorm(n = nrow(het),mean = 0 + beta,sd=1)
hom$pheno = rnorm(n = nrow(hom),mean = 0 ,sd=1)

# combine
dat = bind_rows(rare_hom,het,hom)

# model
summ = summary(lm(data = dat,pheno ~ genos))$coefficients
pval = data.frame(summ)$`Pr...t..`[2]
pvals[[i]] = pval
}

power = sum(unlist(pvals)<5e-8) / 1000
power
}

betas = c(0,-0.1,-0.2,-0.3,-0.4,-0.5)
powers = sapply(betas,find_power)

power_dat = data.frame(betas,powers)
ggplot(power_dat,aes(betas,powers))

```

#### Analyse HbA1c conditional GWAS
```R
library(tidyverse)
hba1c_cond = read_table("../outputs/quant_trait_gwas_conditioned_hba1c_HbA1c.mean_primary_care.regenie")
hba1c_original = read_table("../outputs/quant_trait_gwas_HbA1c.mean_primary_care.regenie")

# look at piezo1
chr = 16
pos = 88716656

hba1c_cond = hba1c_cond %>% mutate(GWAS = "Conditioned") %>% filter(CHROM == chr & GENPOS > pos - 5e+5  & GENPOS < pos + 5e+5)
hba1c_original = hba1c_original %>% mutate(GWAS = "Unconditioned") %>% filter(CHROM == chr & GENPOS > pos - 5e+5  & GENPOS < pos + 5e+5)

p = ggplot(hba1c_cond,aes(GENPOS,LOG10P))+geom_point()+
  theme_bw()+
  ggtitle("HbA1c GWAS: chr16:88250000 - 89250000\nConditioned on rs563555492")+
  geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  labs(x="Position on chromosome 16",y="-log10(P)")
p2 = ggplot(hba1c_original,aes(GENPOS,LOG10P))+geom_point()+
  theme_bw()+
  ggtitle("HbA1c GWAS: chr16:88250000 - 89250000\nUnconditional analysis")+
  geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  labs(x="Position on chromosome 16",y="-log10(P)")+
  ggrepel::geom_text_repel(data = hba1c_original %>% filter(GENPOS == pos),aes(label="rs563555492"),min.segment.length = 0,nudge_y=1)+
  ggrepel::geom_text_repel(data = hba1c_original %>% filter(GENPOS == 88649138),aes(label="rs72549200"),min.segment.length = 0,nudge_y = 1)



png("../outputs/piezo1_condition.png",res=900,height=6,width=8,units="in")
cowplot::plot_grid(p2,p,align="v",ncol=1)
dev.off()

```

### Post-GWAS analysis
```R
library(tidyverse)
library(qqman)

path = "/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/"
# functions
filter_sumstats = function(x){
  x %>%
    rename("CHR" = CHROM, "BP"=GENPOS) %>%
    mutate(CHR = as.numeric(CHR),BP = as.numeric(BP)) %>%
    filter(!is.na(CHR) & !is.na(BP)) %>%
    filter(A1FREQ >= 0.01 & A1FREQ <= 0.99)
}
read_gwas = function(x){
  read_table(x,col_types = cols_only(
  CHROM = col_double(),
  GENPOS = col_double(),
  ID = col_character(),
  ALLELE0 = col_character(),
  ALLELE1 = col_character(),
  A1FREQ = col_double(),
  N = col_double(),
  BETA = col_double(),
  SE = col_double(),
  CHISQ = col_double(),
  LOG10P = col_double())) %>%
    mutate(pheno = x, P = 10^-LOG10P) %>%
    dplyr::select(-CHISQ,-LOG10P)
}
make_manhattan = function(x){
  message("Filtering GWAS results")
  x = x %>% filter(P<0.1)
  topsnps = x %>% filter(P<1e-5)  
  message("plotting")
  manhattan(x, snp = "ID",highlight = topsnps$ID)
}
print_top_snps = function(x){
  topsnps = x %>% filter(P<1e-5)  
  print(x %>% arrange(P),n=100)  
}


# list gwas results

# run for traits
files = list.files("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/",
           pattern="\\.regenie",full.names = T)
files = files[grepl("quant_trait_gwas_",files)]

# read in
library(doParallel)
registerDoParallel(cl = detectCores()/2)
gwas = foreach(i=1:length(files))%dopar%{
  message("Reading in trait ",i)
  x = files[i]
  trait_name = str_remove_all(str_remove_all(
    str_remove_all(
      x,"/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs//quant_trait_gwas_"),
    ".regenie"),"\\.mean")

  y = x %>%
    read_gwas() %>%
    filter_sumstats() %>%
    mutate(pheno = trait_name)

  # modify trait name
  y = y %>%
    mutate(analysis = ifelse(grepl("primary_care",pheno),"primary_care","all")) %>%
    mutate(pheno = ifelse(grepl("primary_care",pheno),str_remove_all(pheno,"_primary_care"),pheno))

  y
}
# combine into mega-file
gwas_all = do.call("bind_rows",gwas)

# split into main analysis and primary care analysis
primary_care_analysis = gwas_all %>% filter(analysis=="primary_care")
main_analysis = gwas_all %>% filter(analysis!="primary_care")

# save for faster I/O
saveRDS(main_analysis,"/home/ivm/all_gh_gwas_19_12_23.rds",compress=F)
saveRDS(primary_care_analysis,"/home/ivm/all_gh_gwas_primary_care_19_12_23.rds",compress=F)

primary_care_analysis = readRDS("/home/ivm/all_gh_gwas_primary_care_19_12_23.rds")
main_analysis = readRDS("/home/ivm/all_gh_gwas_19_12_23.rds")

traits = unique(main_analysis$pheno)
for(i in c(1:length(traits))){
#foreach(i=1:length(traits))%dopar%{

  this_trait = traits[i]
  message("Doing ",this_trait)
  main = main_analysis %>% filter(pheno == this_trait)
  pc = primary_care_analysis %>% filter(pheno == this_trait)
  dat = bind_rows(main,pc)  

  # spread
  dat_wide = dat %>% pivot_wider(id_cols = c(CHR,BP,ID,ALLELE0,ALLELE1,pheno),
                                 values_from = c(A1FREQ,N,BETA,SE,P),names_from = analysis)

  # corr test
  cor = cor.test(dat_wide$BETA_all,dat_wide$BETA_primary_care)
  p=ggplot(dat_wide %>% filter(P_all<0.05),aes(BETA_all,BETA_primary_care))+
    geom_point()+
    theme_minimal()+
    labs(x = "Beta - main analysis",y="Beta - primary care only")+
    geom_abline(intercept = 0,slope=1,linetype="dashed")+
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed")+
    ggtitle(paste0(this_trait,"\nr = ",round(cor$estimate,1)))
  png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/main_analysis_vs_primary_care/comparison_plot",this_trait,".png"),res=900,units="in",width=6,height=6)
  print(p)
  dev.off()

  # make qq plot
  dat_qq = dat_wide %>% filter(!is.na(P_all) & !is.na(P_primary_care))
  p_all = -log10(sort(dat_qq$P_all,decreasing = T))
  p_pc = -log10(sort(dat_qq$P_primary_care,decreasing = T))
  qq_dat = data.frame(cuts = (seq_along(along.with = p_all))/length(p_all) ) %>%
    arrange(desc(cuts)) %>%
    mutate(expected = -log10(cuts),p_all = p_all,p_pc = p_pc) %>%
    pivot_longer(cols = c(p_all,p_pc)) %>%
    mutate(analysis = ifelse(name == "p_all","Main analysis","Primary care only"))

  p2=ggplot(qq_dat,
         aes(expected,value,col=analysis))+
    geom_line()+
    labs(x="Observed -log10(P)",y="Expected -log10(P)",color="Analysis")+
    theme_minimal()+
    scale_color_brewer(palette="Set1")+
    geom_abline(intercept = 0,slope=1,linetype="dashed")

  png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/main_analysis_vs_primary_care/qqplot_plot",this_trait,".png"),res=900,units="in",width=6,height=6)
  print(p2)
  dev.off()

  # manhattans
    png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/main_analysis_vs_primary_care/manhattan_main_",this_trait,".png"),res=900,units="in",width=10,height=6)
  print(qqman::manhattan(dat_wide %>% filter(!is.na(P_all) & P_all < 0.05) %>% dplyr::rename("P" = P_all,"SNP" = ID) %>% mutate(P = P+1e-300),main=paste0(this_trait," Main analysis")))
  dev.off()
    png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/main_analysis_vs_primary_care/manhattan_pc_",this_trait,".png"),res=900,units="in",width=10,height=6)
print(qqman::manhattan(dat_wide %>% filter(!is.na(P_primary_care) &  P_primary_care < 0.05) %>% dplyr::rename("P" = P_primary_care,"SNP" = ID)%>% mutate(P = P+1e-300),main=paste0(this_trait," Primary care only analysis")))
  dev.off()
  }

  library(tidyverse)

  primary_care_analysis = readRDS("/home/ivm/all_gh_gwas_primary_care_19_12_23.rds")

  # main_analysis = readRDS("/home/ivm/all_gh_gwas_19_12_23.rds")
  traits = unique(primary_care_analysis$pheno)

  # just use primary care-only gwas
  gwas_all = primary_care_analysis
  rm(primary_care_analysis)

  # N SNP
  unique(gwas_all$ID) %>% length

  # set p threshold
  p_bonf = 0.05 / nrow(gwas_all)

  # get sig SNP counts
  sig_snps = gwas_all %>%
    filter(P < p_bonf)
  message(nrow(sig_snps)," sig SNPS")

  sig_snps %>%
    dplyr::count(pheno) %>%
    arrange(desc(n)) %>%
    print(n=50)

  traits = unique(gwas_all$pheno)

  # LD-clump
  #for(i in c(1:length(traits))){
  # define p threshold for clumping as 5e-8/29 (as 29 traits overlap with UKB)
  p_for_clump = 5e-8/29
  foreach(i=1:length(traits))%dopar%{

    this_trait = traits[i]
    gwas_this_trait = gwas_all %>% filter(pheno==this_trait)
    sig_hits_this_trait = gwas_this_trait %>% filter(P < p_bonf)
    if(nrow(sig_hits_this_trait)>0){
    cmd = paste0("plink2 --pfile /home/ivm/genetics/chr_all_filtered3 --clump ",
    "/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/quant_trait_gwas_",this_trait,".mean_primary_care.regenie ",
    "--clump-log10 --clump-p-field LOG10P --clump-kb 1000 --clump-p1 ",p_for_clump," --clump-r2 0.001 ",
    "--out /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/clump_results/clumped_",this_trait)
    system(cmd)
    }
  }

  # read back in
  sig_snps = list()
  for(i in c(1:length(traits))){
    this_trait = traits[i]
    file = paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/clump_results/clumped_",this_trait,".clumps")
    if(file.exists(file)){
    dat = read_table(file) %>%
      mutate(pheno = this_trait)
    sig_snps[[length(sig_snps)+1]] = dat
    }
  }

  # look at number of independent loci
  snps_to_keep = do.call("bind_rows",sig_snps)
  sig_snp_counts$n %>% sum()
  nrow(snps_to_keep)
  snps_to_keep %>%
    dplyr::count(pheno) %>%
    summarise(
      median(n),
      range(n)
    )

  snps_to_keep %>%
    dplyr::count(pheno) %>%
    arrange(n)

  # discovery ~ n
  gwas_n = gwas_all %>% distinct(pheno,N)
  sig_snp_counts = snps_to_keep %>%
    group_by(pheno) %>%
    dplyr::count() %>%
    left_join(gwas_n,by="pheno")

  cor.test(sig_snp_counts$N,sig_snp_counts$n)

  # bring in full trait names
  shortlist = read_csv("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/inputs/trait_shortlist.csv")
  sig_snp_counts = sig_snp_counts %>% left_join(shortlist %>%
                                 dplyr::rename("pheno" = trait),
                               by="pheno")
  p = ggplot(sig_snp_counts,
         aes(N,n,label=`Trait (full name)`))+
    geom_point()+
    ggrepel::geom_text_repel(min.segment.length = 0,force_pull = 0)+
    theme_bw()+
    labs(x="N individuals",y="N significant loci")
  png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/sig_loci.png"),res=900,units="in",width=10,height=10)
  print(p)
  dev.off()

  # inflation

  chisq = gwas_all %>%
    dplyr::select(ID,pheno,P) %>%
    mutate(chisq = qchisq(1-P,1))

  inflation = chisq %>%
    group_by(pheno) %>%
    summarise(med_chisq = median(chisq)) %>%
    mutate(lambda = med_chisq / qchisq(0.5,1))

  inflation = inflation %>%
    arrange(desc(lambda))
  inflation = inflation %>%
    left_join(shortlist %>% dplyr::rename("pheno" = trait),
              by="pheno")
  inflation %>% summarise_at(.vars = "lambda",.funs = c("median","range"))
  inflation$`Trait (full name)` = factor(inflation$`Trait (full name)`,levels = inflation$`Trait (full name)`,ordered = T)

  p = ggplot(inflation,aes(lambda,`Trait (full name)`))+
    geom_point()+
    theme_minimal()+
    geom_vline(xintercept = 1,linetype="dashed")+
    labs(x="Genomic inflation (lambda)",y="Trait")
  png(paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/inflation"),res=900,units="in",width=10,height=10)
  print(p)
  dev.off()

```
### PIEZO1 analysis
#### Get hard-called genotypes
```unix
cd /genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/
plink2 --pfile /home/ivm/genetics/chr_all_filtered3 \
--snp 16:88716656G:T \
--recode AD \
--out piezo1_genotypes_snp1_55k
```

#### Explore relation to phenotypes
```R
library(tidyverse)
piezo1 = read_table("../outputs/piezo1_genotypes_snp1_55k.raw")


# join with pheno & hard-call
phenotype_data = read_tsv("../outputs/pheno_modified.tsv") %>%
  left_join(piezo1,by="IID") %>%
  mutate(`16:88716656_G` = round(as.numeric(`16:88716656G:T_G`),0)) %>%
  mutate("rs563555492 genotype" =
           case_when(`16:88716656_G` == 2 ~ "G/G",
                     `16:88716656_G` == 1 ~ "G/T",
                     `16:88716656_G` == 0 ~ "T/T"
                     ))

# age at dx
# codes = c("E08","E09","E10","E11","E13")
codes = c("E11")
dm_age_at_dx = purrr::map(codes,function(x){
  read_csv(paste0(
    "/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/3digitICD10/3-digit-ICD/",
    x,"/",x,"_summary_report.csv"
  ),col_types="cccdcc")
})

dm_age_at_dx = do.call("bind_rows",dm_age_at_dx)

# get earliest age at dx per person
dm_age_at_dx = dm_age_at_dx %>%
  group_by(nhs_number) %>%
  slice_min(age_at_event,with_ties = F) %>%
  ungroup()

# read nhs <> GSA ID codex
link_file = read_tsv("/genesandhealth/library-red/genesandhealth/2024_01_19_pseudoNHS_oragene_withmissing.tab") %>%
  dplyr::select(3,1)
colnames(link_file) = c("nhs_number","IID")
link_file = link_file %>% filter(nhs_number %in% dm_age_at_dx$nhs_number)
link_file = link_file %>% filter(IID %in% piezo1$IID)

# add covars
covars = read_tsv("../outputs/covars.tsv")

# age at dx
dm_age_at_dx = dm_age_at_dx %>%
  left_join(link_file,by="nhs_number") %>%
  left_join(piezo1,by="IID") %>%
  left_join(covars,by="IID")
dm_age_at_dx = dm_age_at_dx %>%
  filter(!is.na(`16:88716656G:T_G`) & !is.na(age_at_event))
dm_age_at_dx = dm_age_at_dx %>%
  mutate(`16:88716656_G` = round(as.numeric(`16:88716656G:T_G`),0)) %>%
  mutate("rs563555492 genotype" =
           case_when(`16:88716656_G` == 2 ~ "G/G",
                     `16:88716656_G` == 1 ~ "G/T",
                     `16:88716656_G` == 0 ~ "T/T"
                     ))


dm_age_at_dx = dm_age_at_dx %>%
  filter(age_at_event > 18)
p=ggplot(dm_age_at_dx,
         aes(`rs563555492 genotype`,age_at_event,fill=`rs563555492 genotype`))+
  geom_violin(alpha=0.5,position = position_dodge(width=0.7),width=0.7)+
  geom_boxplot(alpha=0.5,width=0.1,position = position_dodge(width=0.7),show.legend = F,outlier.shape = NA)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(y = "Age at DM diagnosis")

dm_age_at_dx %>%
  group_by(`rs563555492 genotype`) %>%
  summarise(median(age_at_event))

# encode t allele
summary(dm_age_at_dx$age_at_event)
dm_age_at_dx$rs563555492_T = 2 - dm_age_at_dx$`16:88716656_G`
summ_model = summary(lm(data = dm_age_at_dx,age_at_event ~ rs563555492_T + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
             PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +  S1QST_Gender + pop.GSA.51k))
res_df = data.frame(summ_model$coefficients)
res_df$lower_ci = res_df$Estimate - 1.96 * res_df$Std..Error
res_df$upper_ci = res_df$Estimate + 1.96 * res_df$Std..Error


png("../outputs/piezo_age_at_dx.png",res=900,units="in",width=6,height=4)
p
dev.off()


phenotype_data %>%
  dplyr::count(`rs563555492 genotype`) %>%
  mutate(pct = n/sum(n)*100)

# plot - all
phenotype_data$gluc_bin = Hmisc::cut2(phenotype_data$Glucose_non_fasting.mean_primary_care,g=4)
levels(phenotype_data$gluc_bin) = c("Q1 (<4.7)","Q2 (4.7 - 5.3)","Q3 (5.3 - 6.4)","Q4 (>6.4)")

p=ggplot(phenotype_data %>% filter(!is.na(gluc_bin)),
         aes(gluc_bin,HbA1c.mean_primary_care,fill=`rs563555492 genotype`))+
  geom_violin(alpha=0.5,position = position_dodge(width=0.7),width=0.7)+
  geom_boxplot(alpha=0.5,width=0.1,position = position_dodge(width=0.7),show.legend = F,outlier.shape = NA)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(x="Random glucose quartile (mM)",y = "mean HbA1c")

png("../outputs/gluc_geno_a1c_plot.png",res=900,units="in",width=6,height=4)
p
dev.off()


# repeat in nondiabetics only
dm_cases = read_tsv("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/diabetic_pts_to_exclude.tsv")

phenotype_data_all = phenotype_data %>% mutate(dm_case = ifelse(IID %in% dm_cases$IID,"DM","no_DM"))

phenotype_data_all %>%
  group_by(`16:88716656_G`) %>%
  dplyr::count(dm_case) %>%
  mutate(prop = n/sum(n))
summary(
  glm(data = phenotype_data_all,factor(dm_case) ~ `16:88716656_G`,family=binomial(link="logit"))
)
chisq.test(phenotype_data_all$`16:88716656_G`,phenotype_data_all$dm_case)

phenotype_data = phenotype_data %>% filter(!IID %in% dm_cases$IID)

phenotype_data$gluc_bin = Hmisc::cut2(phenotype_data$Glucose_non_fasting.mean_primary_care,g=4)
levels(phenotype_data$gluc_bin) = c("Q1 (<4.6)","Q2 (4.6 - 5.1)","Q3 (5.1 - 5.7)","Q4 (>5.7)")

# encode dominant T allele
phenotype_data = phenotype_data %>%
  mutate(rs563555492_T_dom = ifelse(phenotype_data$`rs563555492 genotype` == "G/G","G/G","G/T or T/T"))

p=ggplot(phenotype_data %>% filter(!is.na(gluc_bin)),
         aes(gluc_bin,HbA1c.mean_primary_care,fill=rs563555492_T_dom))+
  geom_violin(alpha=0.5,position = position_dodge(width=0.7),width=0.7)+
  geom_boxplot(alpha=0.5,width=0.1,position = position_dodge(width=0.7),show.legend = F,outlier.shape = NA)+
  theme_bw()+
  scale_fill_brewer(palette="Set1")+
  labs(x="Random glucose quartile (mM)",y = "mean HbA1c",fill="rs563555492 genotype")

png("../outputs/gluc_geno_a1c_plot_nondiab.png",res=900,units="in",width=6,height=4)
p
dev.off()


summary(lm(data = phenotype_data %>%
         filter(`rs563555492 genotype` == "G/G"),log10(HbA1c) ~ (Glucose_non_fasting)))
summary(lm(data = phenotype_data %>%
         filter(`rs563555492 genotype` == "G/T"),log10(HbA1c) ~ (Glucose_non_fasting)))


summary(lm(data = phenotype_data %>%
         filter(`rs563555492 genotype` == "G/T"),HbA1c ~ Glucose_non_fasting))

cor.test(phenotype_data$HbA1c,phenotype_data$Glucose_non_fasting,method="spearman")

medians = phenotype_data %>%
  group_by(`16:88716656_G`) %>%
  summarise(med_hba1c = median(HbA1c.mean_primary_care,na.rm = T),
            qlow_hba1c = quantile(HbA1c.mean_primary_care,0.25,na.rm = T),
            qhi_hba1c = quantile(HbA1c.mean_primary_care,0.75,na.rm = T)
            )
medians

# test distros
wt = phenotype_data[phenotype_data$`rs563555492 genotype`=="G/G",]$HbA1c.mean_primary_care
het = phenotype_data[phenotype_data$`rs563555492 genotype`=="G/T",]$HbA1c.mean_primary_care
rarehom = phenotype_data[phenotype_data$`rs563555492 genotype`=="T/T",]$HbA1c.mean_primary_care
ks.test(wt,het)
ks.test(wt,rarehom)
ks.test(het,rarehom)


```

#### Within-ancestry fine-mapping
```R

# read in specific hits
library(tidyverse)

finemap_fx = function(
    phenotype_to_test = "HbA1c",
    snp="16:88649138T:C"){

    # fine map
  filename = paste0(
    "/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/quant_trait_gwas_",
    phenotype_to_test,
    ".mean_primary_care.regenie")
  sumstats_for_finemap = read_table(filename,col_types = cols(.default = "c")) %>%
    mutate(CHROM = as.numeric(CHROM),
           GENPOS = as.numeric(GENPOS),
           P = 10^-as.numeric(LOG10P))




  # keep SNPs within a 1MB window of lead snp
  this_snp = sumstats_for_finemap %>% filter(ID == snp)
  sumstats_for_finemap = sumstats_for_finemap %>%
      filter(
               CHROM == this_snp$CHROM &
               GENPOS > this_snp$GENPOS - 5e5 &
               GENPOS < this_snp$GENPOS + 5e5
             )

  n_snp = nrow(sumstats_for_finemap)
  message("N SNPs: ",n_snp)
  dup_ids = sumstats_for_finemap %>% dplyr::count(ID) %>% filter(n>=2)
  sumstats_for_finemap = sumstats_for_finemap %>% filter(!ID %in% dup_ids$ID)
  n_snp = nrow(sumstats_for_finemap)
  message("N SNPs after removing duplicates: ",n_snp)

  message(phenotype_to_test)
  # finemap  
  message("Doing locus ",this_snp$ID)
  sumstats_for_finemap_z = sumstats_for_finemap %>%
    dplyr::select(ID,CHROM,GENPOS,ALLELE0,ALLELE1,A1FREQ,BETA,SE)
  colnames(sumstats_for_finemap_z) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

  # convert A1FREQ to MAF
  sumstats_for_finemap_z = sumstats_for_finemap_z %>%
    mutate(maf = as.numeric(maf)) %>%
    mutate(maf = ifelse(maf > 0.5,1-maf,maf))

  z_filename = paste0(
                "/home/ivm/sumstats_z",
                phenotype_to_test,"_",snp,"_file.z")
  write.table(sumstats_for_finemap_z,
              z_filename,
              quote = F,row.names = F)

  # compute LD matrix for finemap


  # make bed
  cmd = paste0("plink2 --pfile /home/ivm/genetics/chr_all_filtered3 --extract ",z_filename," --make-bed --out /home/ivm/temp_plink1_ld_snps; wait")
  system(cmd)

  # calculate ld matrix
  cmd2=paste0("plink --bfile /home/ivm/temp_plink1_ld_snps --r square --out /home/ivm/temp_ld_mat; wait")
  system(cmd2)
  system("rm /home/ivm/temp_plink1_ld_snps.bed;rm /home/ivm/temp_plink1_ld_snps.bim;rm /home/ivm/temp_plink1_ld_snps.fam")

  # format LD file
  ld_file = read_tsv("/home/ivm/temp_ld_mat.ld",col_names = F)
  ld_filename = paste0(
            "/home/ivm/sumstats_",
            phenotype_to_test,"_",snp,"_ld_file.ld")

  write.table(ld_file,ld_filename,quote = F,row.names = F,col.names = F)

  n_gwas = mean(as.numeric(sumstats_for_finemap$N),na.rm=T)
  system("echo 'z;ld;snp;config;cred;log;n_samples' > /home/ivm/master_finemap_file.txt")
  system(paste0("echo '",z_filename,";",ld_filename,";/home/ivm/dataset1.snp;/home/ivm/dataset1.config;/home/ivm/dataset1.cred;/home/ivm/dataset1.log;",n_gwas,"' >> /home/ivm/master_finemap_file.txt"))

  # run finemap
  cmd = paste0("finemap --sss --n-causal-snps 1 --in-files /home/ivm/master_finemap_file.txt")
  system(cmd)


  # read in results
  post_probs = read.table("/home/ivm/dataset1.snp", header=T) %>%
    tibble() %>%
    arrange(desc(prob)) %>%
    mutate(cumprob = cumsum(prob)) %>%
    mutate(in_cred_set_99 = ifelse(cumprob <= 0.99,"yes","no")) %>%
    mutate(in_cred_set_99 = ifelse(prob ==1,"yes",in_cred_set_99))

  # plot
  post_probs_for_merge = post_probs %>% dplyr::select(rsid,prob,in_cred_set_99)
  sumstats_for_finemap_z = sumstats_for_finemap_z %>%
    left_join(post_probs_for_merge,
        by="rsid")

  n_cred = sumstats_for_finemap_z %>% filter(in_cred_set_99=="yes") %>% nrow()

  # get extra SNP info    
  extra_info = sumstats_for_finemap %>%
    filter(ID %in% sumstats_for_finemap_z$rsid) %>%
    dplyr::select(ID,P,ALLELE0,ALLELE1,A1FREQ) %>%
    dplyr::rename("rsid" = ID)

  sumstats_for_finemap_z = sumstats_for_finemap_z %>%
    left_join(extra_info,by="rsid")

  sumstats_for_finemap_z = sumstats_for_finemap_z %>% arrange(desc(prob))

  top_snp = sumstats_for_finemap_z$rsid[1]


  p=ggplot(sumstats_for_finemap_z,
         aes(position,-log10(P),col=prob,label=rsid))+
    geom_point()+
    theme_minimal()+
   scale_color_gradient2(low="grey",high="red",mid="blue",midpoint=0.1)+
    labs(y="-log10(P)",col="PIP",x="Genomic position (hg38)")+
    geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  ggrepel::geom_text_repel(data = sumstats_for_finemap_z %>% filter(in_cred_set_99 == "yes"),
                             force_pull = 0,
                             direction="x",
                             angle=0,
                             segment.size=0.2,
                             hjust=0.9,
                             min.segment.length = 0,
                             show.legend = F)+
    ggtitle(paste0(
      phenotype_to_test,
      "\nLead SNP ",top_snp,"\nSNPs in 99% cred set: ",n_cred))


  plot_filename = paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/finemap_plots/",
                         phenotype_to_test,"_",
                         snp,"_finemap_plots.png")
png(plot_filename,res=900,units="in",height=6,width=6)
print(p)
dev.off()

# get pairwise LD with top snp
ld_file = data.frame(ld_file)
colnames(ld_file) = sumstats_for_finemap$ID


ld_file = ld_file %>%
  dplyr::select(all_of(top_snp))
ld_file = ld_file %>%
  mutate(ID = sumstats_for_finemap$ID)
colnames(ld_file)[1] = "r"
colnames(ld_file)[2] = "rsid"
ld_file = ld_file %>%
  mutate(r2 = r^2)

# join
sumstats_for_finemap_z = sumstats_for_finemap_z %>%
  left_join(ld_file,by="rsid")

 p1=ggplot(sumstats_for_finemap_z %>% filter(P < 5e-8),
         aes(position,-log10(P),fill=r2,label=rsid))+
    geom_point(size=3,shape=21,color="black")+
   geom_point(data = sumstats_for_finemap_z %>% filter(P > 5e-8),alpha=0.1)+
    theme_bw()+
   scale_fill_gradient2(low="grey",high="orange",mid="purple",midpoint=0.5)+
    labs(y="-log10(P)",col="R2 with top SNP",x="Genomic position (hg38)")+
    geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  ggrepel::geom_text_repel(data = sumstats_for_finemap_z %>% filter(in_cred_set_99 == "yes"),
                             force_pull = 0,
                             force = 5,
                             min.segment.length = 0,
                             show.legend = F)+
    ggtitle(paste0(
      phenotype_to_test,
      "\nLead SNP ",top_snp,"\nSNPs in 99% cred set: ",n_cred))


  plot_filename = paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/finemap_plots/",
                         phenotype_to_test,"_",
                         snp,"_finemap_plots_r2.png")
png(plot_filename,res=900,units="in",height=6,width=6)
print(p1)
dev.off()

# write results
res_filename = paste0("/genesandhealth/red/quantitative_traits_v2/quant_traits_paper/outputs/finemap_plots/",
                   phenotype_to_test,"_",
                   snp,"_finemap_res.csv")          
write_csv(sumstats_for_finemap_z,res_filename)

}


finemap_fx(phenotype_to_test = "HbA1c",snp="16:88649138T:C")
finemap_fx(phenotype_to_test = "Bilirubin",snp="16:88649138T:C")
finemap_fx(phenotype_to_test = "MCV",snp="11:4960708A:C")
finemap_fx(phenotype_to_test = "MCV",snp="11:5026936T:A")
finemap_fx(phenotype_to_test = "RBC",snp="11:5122523A:T")
finemap_fx(phenotype_to_test = "MCV",snp="11:5207966G:C")
finemap_fx(phenotype_to_test = "Hb",snp="11:5331278G:A")

```
## Process GWAS data
### Download GH GWAS data
Code used to perform GWAS and run analysis in the Genes & Health TRE provided separately.
```unix
cd /data/scratch/hmy117/gwas_raw_results/raw/
~/google-cloud-sdk/bin/gsutil cp gs://fg-qmul-production-sandbox-2_green/forBen-2023-12-20/* ./

# unzip
tar -xvf /data/scratch/hmy117/gwas_raw_results/raw/exports_20_12_23.tar.gz
```

### Download pan-UKB GWAS
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
module load samtools
for file in *.bgz;
  do
    bgzip -d $file &
  done
```

### Liftover UKB GWAS
```unix

# liftover ukb gwas
awk 'NR>1{print "chr"$1,$2-1,$2,$1":"$2}' /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/biomarkers-30600-both_sexes-irnt.tsv > /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg19

# run liftover
/data/Wolfson-UKBB-Dobson/liftover/liftOver /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg19 \
/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
/data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/ukb_gwas_hg38.bed /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/unlifted.bed

```

## Meta-analysis
### Format UKB & GH files
```unix
qsub ~/gh_quant_traits/scripts/prepare_gh_files.sh

qsub ~/gh_quant_traits/scripts/prepare_ukb_files.sh
```
### MR-MEGA
```unix
qsub ~/gh_quant_traits/scripts/mr_mega.sh
```

### Plot ancestral PCs
```R
library(tidyverse)
setwd("/data/scratch/hmy117/gwas_raw_results/")

# read in log to check PCs
mr_mega_log = read_table("./mr_mega/Albumin_mr_mega_res.log",skip=47,n_max = 7)
mr_mega_log$PCs = str_remove_all(mr_mega_log$PCs,"../mr_mega_inputs/")
mr_mega_log$PCs = str_remove_all(mr_mega_log$PCs,".tsv")

a = ggplot(mr_mega_log,aes(PC0,PC1,label=PCs))+geom_point()+ggrepel::geom_text_repel(show.legend = F)+
  theme_minimal()+
    scale_color_brewer(palette="Paired")
b = ggplot(mr_mega_log,aes(PC1,PC2,label=PCs))+geom_point()+ggrepel::geom_text_repel(show.legend = F)+
  theme_minimal()+
  scale_color_brewer(palette="Paired")
png("/data/home/hmy117/gh_quant_traits/outputs/pc_plots.png",res=900,units="in",width=8,height=4)
gridExtra::grid.arrange(a,b,ncol=2)
dev.off()
```

### Prepare 1kg files for clumping
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
awk '{print "chr"$1,$4-1,$4,$2}' g1000_eur.bim > g1000_eur_hg19_bed &
awk '{print "chr"$1,$4-1,$4,$2}' g1000_sas.bim > g1000_sas_hg19_bed &
awk '{print "chr"$1,$4-1,$4,$2}' g1000_afr.bim > g1000_afr_hg19_bed &
awk '{print "chr"$1,$4-1,$4,$2}' g1000_eas.bim > g1000_eas_hg19_bed &
awk '{print "chr"$1,$4-1,$4,$2}' g1000_amr.bim > g1000_amr_hg19_bed &

# liftOver
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_eur_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_eur_hg38_bed unlifted.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_sas_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_sas_hg38_bed unlifted2.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_eas_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_eas_hg38_bed unlifted3.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_afr_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_afr_hg38_bed unlifted4.bed &
/data/Wolfson-UKBB-Dobson/liftover/liftOver g1000_amr_hg19_bed /data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz g1000_amr_hg38_bed unlifted5.bed &

# update SNP IDs
awk '{print $4,$3}' g1000_eur_hg38_bed > eur_vars &
awk '{print $4,$3}' g1000_sas_hg38_bed > sas_vars &
awk '{print $4,$3}' g1000_eas_hg38_bed > eas_vars &
awk '{print $4,$3}' g1000_afr_hg38_bed > afr_vars &
awk '{print $4,$3}' g1000_amr_hg38_bed > amr_vars &

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
    --make-bed &
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

### Process MR-MEGA results
```unix
qsub ~/gh_quant_traits/scripts/process_mr_mega_res.sh
```

### Clump results
```unix
qsub ~/gh_quant_traits/scripts/clump_results.sh
qsub ~/gh_quant_traits/scripts/clump_gh_results.sh
```

### Plots
```unix
qsub ~/gh_quant_traits/scripts/make_manhattans.sh
```

### Combine results for summary tables
```unix
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annotations",files)]

dat = purrr::map(files,function(x){
  read_csv(x,col_types = cols_only(
  CHROMOSOME = col_double(),
  POSITION = col_double(),
  NEA = col_character(),
  EA = col_character(),
  EAF = col_character(),
  N = col_double(),
  BETA = col_double(),
  SE = col_double(),
  P = col_double(),
  trait = col_character(),
  cohort = col_character(),
  MARKERNAME = col_character(),
  sig_snp = col_character(),
  gh_specific_locus = col_character()
)
)
})
dat = do.call("bind_rows",dat)

vep_input = dat %>%
mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
dplyr::select(1,start,end,alleles,strand)
write_tsv(vep_input,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv",col_names=F)

# repeat for meta
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = purrr::map(files,function(x){
  read_csv(x,col_types = cols_only(
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
  P.value_residual_het = col_double()
  )

)
})
dat = do.call("bind_rows",dat)


vep_input = dat %>%
mutate(start = Position,end = Position,alleles = paste0(EA,"/",NEA),strand="+") %>%
dplyr::select(Chromosome,start,end,alleles,strand)
write_tsv(vep_input,"/data/home/hmy117/gh_quant_traits/outputs/all_meta_sig_hits_for_vep.tsv",col_names=F)

```

### VEP
```unix
module load ensembl-vep
cd /data/home/hmy117/gh_quant_traits/outputs/

~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--check_existing \
-o gh_sig_hits_annotated.tsv \
--tab \
--no_check_alleles \
--nearest symbol \
--coding_only \
--sift b \
--polyphen b


~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_for_vep.tsv \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--nearest symbol \
-o nearest_gene.tsv \
--tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Feature,Feature_type,Consequence"

 ~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/all_meta_sig_hits_for_vep.tsv \
 --cache \
 --dir_cache /data/scratch/hmy117/.vep \
 --force_overwrite \
 --nearest symbol \
 -o meta_nearest_gene.tsv \
 --tab --fields "Uploaded_variation,Location,Allele,Gene,NEAREST,Feature,Feature_type,Consequence"

 # get ukb freqs
 awk '{print $1,$2,$3,$4,$16,$17,$18,$19,$20}' /data/scratch/hmy117/gwas_raw_results/pan_ukb_gwas/continuous-30020-both_sexes-irnt.tsv > /data/scratch/hmy117/ukb_freqs.tsv


```

### Process rsids in R (to get freqs)
```R
library(tidyverse)
# read in annotations
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_annotated.tsv",skip=47,col_types = cols(.default="c"))

# get  just rsids
vep_res_rsids = vep_res %>%
dplyr::select(1,Existing_variation) %>%
separate(Existing_variation,sep=",",into=c("rsid","other")) %>%
filter(grepl("^rs",rsid)) %>%
distinct(`#Uploaded_variation`,rsid)
write_tsv(vep_res_rsids %>% dplyr::select(rsid),"/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_rsids.tsv")
```

### VEP
```unix
module load ensembl-vep
cd /data/home/hmy117/gh_quant_traits/outputs/

~/ensembl-vep/vep -i /data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_rsids.tsv \
--cache \
--dir_cache /data/scratch/hmy117/.vep \
--force_overwrite \
--check_existing \
-o gh_sig_hits_annotated_rsids.tsv \
--tab \
--no_check_alleles \
--nearest symbol \
--coding_only \
--sift b \
--polyphen b \
--af_1kg \
--af_gnomadg \
--af_gnomade

````

### Combine results for summary tables (with VEP results)
```r
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annot",files)]
files = files[!grepl("rsids",files)]

# read in annotations
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_annotated.tsv",skip=47,col_types = cols(.default="c")) %>%
dplyr::select(1,Existing_variation)
vep_res_rsids = read_table("/data/home/hmy117/gh_quant_traits/outputs/gh_sig_hits_annotated_rsids.tsv",skip=73,col_types = cols(.default="c"))



# read in gh gwas
dat = read_csv(files,col_types = "ddccdddddccccc") %>%
  mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(CHROMOSOME,"_",POSITION,"_",EA,"/",NEA)) %>%
  left_join(vep_res %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  dplyr::rename("rsid" = Existing_variation) %>%
  left_join(vep_res_rsids %>%
    dplyr::rename("rsid" = `#Uploaded_variation`),
    by="rsid")

# save
write_csv(dat,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_mapped_consequences.csv")


# read VEP results
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/nearest_gene.tsv",skip=32,col_types = cols(.default="c"))

# read in gh gwas
dat = read_csv(files,col_types = "ddccdddddccccc") %>%
  mutate(start = POSITION,end = POSITION,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(CHROMOSOME,"_",POSITION,"_",EA,"/",NEA)) %>%
  left_join(vep_res %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  dplyr::select(1,2,full_snp_id,NEA,EA,EAF,N,BETA,SE,P,trait,sig_snp,gh_specific_locus,NEAREST) %>%
  distinct(full_snp_id,NEAREST,.keep_all=T)

# view gh-specific
write_csv(dat,"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_mapped.csv")
write_csv(dat %>% filter(sig_snp == "lead_sig_snp"),"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_just_lead_snps_mapped.csv")
write_csv(dat %>% filter(!is.na(gh_specific_locus)),"/data/home/hmy117/gh_quant_traits/outputs/all_gh_sig_hits_just_specific_snps_mapped.csv")

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
dat = read_csv(files,col_types = "cddccddddddcccc")

# read VEP results
vep_res = read_table("/data/home/hmy117/gh_quant_traits/outputs/meta_nearest_gene.tsv",skip=32,col_types = cols(.default="c"))

# take one gene per SNP
vep_res = vep_res %>% distinct(`#Uploaded_variation`,NEAREST)

dat = dat %>%
  mutate(start = Position,end = Position,alleles = paste0(EA,"/",NEA),strand="+") %>%
  mutate(full_snp_id = paste0(Chromosome,"_",Position,"_",EA,"/",NEA)) %>%
  left_join(vep_res %>% dplyr::rename("full_snp_id" = `#Uploaded_variation`),
  by="full_snp_id") %>%
  mutate(MarkerName = full_snp_id)

# filter to sig het SNPs

write_csv(dat %>% filter(P.value_ancestry_het < 1.72e-9),"/data/home/hmy117/gh_quant_traits/outputs/all_meta_anno.csv")

```
## Fine-mapping
### Run SuSiEx
```unix
cd /data/home/hmy117/gh_quant_traits/
qsub ./scripts/susie_prep.sh

qsub /data/home/hmy117/gh_quant_traits/scripts/susie.sh

```

### Read in fine mapping results
```R
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

# define specificity
finemap_res = finemap_res %>%
mutate(sas_specific = ifelse(`POST-HOC_PROB_POP1` > 0.8 & `POST-HOC_PROB_POP5` > 0.8 & `POST-HOC_PROB_POP2` < 0.8 & `POST-HOC_PROB_POP3` <0.8 & `POST-HOC_PROB_POP4`<0.8,"yes","no")) %>%
mutate(eur_specific = ifelse(`POST-HOC_PROB_POP1` < 0.8 & `POST-HOC_PROB_POP5` < 0.8 & `POST-HOC_PROB_POP2` > 0.8 & `POST-HOC_PROB_POP3` <0.8 & `POST-HOC_PROB_POP4` <0.8,"yes","no"))  %>%
mutate(afr_specific = ifelse(`POST-HOC_PROB_POP1` < 0.8 & `POST-HOC_PROB_POP5` < 0.8 & `POST-HOC_PROB_POP2` < 0.8 & `POST-HOC_PROB_POP3` >0.8 & `POST-HOC_PROB_POP4` <0.8,"yes","no")) %>%
mutate(eas_specific = ifelse(
  `POST-HOC_PROB_POP1` < 0.8 & `POST-HOC_PROB_POP5` < 0.8 & `POST-HOC_PROB_POP2` < 0.8 & `POST-HOC_PROB_POP3` <0.8 & `POST-HOC_PROB_POP4` >0.8,"yes","no"))
# save to file
write_csv(finemap_res,"~/gh_quant_traits/outputs/ancestry_specific_finemap.csv")

```

## Misc
### Locus plots
```unix
Rscript /data/home/hmy117/gh_quant_traits/scripts/make_locus_plots.R HbA1c
view /data/home/hmy117/gh_quant_traits/scripts/make_locus_plots.sh
```

### Pleiotropy plot GH
```r
library(tidyverse)
files = list.files("/data/home/hmy117/gh_quant_traits/outputs/",pattern="^gh_sig_hits",full.names=T)
files = files[!grepl("annotations",files)]
dat = read_csv(files,col_types = "ddccdddddccccc")

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



  p=ggplot(plot_dat,aes(BETA,trait,size=-log10(P),fill=trait))+
    geom_vline(xintercept=0,linetype="dashed",alpha=0.5)+
    geom_errorbarh(mapping= aes(y=trait,xmin = BETA - 1.96*SE,xmax = BETA+1.96*SE),color="black",height=0.1,linewidth=0.1)+
  geom_point(color="black",shape=21,show.legend=F)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  labs(size = "-log10(P)",
       x = paste0("Effect size (beta)\n of ",snp,"-",ea1," on trait"),y="Trait")
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

### Phewas
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

### Popcorn
#### Setup and calculate scores
````unix
# setup - run once
module load python
virtualenv ~/popcorn
source ~/popcorn/bin/activate

# pull repo
cd ~
git clone https://github.com/brielin/Popcorn
cd ~/Popcorn
pip install .

# after setup run once, run this
qsub ~/gh_quant_traits/scripts/popcorn_prep.sh

````
#### Analyze
````R
library(tidyverse)

input_traits = read_csv("~/gh_quant_traits/inputs/gh_ukb_gwas_paths.csv")

all_res = list()
for(i in c(1:29)){
  trait = input_traits$gh_phenotype[i]
  all_res[[i]] = read_table(paste0("/data/scratch/hmy117/popcorn.sh.o3303433.",i),skip=236,n_max=3,col_types = "cdddd") %>%
  mutate(phenotype = trait) %>%
  filter(Val == "pgi") %>%
  dplyr::rename(rg = `(obs)`)
}

all_res = do.call("bind_rows",all_res)

# order and plot
all_res = all_res %>%
  arrange(desc(rg))
all_res$phenotype = factor(all_res$phenotype,levels = all_res$phenotype,ordered=T)

png("~/gh_quant_traits/outputs/popcorn_res.png",res=900,units="in",width=4,height=5)
ggplot(all_res,aes(rg,phenotype))+
geom_point(size=3,color="black")+
geom_errorbarh(mapping = aes(y=phenotype,xmin = rg-SE, xmax = rg+SE),height=0.1)+
theme_bw()+
labs(x="Cross-ancestry genetic correlation",y="Trait")+
geom_vline(xintercept = 1,linetype="dashed")+
geom_vline(xintercept = 0,linetype="dashed")
dev.off()

write_csv(all_res,"~/gh_quant_traits/outputs/popcorn_res.csv")

````
