# This script performs differential abundance analysis at genus level

library(phyloseq)
library(ANCOMBC)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")

tax = taxonomy %>% select(-size) %>% as.data.frame()
rownames(tax) = tax$otu
tax = tax %>% 
  select(-otu) %>% 
  as.matrix(.)
tax = tax_table(tax)


## non-lesional v.s. lesional areas at baseline, controlling for individual effect

map_temp = map_sample %>% filter(visit == "V3")
keep = names(which(tapply(map_temp$lesional, map_temp$subjid, function(x){length(unique(x))}) == 2))
count_data = shared %>% 
  mutate(sample = str_replace_all(string=sample, pattern="non_lesional", replacement="non-lesional")) %>%
  separate(col="sample", into = c("unique", "subjid", "lesional","visit"), sep="_", remove=FALSE) %>% 
  filter(subjid %in% keep & visit == "V3") %>% 
  arrange(sample) %>% 
  select(-sample, -unique, -subjid, -lesional, -visit) %>% 
  as.data.frame() %>% 
  t()
col_data = map_sample %>% filter(subjid %in% keep & visit == "V3") %>% arrange(sample) 

count_data = otu_table(count_data, taxa_are_rows = TRUE)
col_data = sample_data(col_data)
data_phyloseq = phyloseq(count_data, col_data, tax)

set.seed(123)
res = ancombc(
  data_phyloseq, formula = "lesional + subjid", tax_level = "genus",
  group = NULL, p_adj_method = "fdr", max_iter = 100
  )
results = res$res
p_adjust = results$q_val %>% 
  select(taxon, `lesionallesional`) %>% 
  rename(p.adj = `lesionallesional`) 
lfc = results$lfc %>% 
  select(taxon, `lesionallesional`) %>% 
  rename(lfc = `lesionallesional`)
de_genus_lesional = inner_join(p_adjust, lfc)
de_genus_lesional = taxonomy %>% 
  select(otu, genus) %>% 
  inner_join(., de_genus_lesional, by = c("genus" = "taxon"))
write.csv(
  de_genus_lesional, 
  file = "ANCOMBC results/de_genus_lesional.csv"
  )


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in treatment responded groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional, nrs_w12)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "Y") %>% ### 
  arrange(sample) %>% 
  select(-subjid, -lesional, -trt01p, -nrs_w12) %>% 
  as.data.frame() 
rownames(count_data) = count_data$sample
count_data = count_data %>% select(-sample) %>% t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "Y") %>% 
  arrange(sample) ###

count_data = otu_table(count_data, taxa_are_rows = TRUE)
col_data = sample_data(col_data)
data_phyloseq = phyloseq(count_data, col_data, tax)

set.seed(123)
res = ancombc(
  data_phyloseq, formula = "visit + subjid", tax_level = "genus",
  group = NULL, p_adj_method = "fdr", max_iter = 100
  )
results = res$res
p_adjust = results$q_val %>% 
  select(taxon, `visitV8`) %>% 
  rename(p.adj = `visitV8`)
lfc = results$lfc %>% 
  select(taxon, `visitV8`) %>% 
  rename(lfc = `visitV8`)
de_genus_visit_treatment = inner_join(p_adjust, lfc)
de_genus_visit_treatment = taxonomy %>% 
  select(otu, genus) %>% 
  inner_join(., de_genus_visit_treatment, by = c("genus" = "taxon"))
write.csv(
  de_genus_visit_treatment, 
  file = "ANCOMBC results/de_genus_visit_treatment.responded.csv"
  )


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in treatment not responded groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional, nrs_w12)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "N") %>% ### 
  arrange(sample) %>% 
  select(-subjid, -lesional, -trt01p, -nrs_w12) %>% 
  as.data.frame() 
rownames(count_data) = count_data$sample
count_data = count_data %>% select(-sample) %>% t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "N") %>% 
  arrange(sample) ###

count_data = otu_table(count_data, taxa_are_rows = TRUE)
col_data = sample_data(col_data)
data_phyloseq = phyloseq(count_data, col_data, tax)

set.seed(123)
res = ancombc(
  data_phyloseq, formula = "visit + subjid", tax_level = "genus",
  group = NULL, p_adj_method = "fdr", max_iter = 100
  )
results = res$res
p_adjust = results$q_val %>% 
  select(taxon, `visitV8`) %>% 
  rename(p.adj = `visitV8`)
lfc = results$lfc %>% 
  select(taxon, `visitV8`) %>% 
  rename(lfc = `visitV8`)
de_genus_visit_treatment = inner_join(p_adjust, lfc)
de_genus_visit_treatment = taxonomy %>% 
  select(otu, genus) %>% 
  inner_join(., de_genus_visit_treatment, by = c("genus" = "taxon"))
write.csv(
  de_genus_visit_treatment, 
  file = "ANCOMBC results/de_genus_visit_treatment.not.responded.csv"
  )


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in treatment groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% ### 
  arrange(sample) %>% 
  select(-subjid, -lesional, -trt01p) %>% 
  as.data.frame() 
rownames(count_data) = count_data$sample
count_data = count_data %>% select(-sample) %>% t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  arrange(sample) ###

count_data = otu_table(count_data, taxa_are_rows = TRUE)
col_data = sample_data(col_data)
data_phyloseq = phyloseq(count_data, col_data, tax)

set.seed(123)
res = ancombc(
  data_phyloseq, formula = "visit + subjid", tax_level = "genus",
  group = NULL, p_adj_method = "fdr", max_iter = 100
  )
results = res$res
p_adjust = results$q_val %>% 
  select(taxon, `visitV8`) %>% 
  rename(p.adj = `visitV8`)
lfc = results$lfc %>% 
  select(taxon, `visitV8`) %>% 
  rename(lfc = `visitV8`)
de_genus_visit_treatment = inner_join(p_adjust, lfc)
de_genus_visit_treatment = taxonomy %>% 
  select(otu, genus) %>% 
  inner_join(., de_genus_visit_treatment, by = c("genus" = "taxon"))
write.csv(
  de_genus_visit_treatment, 
  file = "ANCOMBC results/de_genus_visit_treatment.csv"
  )


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in placebo groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% ### 
  arrange(sample) %>% 
  select(-subjid, -lesional, -trt01p) %>% 
  as.data.frame() 
rownames(count_data) = count_data$sample
count_data = count_data %>% select(-sample) %>% t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% 
  arrange(sample) ###

count_data = otu_table(count_data, taxa_are_rows = TRUE)
col_data = sample_data(col_data)
data_phyloseq = phyloseq(count_data, col_data, tax)

set.seed(123)
res = ancombc(
  data_phyloseq, formula = "visit + subjid", tax_level = "genus",
  group = NULL, p_adj_method = "fdr", max_iter = 100
  )
results = res$res
p_adjust = results$q_val %>% 
  select(taxon, `visitV8`) %>% 
  rename(p.adj = `visitV8`)
lfc = results$lfc %>% 
  select(taxon, `visitV8`) %>% 
  rename(lfc = `visitV8`)
de_genus_visit_placebo = inner_join(p_adjust, lfc)
de_genus_visit_placebo = taxonomy %>% 
  select(otu, genus) %>% 
  inner_join(., de_genus_visit_placebo, by = c("genus" = "taxon"))
write.csv(
  de_genus_visit_placebo, 
  file = "ANCOMBC results/de_genus_visit_placebo.csv"
  )



