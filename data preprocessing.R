# This script preprocess mothur-formatted data, including phylotype filtering and an overview for metadata
# Note: All "abundant genus" appear in R objects refer to common genus in the manuscript

library(readxl)
library(purrr)
library(broom)
library(tidyverse)

wd = "/Users/zrayw/Desktop/PN_microbiome/analysis/dat"
setwd(wd)

## Loading phylotype data
map_dir = "PN2_subject_annotation.xlsx"
tax_dir = "galderma.tx.1.cons.taxonomy"
shared_dir = "galderma.tx.shared"

## Import data as data.frame
map = read_excel(path = map_dir) %>% 
  rename_all(tolower) %>% 
  mutate(subjid = as.character(subjid)) %>% 
  filter(! subjid %in% c(5433004, 5471012, 5891001))
taxonomy = read.table(file = tax_dir, header = T) %>% 
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy = str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")
shared = read.table(file = shared_dir, header=T) %>% 
  dplyr::rename(sample = Group) %>% 
  select(-label, -numPhylos)
map_sample = shared %>%
  select(sample) %>%
  mutate(sample = str_replace_all(string=sample, pattern="non_lesional", replacement="non-lesional")) %>%
  separate(col="sample", into = c("unique", "subjid", "lesional","visit"), sep="_", remove=FALSE) %>% 
  select(c("sample", "subjid", "lesional","visit")) %>% 
  inner_join(., map, by="subjid") %>% 
  mutate(sample = str_replace_all(string=sample, pattern="non-lesional", replacement="non_lesional"),
         lesional = ifelse(lesional == "lesional", "lesional", "non_lesional"),
         lesional_visit = paste(lesional, visit, sep = "_"),
         pnrs = ifelse(visit == "V3", base, aval),
         group1 = ifelse(trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "Y", "responded treatments", 
                         ifelse(trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "N", "non-resonded treatments", "placebos")),
         lesional = factor(lesional, levels = c("non_lesional", "lesional")),
         visit = factor(visit, levels = c("V3", "V8")),
         trt01p = factor(trt01p, levels = c("Placebo", "Nemolizumab 0.5mg/kg"))) %>% 
  arrange(sample)
rownames(map_sample) = map_sample$sample

shared = shared %>% 
  filter(sample %in% map_sample$sample) %>% 
  arrange(sample)
rownames(shared) = shared$sample


## Phylotype filtering: filter rare genera
keep_otu = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(otu) %>% 
  dplyr::summarize(n = n(),
            n0 = sum(count == 0),
            n1 = n - n0,
            count = sum(count)) %>% 
  select(otu, count, n1) %>% 
  filter(n1 >= 18 & count >= 1000) %>% 
  ungroup() %>% 
  pull(otu)

shared = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  filter(otu %in% keep_otu) %>% 
  pivot_wider(names_from = otu, values_from = count) %>% 
  as.data.frame()
rownames(shared) = shared$sample

taxonomy = taxonomy %>% 
  filter(otu %in% keep_otu) %>% 
  as.data.frame()

classified_genus = taxonomy %>% 
  mutate(unclassified = str_detect(genus, pattern = "unclassified")) %>% 
  filter(!unclassified) %>% 
  pull(otu)


## abundant genus
abund_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  filter(mean >= 0.01) %>% 
  pull(otu)
abund_genus_name = taxonomy %>% 
  filter(otu %in% abund_genus) %>% 
  pull(genus)
  
