# This script preprocess mothur-formatted data, including phylotype filtering and an overview for metadata

library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")

## Import data as data.frame

### metadata
meta = read.table("dat/metaDataGaldermaNoOutlirs.txt", header = T) %>% 
  dplyr::rename("subjid" = "Subject") %>% 
  group_by(subjid) %>% 
  summarize(
    age = mean(Age)
  ) %>% 
  mutate(subjid = as.character(subjid))
location = read.csv(file = "dat/115828_238637.csv") %>% 
  mutate(
    subjid = gsub(pattern = "-", replacement = "", x = Subject)
  ) %>% 
  select(subjid, SQLA, SQNLA, WDSPL, NPLA, OTHLESAR) %>% 
  filter(SQLA != "") %>% 
  mutate(
    W12SQLA = ifelse(NPLA == "", SQLA, NPLA)
  )
map = read_excel(path = "dat/PN2_subject_annotation.xlsx") %>% 
  rename_all(tolower) %>% 
  mutate(subjid = as.character(subjid)) %>% 
  filter(! subjid %in% c(5433004, 5471012, 5891001)) %>% 
  inner_join(., location, by = "subjid") 

### taxonomy data
taxonomy = read.table(file = "dat/galderma.tx.1.cons.taxonomy", header = T) %>% 
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(string = taxonomy, pattern="\\(\\d*\\)", replacement = "")) %>%
  mutate(taxonomy = str_replace_all(string = taxonomy, pattern=";$", replacement = "")) %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")
chloroplast_otu = taxonomy %>% 
  filter(str_detect(phylum, "Chloroplast")) %>% 
  pull(otu)
chloroplast_size = taxonomy %>% 
  filter(str_detect(phylum, "Chloroplast")) %>% 
  group_by(.) %>% 
  summarize(sum = sum(size))
chloroplast_size / sum(taxonomy$size) # 0.0007567522; 0.076%

### removal of Chloroplast sequences
taxonomy = taxonomy %>% 
  filter(!otu %in% chloroplast_otu)

### otu table
shared = read.table(file = "dat/galderma.tx.shared", header=T) %>% 
  dplyr::rename(sample = Group) %>% 
  select(-label, -numPhylos) %>% 
  select(-all_of(chloroplast_otu))

### sample meta data 
map_sample = shared %>%
  select(sample) %>%
  mutate(sample = str_replace_all(string=sample, pattern="non_lesional", replacement="non-lesional")) %>%
  separate(col="sample", into = c("unique", "subjid", "lesional","visit"), sep="_", remove=FALSE) %>% 
  select(c("sample", "subjid", "lesional","visit")) %>% 
  inner_join(., map, by="subjid") %>% 
  mutate(
    sample = str_replace_all(string=sample, pattern="non-lesional", replacement="non_lesional"),
    lesional = ifelse(lesional == "lesional", "lesional", "non_lesional"),
    lesional_visit = paste(lesional, visit, sep = "_"),
    pnrs = ifelse(visit == "V3", base, aval),
    sites = case_when(
      lesional_visit == "lesional_V3" ~ SQLA, 
      lesional_visit == "non_lesional_V3" ~ SQNLA,
      lesional_visit == "lesional_V8" ~ W12SQLA
    ),
    group1 = ifelse(
      trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "Y", "responded treatments", 
      ifelse(trt01p == "Nemolizumab 0.5mg/kg" & nrs_w12 == "N", "non-resonded treatments", "placebos")),
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    visit = factor(visit, levels = c("V3", "V8")),
    trt01p = factor(trt01p, levels = c("Placebo", "Nemolizumab 0.5mg/kg"))
    ) %>% 
  arrange(sample)

### batch data
run1 = read.table(file = "dat/run1_fastq-gz_md5sum.md5", header = FALSE) %>% 
  mutate(V2 = str_replace_all(V2, pattern = "non-lesional", replacement = "non_lesional")) %>% 
  separate(V2, into = c("unique", "subjid", "lesional", "visit"), sep = "-", remove = FALSE) %>% 
  drop_na() %>% 
  mutate(visit = substr(visit, start = 1, stop = 2), 
         sample = paste(unique, subjid, lesional, visit, sep = "_")) %>% 
  select(sample) %>%
  mutate(batch = 1) %>% 
  group_by(sample) %>% 
  summarize(batch = max(batch))
run2 = read.table(file = "dat/run2_fastq-gz_md5sum.md5", header = FALSE) %>% 
  mutate(V2 = str_replace_all(V2, pattern = "non-lesional", replacement = "non_lesional")) %>% 
  separate(V2, into = c("unique", "subjid", "lesional", "visit"), sep = "-", remove = FALSE) %>% 
  drop_na() %>% 
  mutate(visit = substr(visit, start = 1, stop = 2), 
         sample = paste(unique, subjid, lesional, visit, sep = "_")) %>% 
  select(sample) %>%
  mutate(batch = 2) %>% 
  group_by(sample) %>% 
  summarize(batch = max(batch))
run3 = read.table(file = "dat/run3_fastq-gz_md5sum.md5", header = FALSE) %>% 
  mutate(V2 = str_replace_all(V2, pattern = "non-lesional", replacement = "non_lesional")) %>% 
  separate(V2, into = c("unique", "subjid", "lesional", "visit"), sep = "-", remove = FALSE) %>% 
  drop_na() %>% 
  mutate(visit = substr(visit, start = 1, stop = 2), 
         sample = paste(unique, subjid, lesional, visit, sep = "_")) %>% 
  select(sample) %>%
  mutate(batch = 3) %>% 
  group_by(sample) %>% 
  summarize(batch = max(batch))
batch = rbind(run1, run2, run3)
rm("run1", "run2", "run3")
map_sample = inner_join(map_sample, batch, by = "sample")
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

### rationale for threshold
shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(otu) %>% 
  dplyr::summarize(
    n = n(),
    n0 = sum(count == 0),
    n1 = n - n0,
    count = sum(count)
    ) %>% 
  select(otu, count, n1) %>% 
  arrange(count) %>% 
  ggplot(aes(x = 1:1582, y = count)) +
  geom_point() +
  geom_line() +
  coord_cartesian(ylim = c(0, 2000)) +
  labs(x = "Genus", y = "Counts") +
  theme_bw()
shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(otu) %>% 
  dplyr::summarize(
    n = n(),
    n0 = sum(count == 0),
    n1 = n - n0,
    count = sum(count)
  ) %>% 
  select(otu, count, n1) %>% 
  arrange(n1) %>% 
  ggplot(aes(y = 1:1582, x = n1)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 18, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 20)) +
  labs(y = "Genus", x = "# of sample occurence in 159 samples") +
  theme_bw()
t1 = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(abundance = count/sum(count)) %>% 
  group_by(otu) %>% 
  summarize(abundance = mean(abundance)) %>% 
  filter(otu %in% keep_otu)
summary(t1$abundance)
t2 = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(abundance = count/sum(count)) %>% 
  group_by(otu) %>% 
  summarize(abundance = mean(abundance)) %>% 
  filter(!otu %in% keep_otu)
summary(t2$abundance)

### Filtering
shared = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  filter(otu %in% keep_otu) %>% 
  pivot_wider(names_from = otu, values_from = count) %>% 
  as.data.frame()
rownames(shared) = shared$sample
taxonomy = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(otu) %>% 
  summarize(size = sum(count)) %>% 
  inner_join(., taxonomy %>% select(-size), by = "otu") %>% 
  as.data.frame()
classified_genus = taxonomy %>% 
  mutate(unclassified = str_detect(genus, pattern = "unclassified")) %>% 
  filter(!unclassified) %>% 
  pull(otu)

## common genus

common_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  filter(mean >= 0.01) %>% 
  pull(otu)
common_genus_name = taxonomy %>% 
  filter(otu %in% common_genus) %>% 
  pull(genus)

