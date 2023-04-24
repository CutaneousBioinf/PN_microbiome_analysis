library(phyloseq)
library(ANCOMBC)

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
res = ancombc(data_phyloseq, formula = "lesional + subjid", tax_level = "genus",
              group = NULL, p_adj_method = "fdr", max_iter = 100)
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
write_csv(de_genus_lesional, file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_lesional.csv")



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
res = ancombc(data_phyloseq, formula = "visit + subjid", tax_level = "genus",
              group = NULL, p_adj_method = "fdr", max_iter = 100)
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
write_csv(de_genus_visit_treatment, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_visit_treatment.csv")


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
res = ancombc(data_phyloseq, formula = "visit + subjid", tax_level = "genus",
              group = NULL, p_adj_method = "fdr", max_iter = 100)
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
write_csv(de_genus_visit_placebo, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_visit_placebo.csv")


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

de_genus_lesional = read_csv(file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_lesional.csv")
de_genus_visit_treatment = read_csv(file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_visit_treatment.csv")
de_genus_visit_placebo = read_csv(file = "/Users/zrayw/Desktop/PN_microbiome/analysis/ANCOMBC results/de_genus_visit_placebo.csv")

### treatments

visit_lesional_treatment = inner_join(de_genus_lesional, de_genus_visit_treatment, 
                                      by = c("otu", "genus"), suffix = c(".1", ".2")) 

cor.test(visit_lesional_treatment$lfc.1, 
         visit_lesional_treatment$lfc.2, 
         method = "spearman")

visit_lesional_treatment %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0))

visit_lesional_treatment %>% 
  mutate(sig.1 = (p.adj.1 <= 0.1)*1,
         sig.2 = (p.adj.2 <= 0.1)*1,
         group = sig.1 + sig.2)  %>%
  ggplot(aes(x = lfc.1, y = lfc.2, 
             color = as.character(sig.1), shape = as.character(sig.2),
             size = as.character(group))) +
  geom_point() +
  geom_abline(intercept = 0, slope = -1, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  labs(title = NULL,
       x = "non-lesional vs\nlesional at baseline", 
       y = "baseline vs week-\n12 in lesional skin") +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  scale_shape_manual(breaks = c("0", "1"),
                     values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"),
                     values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"),
                    values = c(2, 3, 5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 13, face = "bold", family = "Times"))
ggsave("/Users/zrayw/Desktop/PN_microbiome/analysis/manuscript figures/de_genus_1.tiff", width = 4, height = 4, dpi = 500)


### placebos

visit_lesional_placebo = inner_join(de_genus_lesional, de_genus_visit_placebo, 
                                    by = c("otu", "genus"), suffix = c(".1", ".2")) 

cor.test(visit_lesional_placebo$lfc.1, 
         visit_lesional_placebo$lfc.2, 
         method = "spearman")

visit_lesional_placebo %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0))

visit_lesional_placebo %>% 
  mutate(sig.1 = (p.adj.1 <= 0.1)*1,
         sig.2 = (p.adj.2 <= 0.1)*1,
         group = sig.1 + sig.2)  %>%
  ggplot(aes(x = lfc.1, y = lfc.2, 
             color = as.character(sig.1), shape = as.character(sig.2),
             size = as.character(group))) +
  geom_point() +
  geom_abline(intercept = 0, slope = -1, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  labs(title = NULL,
       x = "non-lesional vs\nlesional at baseline", 
       y = "baseline vs week-\n12 in lesional skin") +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  scale_shape_manual(breaks = c("0", "1"),
                     values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"),
                     values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"),
                    values = c(2, 3, 5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 13, face = "bold", family = "Times"))
ggsave("/Users/zrayw/Desktop/PN_microbiome/analysis/manuscript figures/de_genus_2.tiff", width = 4, height = 4, dpi = 500)

