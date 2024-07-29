# This script perform genus-level DE analysis by DESeq2 

library(DESeq2)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")

## non-lesional v.s. lesional areas at baseline, controlling for individual effect

map_temp = map_sample %>% filter(visit == "V3")
keep = names(which(tapply(map_temp$lesional, map_temp$subjid, function(x){length(unique(x))}) == 2))
count_data = shared %>% 
  mutate(sample = str_replace_all(string=sample, pattern="non_lesional", replacement="non-lesional")) %>%
  separate(col="sample", into = c("unique", "subjid", "lesional","visit"), sep="_", remove=FALSE) %>% 
  filter(subjid %in% keep & visit == "V3") %>% 
  arrange(sample) %>% 
  select(-sample, -unique, -subjid, -lesional, -visit) %>% 
  t()
col_data = map_sample %>% filter(subjid %in% keep & visit == "V3") %>% arrange(sample)  

dds = DESeqDataSetFromMatrix(
  countData = count_data, colData = col_data, 
  design = ~ lesional + subjid
  )
dds$lesional = relevel(dds$lesional, ref = "non_lesional")
dds = DESeq(dds)
res = results(dds, name = "lesional_lesional_vs_non_lesional", pAdjustMethod="fdr" )
de_genus_lesional = res@listData %>% 
  as.data.frame() %>% 
  mutate(otu = res@rownames) %>% 
  inner_join(., taxonomy) %>% 
  select(otu, log2FoldChange, stat, pvalue, padj, kingdom, phylum, class, order, family, genus)


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in treatment groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% ### 
  arrange(sample) %>% 
  select(-sample, -subjid, -lesional, -trt01p) %>% 
  t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  arrange(sample) ###

dds = DESeqDataSetFromMatrix(
  countData = count_data, colData = col_data, 
  design = ~ visit + subjid
  )
dds = DESeq(dds)
res = results(dds, name = "visit_V8_vs_V3", pAdjustMethod="fdr" )
de_genus_visit_treatment = res@listData %>% 
  as.data.frame() %>% 
  mutate(otu = res@rownames) %>% 
  inner_join(., taxonomy) %>% 
  select(otu, log2FoldChange, lfcSE, stat, pvalue, padj, kingdom, phylum, class, order, family, genus)


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in placebo groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>%   ### 
  arrange(sample) %>% 
  select(-sample, -subjid, -lesional, -trt01p) %>% 
  t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% ###
  arrange(sample)
  
dds = DESeqDataSetFromMatrix(
  countData = count_data, colData = col_data,
  design = ~ visit + subjid
  )
dds = DESeq(dds)
res = results(dds, name = "visit_V8_vs_V3", pAdjustMethod="fdr" )
de_genus_visit_placebo = res@listData %>% 
  as.data.frame() %>% 
  mutate(otu = res@rownames) %>% 
  inner_join(., taxonomy) %>% 
  select(otu, log2FoldChange, lfcSE, stat, pvalue, padj, kingdom, phylum, class, order, family, genus)


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

### treatments
data = inner_join(
  de_genus_lesional %>% select(otu, log2FoldChange, padj), 
  de_genus_visit_treatment %>% select(otu, log2FoldChange, padj), 
  by = "otu", suffix = c(".1", ".2")
  ) 
cor.test(
  data$log2FoldChange.1, 
  data$log2FoldChange.2, 
  method = "spearman"
  ) # -0.181157 p = 0.001356
data %>% 
  summarize(con.t = sum(log2FoldChange.1*log2FoldChange.2 < 0)) # 178
data %>% 
    mutate(
      sig.1 = (padj.1 <= 0.1)*1,
      sig.2 = (padj.2 <= 0.1)*1,
      group = sig.1 + sig.2
    )  %>%
    ggplot(aes(
      x = log2FoldChange.1, y = log2FoldChange.2, 
      color = as.character(sig.1), 
      shape = as.character(sig.2),
      size = as.character(group)
    )) +
    geom_point() +
    geom_abline(intercept = 0, slope = -1, lty = 2) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    labs(
      title = NULL,
      x = "non-lesional vs\nlesional at baseline", 
      y = "baseline vs week-\n12 in lesional skin"
      ) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
    scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
    scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
    scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18, face = "bold", family = "Times"),
      axis.text = element_text(size = 13, face = "bold", family = "Times")
      )
ggsave("manuscript figures/scatter.logfc_deseq2_genus_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebos
data = inner_join(
  de_genus_lesional %>% select(otu, log2FoldChange, padj), 
  de_genus_visit_placebo %>% select(otu, log2FoldChange, padj), 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$log2FoldChange.1, 
  data$log2FoldChange.2, 
  method = "spearman"
) # -0.16408 p = 0.003748
data %>% 
  summarize(con.t = sum(log2FoldChange.1*log2FoldChange.2 < 0)) # 163
data %>% 
  mutate(
    sig.1 = (padj.1 <= 0.1)*1,
    sig.2 = (padj.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = log2FoldChange.1, y = log2FoldChange.2, 
    color = as.character(sig.1), 
    shape = as.character(sig.2),
    size = as.character(group)
  )) +
  geom_point() +
  geom_abline(intercept = 0, slope = -1, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  labs(
    title = NULL,
    x = "non-lesional vs\nlesional at baseline", 
    y = "baseline vs week-\n12 in lesional skin"
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_deseq2_genus_placebo.tiff", width = 4, height = 4, dpi = 500)





