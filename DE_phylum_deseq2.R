# This script perform phylum level DE analysis

library(DESeq2)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")

## We first summarize shared file and taxonomy file to phylum level

shared_phylum = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., taxonomy %>% select(otu, phylum)) %>% 
  group_by(sample, phylum) %>% 
  summarize(count = sum(count)) %>% 
  pivot_wider(names_from = "phylum", values_from = "count") %>% 
  as.data.frame()
rownames(shared_phylum) = shared_phylum$sample


## non-lesional v.s. lesional areas at baseline, controlling for individual effect

map_temp = map_sample %>% filter(visit == "V3")
keep = names(which(tapply(map_temp$lesional, map_temp$subjid, function(x){length(unique(x))}) == 2))
count_data = shared_phylum %>% 
  mutate(sample = str_replace_all(string=sample, pattern="non_lesional", replacement="non-lesional")) %>%
  separate(col="sample", into = c("unique", "subjid", "lesional","visit"), sep="_", remove=FALSE) %>% 
  filter(subjid %in% keep & visit == "V3") %>% 
  arrange(sample) %>% 
  select(-sample, -unique, -subjid, -lesional, -visit) %>% 
  t()
col_data = map_sample %>% filter(subjid %in% keep & visit == "V3") %>% arrange(sample)

dds = DESeqDataSetFromMatrix(countData = count_data, 
                             colData = col_data,
                             design = ~ lesional + subjid)
dds$lesional = relevel(dds$lesional, ref = "non_lesional")
dds = DESeq(dds)
res = results(dds, name = "lesional_lesional_vs_non_lesional", pAdjustMethod="fdr" )
de_phylum_lesional = res@listData %>% 
  as.data.frame() %>% 
  mutate(phylum = res@rownames) 


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in treatment groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared_phylum %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% ### 
  arrange(sample) %>% 
  select(-sample, -subjid, -lesional, -trt01p) %>% 
  t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% ###
  arrange(sample)

dds = DESeqDataSetFromMatrix(countData = count_data, 
                             colData = col_data,
                             design = ~ visit + subjid)
dds = DESeq(dds)
res = results(dds, name = "visit_V8_vs_V3", pAdjustMethod="fdr" )
de_phylum_visit_treatment = res@listData %>% 
  as.data.frame() %>% 
  mutate(phylum = res@rownames)


## baseline v.s. week 12 at lesional areas, controlling for individual effect, in placebo groups

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2)) ###
count_data = shared_phylum %>% 
  inner_join(., map_sample %>% select(sample, subjid, trt01p, lesional)) %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% ### 
  arrange(sample) %>% 
  select(-sample, -subjid, -lesional, -trt01p) %>% 
  t()
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% ###
  arrange(sample)

dds = DESeqDataSetFromMatrix(countData = count_data, 
                             colData = col_data,
                             design = ~ visit + subjid)
dds = DESeq(dds)
res = results(dds, name = "visit_V8_vs_V3", pAdjustMethod="fdr" )
de_phylum_visit_placebo = res@listData %>% 
  as.data.frame() %>% 
  mutate(phylum = res@rownames)

rbind(de_phylum_lesional %>% 
        select(phylum, log2FoldChange, padj) %>% 
        mutate(var = "l"), 
      de_phylum_visit_treatment %>% 
        select(phylum, log2FoldChange, padj) %>% 
        mutate(var = "t")) %>% 
  rbind(., 
        de_phylum_visit_placebo %>% 
          select(phylum, log2FoldChange, padj) %>% 
          mutate(var = "p")) %>% 
  mutate(sig = ifelse(padj <= 0.1, "Y", "N")) %>% 
  ggplot(aes(x = var, y = log2FoldChange, fill = var, color = sig)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ phylum, scales = "free") +
  scale_fill_manual(name = NULL,
                    breaks = c("l", "p", "t"),
                    labels = c("non-lesional vs lesional",
                               "baseline vs week 12 in placebo",
                               "baseline vs week 12 in treatment"),
                    values = c("grey", "blue", "red")) +
  scale_color_manual(name = NULL,
                     breaks = c("Y", "N"),
                     values = c("black", "white")) +
  labs(title = NULL, 
       x = NULL,
       y = "log-fold change") +
  guides(color = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 25, face = "bold", family = "Times"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 13, family = "Times"),
        axis.line = element_line(),
        legend.text = element_text(size = 25, face = "bold", family = "Times"),
        strip.text = element_text(size = 18, face = "bold.italic", family = "Times"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        )


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

### treatments
data = inner_join(
  de_phylum_lesional %>% select(phylum, log2FoldChange, padj), 
  de_phylum_visit_treatment %>% select(phylum, log2FoldChange, padj), 
  by = "phylum", suffix = c(".1", ".2")
) 
cor.test(
  data$log2FoldChange.1, 
  data$log2FoldChange.2, 
  method = "spearman"
) # -0.6607143  p = 0.009043
data %>% 
  summarize(con.t = sum(log2FoldChange.1*log2FoldChange.2 < 0)) # 14
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
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_deseq2_phylumn_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebos
data = inner_join(
  de_phylum_lesional %>% select(phylum, log2FoldChange, padj), 
  de_phylum_visit_placebo %>% select(phylum, log2FoldChange, padj), 
  by = "phylum", suffix = c(".1", ".2")
) 
cor.test(
  data$log2FoldChange.1, 
  data$log2FoldChange.2, 
  method = "spearman"
) # 0.2321429 p = 0.4039
data %>% 
  summarize(con.t = sum(log2FoldChange.1*log2FoldChange.2 < 0)) # 7
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
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_deseq2_phylum_placebo.tiff", width = 4, height = 4, dpi = 500)



