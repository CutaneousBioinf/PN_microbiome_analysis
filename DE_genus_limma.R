# This script perform genus-level DE analysis by limma-voom

library(edgeR)
library(limma)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")

## calculate normalize factors in edgeR

dge = shared %>% 
  select(-sample) %>% 
  t() %>% 
  DGEList(.)
dge = calcNormFactors(dge, method = "TMM")


## Comparing between non-lesional and lesional at baseline

map_temp = map_sample %>% filter(visit == "V3")
keep = names(which(tapply(map_temp$lesional, map_temp$subjid, function(x){length(unique(x))}) == 2))
col_data = map_sample %>% filter(subjid %in% keep & visit == "V3") %>% arrange(sample) 
keep_sample = col_data$sample
count_data = dge[,keep_sample]

lesional = factor(col_data$lesional, levels = c("non_lesional", "lesional"))
subjid = factor(col_data$subjid)
batch = factor(col_data$batch)
design =  model.matrix(~ 1 + lesional + subjid) 
voom = voom(counts = count_data, design = design, plot=T)
# corfit = duplicateCorrelation(voom, design, block = batch)
# fit = lmFit(voom, design, correlation = corfit$consensus)
fit = lmFit(voom, design)
fit = eBayes(fit)
res = topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_genus_lesional = res %>% 
  as_tibble(rownames = "otu") %>% 
  inner_join(., taxonomy) %>% 
  select(otu, logFC, t, P.Value, adj.P.Val, kingdom, phylum, class, order, family, genus)


## Comparing between baseline and week 12 in treatments

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2))
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  arrange(sample)
keep_sample = col_data$sample
count_data = dge[,keep_sample]

visit = factor(col_data$visit, levels = c("V3", "V8"))
subjid = factor(col_data$subjid)
batch = factor(col_data$batch)
design =  model.matrix(~ 1 + visit + subjid) 
voom = voom(counts = count_data, design = design, plot=T)
# corfit = duplicateCorrelation(voom, design, block = batch)
# fit = lmFit(voom, design, correlation = corfit$consensus)
fit = lmFit(voom, design)
fit = eBayes(fit)
res = topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_genus_visit_treatment = res %>% 
  as_tibble(rownames = "otu") %>% 
  inner_join(., taxonomy) %>% 
  select(otu, logFC, t, P.Value, adj.P.Val, kingdom, phylum, class, order, family, genus)


## Comparing between baseline and week 12 in placebos

map_temp = map_sample %>% filter(lesional == "lesional") ###
keep = names(which(tapply(map_temp$visit, map_temp$subjid, function(x){length(unique(x))}) == 2))
col_data = map_sample %>% 
  filter(subjid %in% keep & lesional == "lesional" & trt01p == "Placebo") %>% 
  arrange(sample)
keep_sample = col_data$sample
count_data = dge[,keep_sample]

visit = factor(col_data$visit, levels = c("V3", "V8"))
subjid = factor(col_data$subjid)
batch = factor(col_data$batch)
design =  model.matrix(~ 1 + visit + subjid) 
voom = voom(counts = count_data, design = design, plot=T)
# corfit = duplicateCorrelation(voom, design, block = batch)
# fit = lmFit(voom, design, correlation = corfit$consensus)
fit = lmFit(voom, design)
fit = eBayes(fit)
res = topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_genus_visit_placebo = res %>% 
  as_tibble(rownames = "otu") %>% 
  inner_join(., taxonomy) %>% 
  select(otu, logFC, t, P.Value, adj.P.Val, kingdom, phylum, class, order, family, genus)


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

### treatments
data = inner_join(
  de_genus_lesional %>% select(otu, logFC, adj.P.Val), 
  de_genus_visit_treatment %>% select(otu, logFC, adj.P.Val), 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$logFC.1, 
  data$logFC.2, 
  method = "spearman"
) # -0.2671902 p = 1.916e-06
data %>% 
  summarize(con.t = sum(data$logFC.1*data$logFC.2 < 0)) # 179
data %>% 
  mutate(
    sig.1 = (adj.P.Val.1 <= 0.1)*1,
    sig.2 = (adj.P.Val.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = logFC.1, y = logFC.2, 
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
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_limma_genus_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebos
data = inner_join(
  de_genus_lesional %>% select(otu, logFC, adj.P.Val), 
  de_genus_visit_placebo %>% select(otu, logFC, adj.P.Val), 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$logFC.1, 
  data$logFC.2, 
  method = "spearman"
) # -0.1923149 p = 0.000664
data %>% 
  summarize(con.t = sum(data$logFC.1*data$logFC.2 < 0)) # 188
data %>% 
  mutate(
    sig.1 = (adj.P.Val.1 <= 0.1)*1,
    sig.2 = (adj.P.Val.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = logFC.1, y = logFC.2, 
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
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_limma_genus_placebo.tiff", width = 4, height = 4, dpi = 500)
