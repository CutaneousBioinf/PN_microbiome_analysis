# This script perform phylum level DE analysis

library(edgeR)
library(limma)

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
## calculate normalize factors in edgeR
dge = shared_phylum %>% 
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
de_phylum_lesional = res %>% 
  as_tibble(rownames = "phylum") 


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
de_phylum_visit_treatment = res %>% 
  as_tibble(rownames = "phylum") 


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
de_phylum_visit_placebo = res %>% 
  as_tibble(rownames = "phylum") 

rbind(de_phylum_lesional %>% 
        select(phylum, logFC, adj.P.Val) %>% 
        mutate(var = "l"), 
      de_phylum_visit_treatment %>% 
        select(phylum, logFC, adj.P.Val) %>% 
        mutate(var = "t")) %>% 
  rbind(., 
        de_phylum_visit_placebo %>% 
          select(phylum, logFC, adj.P.Val) %>% 
          mutate(var = "p")) %>% 
  mutate(sig = ifelse(adj.P.Val <= 0.1, "Y", "N")) %>% 
  ggplot(aes(x = var, y = logFC, fill = var, color = sig)) +
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
  theme_classic() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 13),
        legend.position = "bottom")


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

### treatments
data = inner_join(
  de_phylum_lesional %>% select(phylum, logFC, adj.P.Val), 
  de_phylum_visit_treatment %>% select(phylum, logFC, adj.P.Val), 
  by = "phylum", suffix = c(".1", ".2")
) 
cor.test(
  data$logFC.1, 
  data$logFC.2, 
  method = "spearman"
) # -0.7142857 p =  0.003806
data %>% 
  summarize(con.t = sum(data$logFC.1*data$logFC.2 < 0)) # 9
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
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_limma_phylum_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebos
data = inner_join(
  de_phylum_lesional %>% select(phylum, logFC, adj.P.Val), 
  de_phylum_visit_placebo %>% select(phylum, logFC, adj.P.Val), 
  by = "phylum", suffix = c(".1", ".2")
) 
cor.test(
  data$logFC.1, 
  data$logFC.2, 
  method = "spearman"
) # 0.1 p = 0.7241
data %>% 
  summarize(con.t = sum(data$logFC.1*data$logFC.2 < 0)) # 5
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
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +    
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_limma_phylum_placebo.tiff", width = 4, height = 4, dpi = 500)



