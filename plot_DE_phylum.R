## This script plots scatter plot of phylum level logFC from non-lesional to lesional at baseline and logFC from baseline to week 12 in lesional skin 

library(tibble)
library(forcats)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpattern)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")


## Reading results

de_phylum_lesional = read.csv(file = "ANCOMBC results/de_phylum_lesional.csv")
de_phylum_visit_treatment.responded = read.csv(file = "ANCOMBC results/de_phylum_visit_treatment.responded.csv")
de_phylum_visit_treatment.not.responded = read.csv(file = "ANCOMBC results/de_phylum_visit_treatment.not.responded.csv")
de_phylum_visit_treatment = read.csv(file = "ANCOMBC results/de_phylum_visit_treatment.csv")
de_phylum_visit_placebo = read.csv(file = "ANCOMBC results/de_phylum_visit_placebo.csv")


## log fold change correlation with lesional v.s. non-lesional at baseline (treatments/ placebos)

### responded treatments
data = inner_join(
  de_phylum_lesional, 
  de_phylum_visit_treatment.responded, 
  by = "taxon", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.536 p = 0.042
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 12
data %>% 
  mutate(
    sig.1 = (p.adj.1 <= 0.1)*1,
    sig.2 = (p.adj.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = lfc.1, y = lfc.2, color = as.character(sig.1), 
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
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_phylum_treatment.responded.tiff", width = 4, height = 4, dpi = 500)

### not responded treatments
data = inner_join(
  de_phylum_lesional, 
  de_phylum_visit_treatment.not.responded, 
  by = "taxon", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.536 p = 0.0422
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 11
data %>% 
  mutate(
    sig.1 = (p.adj.1 <= 0.1)*1,
    sig.2 = (p.adj.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = lfc.1, y = lfc.2, color = as.character(sig.1), 
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
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_phylum_treatment.not.responded.tiff", width = 4, height = 4, dpi = 500)

### treatments
data = inner_join(
  de_phylum_lesional, 
  de_phylum_visit_treatment, 
  by = "taxon", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.6357 p = 0.01287
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 12
data %>% 
  mutate(
    sig.1 = (p.adj.1 <= 0.1)*1,
    sig.2 = (p.adj.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = lfc.1, y = lfc.2, color = as.character(sig.1), 
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
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_phylum_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebos
data = inner_join(
  de_phylum_lesional, 
  de_phylum_visit_placebo, 
  by = "taxon", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.0321 p = 0.9132
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 10
data %>% 
  mutate(
    sig.1 = (p.adj.1 <= 0.1)*1,
    sig.2 = (p.adj.2 <= 0.1)*1,
    group = sig.1 + sig.2
  )  %>%
  ggplot(aes(
    x = lfc.1, y = lfc.2, color = as.character(sig.1), 
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
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  scale_shape_manual(breaks = c("0", "1"), values = c(16, 17)) +
  scale_color_manual(breaks = c("0", "1"), values = c("grey", "black")) +
  scale_size_manual(breaks = c("0", "1", "2"), values = c(2, 3, 5)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 13, face = "bold", family = "Times")
  )
ggsave("manuscript figures/scatter.logfc_phylum_placebo.tiff", width = 4, height = 4, dpi = 500)

rbind(
  de_phylum_lesional %>% mutate(var = "l") %>% rownames_to_column("otu"), 
  de_phylum_visit_treatment %>% mutate(var = "t") %>% rownames_to_column("otu")
) %>% 
  rbind(., de_phylum_visit_placebo %>% mutate(var = "p") %>% rownames_to_column("otu")) %>% 
  mutate(
    sig = (p.adj <= 0.1)*1,
    taxon = str_replace_all(taxon, "Bacteria_unclassified", "Unclassified"),
    taxon = str_replace_all(taxon, "Candidatus_Saccharibacteria", "Candi._Saccha."),
    taxon = str_replace_all(taxon, "Cyanobacteria/Chloroplast", "Cyano./Chloro."),
    taxon = str_replace_all(taxon, "Deinococcus-Thermus", "Deino.-Thermus"),
    taxon = str_replace_all(taxon, "Gemmatimonadetes", "Gemma."),
    taxon = fct_reorder(taxon, as.numeric(otu))
  ) %>% 
  ggplot(aes(x = var, y = lfc, fill = var, pattern = as.factor(sig))) +
  geom_bar_pattern(
    stat = "identity", 
    position = "dodge", 
    pattern_angle = 45,
    pattern_density = 0.001
    ) +
  geom_abline(intercept = 0, slope = 0) +
  facet_wrap(~ taxon, scales = "free") +
  scale_fill_manual(
    name = NULL,
    breaks = c("l", "p", "t"),
    labels = c(
      "non-lesional vs lesional", 
      "baseline vs week 12 in placebo", 
      "baseline vs week 12 in treatment"
      ),
    values = c("grey", "blue", "red")) +
  scale_pattern_manual(
    name = NULL,
    breaks = c("0", "1"), labels = c("Non Significant", "Significant"),
    values = c("none", "stripe"), guide = "none"
    ) +
  labs(title = NULL,x = NULL, y = "log(Fold Change)") +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE)
    ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(size = 50, face = "bold", family = "Times"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 13, family = "Times"),
    axis.text.x = element_blank(),
    axis.line = element_line(),
    legend.text = element_text(size = 30, face = "bold", family = "Times"),
    strip.text = element_text(size = 28, face = "bold.italic", family = "Times"),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )
ggsave("manuscript figures/bar.logfc_phylum_nonles.les_treatment_placebo.pdf", width = 15, height = 12, dpi = 500)
