## This script plots scatter plot of genus-level logFC from non-lesional to lesional at baseline and logFC from baseline to week 12 in lesional skin 

library(tibble)
library(forcats)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpattern)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")


## Reading data

de_genus_lesional = read.csv(file = "ANCOMBC results/de_genus_lesional.csv")
de_genus_visit_placebo = read.csv(file = "ANCOMBC results/de_genus_visit_placebo.csv")
de_genus_visit_treatment.responded = read.csv(file = "ANCOMBC results/de_genus_visit_treatment.responded.csv")
de_genus_visit_treatment.not.responded = read.csv(file = "ANCOMBC results/de_genus_visit_treatment.not.responded.csv")
de_genus_visit_treatment = read.csv(file = "ANCOMBC results/de_genus_visit_treatment.csv")


######### Including all genus #########

### responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.responded, 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.206 p = 0.0002985
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 183
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
ggsave("manuscript figures/scatter.logfc_genus_treatment.responded.tiff", width = 4, height = 4, dpi = 500)

### not responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.not.responded, 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.19235 p = 0.0007477
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 160
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
ggsave("manuscript figures/scatter.logfc_genus_treatment.not.responded.tiff", width = 4, height = 4, dpi = 500)

### treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment, 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.2927835 p = 2.108e-07
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 180
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
ggsave("manuscript figures/scatter.logfc_genus_treatment.tiff", width = 4, height = 4, dpi = 500)


### placebo
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_placebo, 
  by = "otu", suffix = c(".1", ".2")
) 
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.1576817 p = 0.005625
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 165
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
ggsave("manuscript figures/scatter.logfc_genus_placebo.tiff", width = 4, height = 4, dpi = 500)


######### Including top 5 genus #########

top_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  top_n(mean, n = 5) %>% 
  pull(otu)

### responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.3 p = 0.6833
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 2
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
ggsave("manuscript figures/scatter.logfc_genus.top5_treatment.responded.tiff", width = 4, height = 4, dpi = 500)

### not responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.not.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.1 p = 0.95
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 1
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
ggsave("manuscript figures/scatter.logfc_genus.top5_treatment.not.responded.tiff", width = 4, height = 4, dpi = 500)


### treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.3 p = 0.6833
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 3
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
ggsave("manuscript figures/scatter.logfc_genus.top5_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebo
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_placebo, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.1 p = 0.95
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 4
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
ggsave("manuscript figures/scatter.logfc_genus.top5_placebo.tiff", width = 4, height = 4, dpi = 500)

rbind(
  de_genus_lesional %>% mutate(var = "l"), 
  de_genus_visit_treatment %>% mutate(var = "t")
) %>% 
  rbind(., de_genus_visit_placebo %>% mutate(var = "p")) %>% 
  filter(otu %in% top_genus) %>%
  mutate(
    genus = fct_reorder(genus, as.numeric(factor(otu))),
    sig = (p.adj <= 0.1)*1
  ) %>% 
  ggplot(aes(x = var, y = lfc, fill = var, pattern = as.factor(sig))) +
  geom_bar_pattern(
    stat = "identity", 
    position = "dodge", 
    pattern_angle = 45,
    pattern_density = 0.001
  ) +
  geom_abline(intercept = 0, slope = 0) +
  facet_wrap(~ genus, scales = "free") +
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
ggsave("manuscript figures/bar.logfc_genus.top5_nonles.les_treatment_placebo.pdf", width = 15, height = 12, dpi = 500)

rbind(
  de_genus_lesional %>% mutate(var = "l"), 
  de_genus_visit_treatment.responded %>% mutate(var = "tr")
) %>% 
  rbind(., de_genus_visit_treatment.not.responded %>% mutate(var = "tnr")) %>% 
  rbind(., de_genus_visit_placebo %>% mutate(var = "p")) %>% 
  filter(otu %in% top_genus) %>%
  mutate(
    genus = fct_reorder(genus, as.numeric(factor(otu))),
    sig = (p.adj <= 0.1)*1
  ) %>% 
  ggplot(aes(x = var, y = lfc, fill = var, pattern = as.factor(sig))) +
  geom_bar_pattern(
    stat = "identity", 
    position = "dodge", 
    pattern_angle = 45,
    pattern_density = 0.001
  ) +
  geom_abline(intercept = 0, slope = 0) +
  facet_wrap(~ genus, scales = "free") +
  scale_fill_manual(
    name = NULL,
    breaks = c("l", "p", "tnr", "tr"),
    labels = c(
      "non-lesional vs lesional", 
      "baseline vs week 12 in placebo", 
      "baseline vs week 12 in non-responded treatment",
      "baseline vs week 12 in responded treatment"
    ),
    values = c("grey", "blue", "yellow", "red")) +
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
ggsave("manuscript figures/bar.logfc_genus.top5_nonles.les_responded.treatment_not.responded.treatment_placebo.pdf", width = 15, height = 12, dpi = 500)


######### Including top 10 genus #########

top_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  top_n(mean, n = 10) %>% 
  pull(otu)

### responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # 0.054545 p = 0.8916
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 6
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
ggsave("manuscript figures/scatter.logfc_genus.top10_treatment.responded.tiff", width = 4, height = 4, dpi = 500)

### not responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.not.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.455 p = 0.1909
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 3
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
ggsave("manuscript figures/scatter.logfc_genus.top10_treatment.not.responded.tiff", width = 4, height = 4, dpi = 500)

### treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.321 p = 0.3677
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 5
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
ggsave("manuscript figures/scatter.logfc_genus.top10_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebo
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_placebo, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.0788 p = 0.838
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 6
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
ggsave("manuscript figures/scatter.logfc_genus.top10_placebo.tiff", width = 4, height = 4, dpi = 500)

rbind(
  de_genus_lesional %>% mutate(var = "l"), 
  de_genus_visit_treatment %>% mutate(var = "t")
) %>% 
  rbind(., de_genus_visit_placebo %>% mutate(var = "p")) %>% 
  filter(otu %in% top_genus) %>%
  mutate(
    genus = fct_reorder(genus, as.numeric(factor(otu))),
    sig = (p.adj <= 0.1)*1
  ) %>% 
  ggplot(aes(x = var, y = lfc, fill = var, pattern = as.factor(sig))) +
  geom_bar_pattern(
    stat = "identity", 
    position = "dodge", 
    pattern_angle = 45,
    pattern_density = 0.001
  ) +
  geom_abline(intercept = 0, slope = 0) +
  facet_wrap(~ genus, scales = "free") +
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
ggsave("manuscript figures/bar.logfc_genus.top10_nonles.les_treatment_placebo.pdf", width = 15, height = 12, dpi = 500)

rbind(
  de_genus_lesional %>% mutate(var = "l"), 
  de_genus_visit_treatment.responded %>% mutate(var = "tr")
) %>% 
  rbind(., de_genus_visit_treatment.not.responded %>% mutate(var = "tnr")) %>% 
  rbind(., de_genus_visit_placebo %>% mutate(var = "p")) %>% 
  filter(otu %in% top_genus) %>%
  mutate(
    genus = fct_reorder(genus, as.numeric(factor(otu))),
    sig = (p.adj <= 0.1)*1
  ) %>% 
  ggplot(aes(x = var, y = lfc, fill = var, pattern = as.factor(sig))) +
  geom_bar_pattern(
    stat = "identity", 
    position = "dodge", 
    pattern_angle = 45,
    pattern_density = 0.001
  ) +
  geom_abline(intercept = 0, slope = 0) +
  facet_wrap(~ genus, scales = "free") +
  scale_fill_manual(
    name = NULL,
    breaks = c("l", "p", "tnr", "tr"),
    labels = c(
      "non-lesional vs lesional", 
      "baseline vs week 12 in placebo", 
      "baseline vs week 12 in non-responded treatment",
      "baseline vs week 12 in responded treatment"
    ),
    values = c("grey", "blue", "yellow", "red")) +
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
ggsave("manuscript figures/bar.logfc_genus.top10_nonles.les_responded.treatment_not.responded.treatment_placebo.pdf", width = 15, height = 12, dpi = 500)


######### Including top 11 genus (with mean relative abundance >= 0.01) #########

top_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  top_n(mean, n = 11) %>% 
  pull(otu)

### responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # 0.127 p = 0.7138
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 7
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
ggsave("manuscript figures/scatter.logfc_genus.top11_treatment.responded.tiff", width = 4, height = 4, dpi = 500)


### not responded treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment.not.responded, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.318, p = 0.3414
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 4
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
ggsave("manuscript figures/scatter.logfc_genus.top11_treatment.not.responded.tiff", width = 4, height = 4, dpi = 500)


### treatments
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_treatment, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.163 p = 0.6339
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 6
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
ggsave("manuscript figures/scatter.logfc_genus.top11_treatment.tiff", width = 4, height = 4, dpi = 500)

### placebo
data = inner_join(
  de_genus_lesional, 
  de_genus_visit_placebo, 
  by = "otu", suffix = c(".1", ".2")
) %>% 
  filter(otu %in% top_genus)
cor.test(
  data$lfc.1, data$lfc.2, method = "spearman"
) # -0.1818 p = 0.5952
data %>% 
  summarize(con.t = sum(lfc.1*lfc.2 < 0)) # 6
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
ggsave("manuscript figures/scatter.logfc_genus.top11_placebo.tiff", width = 4, height = 4, dpi = 500)
