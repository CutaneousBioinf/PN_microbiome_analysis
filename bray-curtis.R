# This script calculates bray curtis distance between samples and compare between treatment/placebo, change in pnrs level, also performing nmds

library(vegan)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")


## rarefy and calculate bray-curtis distance between samples

### calculate minimal sample size
sample_count = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  summarize(count = sum(count)) %>% 
  ungroup()
min_seqs = sample_count %>% 
  summarize(min = min(count)) %>% 
  pull(min)

### calculate bray-curtis distance between samples
bray_dist = shared %>%
  select(-sample)  %>% 
  avgdist(iterations = 500, sample = min_seqs) 
dist = bray_dist %>% 
  as.matrix() %>% 
  as_tibble(rownames = "sample") %>% 
  pivot_longer(-sample) %>% 
  filter(sample < name) %>% 
  inner_join(., map_sample %>% select(sample, subjid, lesional, visit)) %>% 
  inner_join(., map_sample %>% select(sample, subjid, lesional, visit), by = c("name" = "sample"), suffix = c(".x", ".y")) 


## nmds

set.seed(1)
nmds = metaMDS(bray_dist, k = 3) 
nmds_sample = scores(nmds, display = "sites") %>% 
  as_tibble(rownames = "sample") %>% 
  inner_join(., map_sample) %>% 
  mutate(var = paste0(lesional_visit, trt01p))

### visualizing batch effect
nmds_sample %>% 
  mutate(trt01p = factor(ifelse(trt01p == "Placebo", "Placebo", "Nemolizumab"), levels = c("Placebo", "Nemolizumab"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = factor(batch))) +
  geom_point() +
  stat_ellipse(show.legend = FALSE) +
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), limits = c(-0.6, 0.6)) +
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), limits = c(-0.6, 0.6)) +
  scale_color_manual(name = NULL, breaks = factor(1:3), values = c("red", "green", "blue")) +
  theme_classic() +
  theme(
    legend.position = c(0.9, 0.2),
    legend.text = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, face = "bold", family = "Times"),
    axis.title = element_text(size = 18, face = "bold", family = "Times")
    )
ggsave("manuscript figures/NMDS_batch.tiff", width = 4, height = 4, dpi = 500)

dist_sample = bray_dist %>% 
  as.matrix() %>% 
  as_tibble(rownames = "sample") %>% 
  column_to_rownames("sample")
test = adonis2(
  dist_sample ~ map_sample$batch,
  strata = map_sample$subjid
)
test$`Pr(>F)`
t = betadisper(bray_dist, group = map_sample$batch)
permutest(t)
fisher.test(map_sample$batch, map_sample$community)

### visualizing non-lesional and lesional at baseline
nmds_sample %>% 
  mutate(trt01p = factor(ifelse(trt01p == "Placebo", "Placebo", "Nemolizumab"), levels = c("Placebo", "Nemolizumab"))) %>% 
  filter(visit == "V3") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = lesional)) +
  geom_point() +
  stat_ellipse(show.legend = FALSE) +
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4)) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
  scale_color_manual(
    name = NULL,
    breaks = c("non_lesional", "lesional"),
    values = c("blue", "red"),
    labels = c("non-lesional", "lesional")
    ) +
  theme_classic() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.text = element_text(size = 18, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, face = "bold", family = "Times"),
    axis.title = element_text(size = 18, face = "bold", family = "Times")
    )
ggsave("manuscript figures/NMDS_lesional.nonlesional.tiff", width = 4, height = 4, dpi = 500)

dist_baseline = bray_dist %>% 
  as.matrix() %>% 
  as_tibble(rownames = "sample") %>% 
  pivot_longer(-sample) %>% 
  inner_join(., map_sample %>% select(sample, visit)) %>% 
  inner_join(., map_sample %>% select(sample, visit)
             , by = c("name" = "sample"), suffix = c(".x", ".y")) %>% 
  filter(visit.x == "V3" & visit.y == "V3") %>% 
  select(-visit.x, -visit.y) %>% 
  pivot_wider(names_from = "name", values_from = "value") %>% 
  column_to_rownames("sample")
map_baseline = map_sample %>% filter(visit == "V3")
test = adonis2(
  dist_baseline ~ map_baseline$lesional, 
  strata = map_baseline$subjid
  )
test$test$`Pr(>F)`

t = betadisper(as.dist(dist_baseline), group = map_baseline$lesional)
permutest(t)


