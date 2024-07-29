## This script performs genus descriptive analysis (Figure 2)

library(ggpattern)

## common genus

common_genus_taxonomy = taxonomy %>% filter(otu %in% common_genus) %>% arrange(desc(size))
common_genus_taxonomy %>% 
  ggplot(aes(x = 1:length(common_genus), y = size)) +
  geom_point(aes(color = phylum), size = 3.5) +
  geom_line() +
  scale_x_continuous(breaks = 1:length(common_genus), labels = common_genus_taxonomy$genus) +
  scale_y_continuous(
    breaks = c(300000, 1:7*1000000),
    labels = c("0.3", "1", "2", "3", "4", "5", "6", "7"),
    limits = c(0, 7100000)
    ) +
  scale_color_manual(name = NULL, values = c("red", "blue", "black")) +
  labs(x = NULL, y = "Counts per million") +
  theme_classic() +
  theme(
    plot.margin = margin(l = 30),
    legend.position = c(0.75, 0.8),
    legend.text = element_text(size = 18, face = "bold.italic", family = "Times"),
    axis.text.x = element_text(
      color = "black", size = 18, face = "bold.italic", family = "Times",
      angle = 45, hjust = 1, vjust = 1
      ), 
    axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 0),
      size = 21, face = "bold", family = "Times")
    )
ggsave("manuscript figures/common_genus.tiff", width = 4.8, height = 3.9, dpi = 500)


## pie chart phylum

shared_phylum = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., taxonomy %>% select(otu, phylum)) %>% 
  group_by(sample, phylum) %>% 
  summarize(count = sum(count)) %>% 
  pivot_wider(names_from = "phylum", values_from = "count") %>% 
  as.data.frame()
rownames(shared_phylum) = shared_phylum$sample
common_phylum = shared_phylum %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund), .groups = "drop") %>% 
  filter(mean >= 0.01) %>% 
  pull(otu)
common_phylum_count = shared_phylum %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(otu) %>% 
  summarize(count = sum(count), .groups = "drop") %>% 
  mutate(total = sum(count)) %>% 
  filter(otu %in% common_phylum) %>% 
  mutate(`Rare phylum` = total - sum(count)) %>% 
  pivot_wider(names_from = otu, values_from = count) %>% 
  pivot_longer(-total, names_to = "otu", values_to = "count") %>% 
  mutate(re_abund = paste(round(count/total*100, 2), "%"))
common_phylum_count %>% 
  mutate(otu = factor(otu, levels = c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Rare phylum"))) %>% 
  ggplot(aes(x = "", y = count, fill = otu)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) + 
  scale_fill_brewer(name = NULL, palette = "Blues", direction = -1) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.position = "none",
    legend.text = element_text(size = 16, face = "bold.italic", family = "Times")
    )
ggsave("manuscript figures/pie_common.phylum.tiff", width = 4, height = 4)


## Firmicutes

Firmicutes_common_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., taxonomy, by = "otu") %>% 
  filter(phylum == "Firmicutes") %>% 
  mutate(total = sum(count)) %>% 
  filter(otu %in% common_genus) %>% 
  group_by(otu, total) %>% 
  summarize(count = sum(count), .groups = "drop") %>% 
  mutate(Others = total - sum(count)) %>% 
  inner_join(., taxonomy %>% select(otu, genus)) %>% 
  select(-otu) %>% 
  pivot_wider(names_from = genus, values_from = count) %>% 
  pivot_longer(-total, names_to = "genus", values_to = "count") %>% 
  mutate(re_abund = round(count/total*100, digits = 2))
Firmicutes_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Staphylococcus", "Others", "Streptococcus", "Anaerococcus", "Finegoldia"))) %>% 
  ggplot(aes(x = "", y = count, fill = genus, width = 0.35)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), aes(label = re_abund)) +
  scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
  theme_void() +
  theme(
    legend.position = c(0.15, 0.5),
    legend.text = element_text(size = 12)
    )
Firmicutes_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Staphylococcus", "Others", "Streptococcus", "Anaerococcus", "Finegoldia"))) %>% 
  ggplot(aes(x = genus, y = re_abund, fill = genus)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.95) +
  scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 82)) +
  labs(title = NULL, x = NULL, y = "Relative Abundance (%)" ) +
  theme_classic() +
  theme(
    plot.margin = margin(t = 30, r = 5, b = 10, l = 15),
    legend.position = "none",
    axis.text.x = element_text(
      color = "black", angle = 45, hjust = 1, vjust = 1,
      size = 18, face = "bold.italic", family = "Times"
      ),
    axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 8, b = 0, l = 0),
      size = 20, face = "bold", family = "Times")
    )
ggsave("manuscript figures/bar.re.abund_Firmicutes_common.genus.tiff", width = 4.2, height = 4.2)


## Actinobacteria

Actinobacteria_common_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., taxonomy, by = "otu") %>% 
  filter(phylum == "Actinobacteria") %>% 
  mutate(total = sum(count)) %>% 
  filter(otu %in% common_genus) %>% 
  group_by(otu, total) %>% 
  summarize(count = sum(count), .groups = "drop") %>% 
  mutate(Others = total - sum(count)) %>% 
  inner_join(., taxonomy %>% select(otu, genus)) %>% 
  select(-otu) %>% 
  pivot_wider(names_from = genus, values_from = count) %>% 
  pivot_longer(-total, names_to = "genus", values_to = "count") %>% 
  mutate(re_abund = round(count/total*100, digits = 2))
Actinobacteria_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Cutibacterium", "Corynebacterium", "Others", "Micrococcus", "Kocuria"))) %>% 
  ggplot(aes(x = "", y = count, fill = genus, width = 0.35)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
  theme_void() +
  theme(
    legend.position = c(0.88, 0.5),
    legend.text = element_text(size = 12)
    )
Actinobacteria_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Cutibacterium", "Corynebacterium", "Others", "Micrococcus", "Kocuria"))) %>% 
  ggplot(aes(x = genus, y = re_abund, fill = genus)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.95) +
    scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
    scale_y_continuous(position = "right", expand = c(0, 0), limits = c(0, 82)) +
    labs(title = NULL, x = NULL, y = "Relative Abundance (%)" ) +
    theme_classic() +
  theme(
    plot.margin = margin(t = 30, r = 10, b = 10, l = 50),
    legend.position = "none",
    axis.text.x = element_text(
      color = "black", angle = 45, hjust = 1, vjust = 1,
      size = 18, face = "bold.italic", family = "Times"
      ),
      axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
      axis.title.y = element_text(
        margin = margin(t = 0, r = 10, b = 0, l = 0),
        size = 22, face = "bold", family = "Times")
    )
ggsave("manuscript figures/bar.re.abund_Actinobacteria_common.genus.tiff", width = 4.5, height = 5, dpi = 500)


## Proteobacteria

Proteobacteria_common_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., taxonomy, by = "otu") %>% 
  filter(phylum == "Proteobacteria") %>% 
  mutate(total = sum(count)) %>% 
  filter(otu %in% common_genus) %>% 
  group_by(otu, total) %>% 
  summarize(count = sum(count), .groups = "drop") %>% 
  mutate(Others = total - sum(count)) %>% 
  inner_join(., taxonomy %>% select(otu, genus)) %>% 
  select(-otu) %>% 
  pivot_wider(names_from = genus, values_from = count) %>% 
  pivot_longer(-total, names_to = "genus", values_to = "count") %>% 
  mutate(re_abund = round(count/total*100, digits = 2))
Proteobacteria_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Others", "Paracoccus", "Acinetobacter", "Enhydrobacter"))) %>% 
  ggplot(aes(x = "", y = count, fill = genus, width = 0.35)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
  theme_void() +
  theme(
    legend.position = c(0.88, 0.5),
    legend.text = element_text(size = 12)
    )
Proteobacteria_common_genus %>% 
  mutate(genus = factor(genus, levels = c("Others", "Paracoccus", "Acinetobacter", "Enhydrobacter"))) %>% 
  ggplot(aes(x = genus, y = re_abund, fill = genus)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.95) +
  scale_fill_brewer(name = NULL, palette = "Reds", direction = -1) +
  scale_y_continuous(position = "right",
                     expand = c(0, 0),
                     limits = c(0, 80)) +
  labs(title = NULL, x = NULL, y = "Relative Abundance (%)" ) +
  theme_classic() +
  theme(
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10),
    legend.position = "none",
    axis.text.x = element_text(
      color = "black", angle = 45, hjust = 1, vjust = 1,
      size = 18, face = "bold.italic", family = "Times"
      ),
    axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 15, b = 0, l = 0),
      size = 22, face = "bold", family = "Times")
    )
ggsave("manuscript figures/bar.re.abund_Proteobacteria_common.genus.tiff", width = 4, height = 4.5)


## comparing between lesional non-lesional abundant genus

de_genus_lesional = read_csv(file = "ANCOMBC results/de_genus_lesional.csv")
indicator_de_genus_lesional = de_genus_lesional %>% 
  mutate(indicator = (p.adj <= 0.1)*1) %>% 
  select(genus, indicator)

common_genus_lesionalvisit = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  inner_join(., map_sample %>% select(sample, lesional_visit)) %>% 
  group_by(otu, lesional_visit) %>% 
  summarize(count = sum(count), .groups = "drop") %>% 
  group_by(lesional_visit) %>% 
  mutate(total = sum(count)) %>% 
  filter(otu %in% common_genus) %>% 
  mutate(`Rare` = total - sum(count)) %>% 
  inner_join(., taxonomy %>% select(otu, genus)) %>% 
  select(-otu) %>% 
  pivot_wider(names_from = genus, values_from = count) %>% 
  pivot_longer(-c(total, lesional_visit), names_to = "genus", values_to = "count") %>% 
  group_by(lesional_visit) %>% 
  mutate(re_abund = count/sum(count))

common_genus_lesionalvisit %>%
  left_join(., indicator_de_genus_lesional, by = "genus") %>% 
  replace(is.na(.), 0) %>% 
  mutate(
    genus = factor(
      genus,
      levels = c(
        "Rare", "Paracoccus", "Kocuria", "Finegoldia", "Acinetobacter", "Enhydrobacter",
        "Anaerococcus", "Micrococcus", "Streptococcus","Corynebacterium", "Cutibacterium", "Staphylococcus"
        )
      ), 
      lesional_visit = factor(lesional_visit, levels = c("non_lesional_V3", "lesional_V3", "lesional_V8")),
      lesional_visit = case_when(
        lesional_visit == "non_lesional_V3" ~ 2.5,
        lesional_visit == "lesional_V3" ~ 5,
        lesional_visit == "lesional_V8" ~ 7.5
        ),
      indicator = factor(indicator)
    ) %>% 
  ggplot(aes(x = lesional_visit, y = re_abund, fill = genus, pattern = indicator)) +
  geom_bar_pattern(
    stat = "identity", position = "stack", width = 1.9, color = "black", 
    pattern_fill = "black", pattern_angle = 45, pattern_density = 0.1,
    pattern_spacing = 0.025, pattern_key_scale_factor = 0.6
    ) +
  scale_y_continuous(
    expand = c(0, 0), breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = c(0, 25, 50, 75, 100)
    ) + 
  scale_x_continuous(
    breaks = c(2.5, 5, 7.5), labels = c("baseline\nnon-lesional", "baseline\nlesional", "week 12\nlesional")
    ) +
  scale_fill_manual(name = NULL, values = rainbow(12)) +
  scale_pattern_manual(breaks = c("0", "1"), values = c("none", "stripe"), guide = "none") +
  labs(title = NULL, x = NULL, y = "Relative Abundance (%)") +
  guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 20, face = "bold.italic", family = "Times"),
    axis.text.x = element_text(color = "black", size = 18, face = "bold", family = "Times"),
    axis.text.y = element_text(size = 15, face = "bold", family = "Times"),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 5, b = 0, l = 0),
      size = 21, face = "bold", family = "Times")
    )
ggsave("manuscript figures/stack.bar.re.abund_common.genus_lesional_non.lesional.tiff", width = 7, height = 4, dpi = 500)

common_genus_lesionalvisit %>% 
  select(lesional_visit, genus, re_abund) %>% 
  pivot_wider(names_from = "genus", values_from = "re_abund") %>% 
  View()


## genus correlation
re_abund = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  select(-count) %>% 
  pivot_wider(names_from = "otu", values_from = "re_abund")
re_abund_sample = inner_join(re_abund, map_sample, by = "sample")

cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0003, method = "pearson")
cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0009, method = "pearson")
cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0028, method = "pearson")

cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0015, method = "pearson")
cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0006, method = "pearson")
cor.test(re_abund_sample$Phylo0001, re_abund_sample$Phylo0007, method = "pearson")


