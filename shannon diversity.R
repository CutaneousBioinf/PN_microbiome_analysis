# This script calculates shannon diversity for all samples after rarefying to the smallest sample size and compare shannon diversity between treatment/placebo, pnrs groups

library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")


## rarefy and calculate Shannon diversity

### calculate minimal sample size
sample_count = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  summarize(count = sum(count)) %>% 
  ungroup()
min_seqs = sample_count %>% 
  summarize(min = min(count)) %>% 
  pull(min) # 28945
shannon_iter = function(){
  shared %>% 
    select(-sample) %>% 
    rrarefy(., min_seqs) %>% 
    diversity(index = "shannon")
}

### calculate Shannon diversity
shannon = replicate(500, shannon_iter()) %>% 
  as_tibble(rownames = "sample", .name_repair = "unique") %>% 
  pivot_longer(-sample) %>% 
  group_by(sample) %>% 
  summarize(shannon = mean(value))
write.csv(shannon, file = "dat/shannon.csv")

### read Shannon diversity
shannon = read.csv(file = "dat/shannon.csv")
shannon_sample =  shannon %>% 
  inner_join(., map_sample)


## compare between non-lesional and lesional

shannon_sample %>%
  filter(visit == "V3") %>%
  select(subjid, lesional, shannon) %>% 
  pivot_wider(names_from = lesional, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-subjid, names_to = "lesional", values_to = "shannon") %>% 
  mutate(lesional = factor(lesional, levels = c("non_lesional", "lesional"))) %>% 
  wilcox.test(shannon ~ lesional, paired = TRUE, data = .) # p = 0.01257
shannon_sample %>%
  filter(visit == "V3") %>%
  select(subjid, lesional, shannon) %>% 
  pivot_wider(names_from = lesional, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-subjid, names_to = "lesional", values_to = "shannon") %>%
  group_by(lesional) %>% 
  summarize(mean = mean(shannon)) # lesional: 2.20 # non-lesional: 2.46
shannon_sample %>%
  filter(visit == "V3") %>%
  select(subjid, lesional, shannon) %>% 
  pivot_wider(names_from = lesional, values_from = shannon) %>%
  drop_na() %>%
  pivot_longer(-subjid, names_to = "lesional", values_to = "shannon") %>%
  mutate(lesional = factor(lesional, levels = c("non_lesional", "lesional"))) %>%
  ggplot(aes(x = lesional, y = shannon)) +
  geom_boxplot() +
  # facet_wrap(sites ~.) +
  labs(title = NULL,
       x = NULL,
       y = "Shannon diversity") +
  scale_x_discrete(breaks = c("non_lesional", "lesional"),
                   labels = c("non-lesional\nskin", "lesional\nskin")) +
  scale_y_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      margin = margin(t = 5, r = 0, b = 0, l = 0), 
      color = "black", size = 16, face = "bold", family = "Times"
      ),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
    axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size = 16, face = "bold", family = "Times")
    ) 
ggsave("manuscript figures/box.shannon_lesional.nonlesional.tiff", width = 4, height = 3.5, dpi = 500)

shannon_sample %>%
  filter(SQLA == SQNLA) %>% 
  filter(visit == "V3") %>%
  select(subjid, lesional, shannon, sites) %>% 
  pivot_wider(names_from = lesional, values_from = shannon) %>% 
  # drop_na() %>% 
  # pivot_longer(-c(subjid, sites), names_to = "lesional", values_to = "shannon") %>%
  ggplot(aes(x = non_lesional, y = lesional)) +
  geom_point() +
  facet_wrap(sites ~.) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(
    title = NULL,
    x = "Shannon diversity in non-lesional skin",
    y = "Shannon diversity in lesional skin"
    ) +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", family = "Times"),
    axis.title = element_text(size = 13, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, face = "bold", family = "Times")
  )
ggsave("manuscript figures/box.shannon_lesional.nonlesional_same.site.tiff", width = 7, height = 7, dpi = 500)


## compare Shannon diversity between responder & non-responders

shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3"))
    ) %>% 
  group_by(comparison) %>% 
  summarize(mean = mean(shannon)) # Placebo: 2.398178, treatment non-responders: 2.770769, treatment responders: 1.930929
shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3"))
  ) %>%
  filter(comparison != "1") %>% 
  wilcox.test(shannon ~ comparison, data = .) # treatment non-responders vs responders: p = 0.02302
shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3"))
  ) %>%
  filter(comparison != "2") %>% 
  wilcox.test(shannon ~ comparison, data = .) # placebo vs responders: p = 0.07145
shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3")),
    comparison = 3*as.numeric(comparison)
    ) %>% 
  ggplot(aes(x = comparison, y = shannon, group = comparison)) +
  geom_boxplot() +
  labs(title = NULL,
       x = NULL,
       y = "Shannon diversity") +
  scale_x_continuous(
    breaks = 1:3*3, 
    labels = c("placebos", "non-\nresponded\ntreatments", "responded\ntreatments")
    ) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      margin = margin(t = 5, r = 0, b = 0, l = 0), color = "black", 
      size = 16, face = "bold", family = "Times"),
      axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
      axis.title.y = element_text(
        margin = margin(t = 0, r = 8, b = 0, l = 0), size = 16, face = "bold", family = "Times")
    )
ggsave("manuscript figures/box.shannon_placebo.responded.nonresponded_lesional_week12.tiff", width = 4.5, height = 4.5, dpi = 500)
shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3")),
    comparison = 3*as.numeric(comparison)
    ) %>% 
  ggplot(aes(x = comparison, y = shannon, group = comparison)) +
  geom_boxplot() +
  facet_wrap(sites ~.) +
  labs(title = NULL,x = NULL, y = "Shannon diversity") +
  scale_x_continuous(breaks = 1:3*3, labels = c("placebos", "non-responded\ntreatments", "responded\ntreatments")) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", family = "Times"),
    axis.text.x = element_text(
      color = "black", size = 13, face = "bold", family = "Times", angle = -90, hjust = 0.5, vjust = 0.5
      ),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
    axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size = 16, face = "bold", family = "Times")
    )

model.diversity.3 = shannon_sample %>% 
  filter(visit == "V8") %>% 
  mutate(
    lesional = factor(lesional, levels = c("non_lesional", "lesional")),
    comparison = ifelse(trt01p == "Placebo", "1", ifelse(nrs_w12 == "N", "2", "3"))
    ) %>%
  lm(shannon ~ agegr1 + asex + adtype + acountry + comparison, data = .) 
car::linearHypothesis(model.diversity.3, hypothesis.matrix = c(rep(0, 7), 1, -1), rhs = 0) # non-responded vs responded: p = 0.006361
summary(model.diversity.3) # non-responded vs responded: beta = -0.96145


## scatter plot: shannon v.s. pnrs

shannon_sample %>% 
  filter(lesional == "lesional") %>% 
  mutate(trt01p = factor(ifelse(trt01p == "Placebo", "Placebo", "Nemolizumab"), levels = c("Placebo", "Nemolizumab")),
         visit = ifelse(visit == "V3", "Baseline", "Week 12")) %>% 
  ggplot(aes(x = pnrs, y = shannon)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, linetype = 1) +
  scale_x_continuous(
    breaks = c(0, 2, 4 , 6, 8, 10),
    limits = c(-0.5, 10.5),
    expand = c(0, 0)
    ) +
  scale_y_continuous(
    limits = c(0, 4), 
    expand = c(0, 0)) +
  facet_grid(visit ~ trt01p) +
  labs(title = NULL, 
       x = "PP NRS level", 
       y = "Shannon diversity") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, family = "Times"),
    strip.text = element_text(size = 16, face = "bold", family = "Times"),
    plot.title = element_text(hjust = 0.5),
    panel.grid=element_blank(),
    strip.background = element_blank(),
    axis.line=element_line()
    )
ggsave("manuscript figures/scatter.ppnrs.shannon_treatment_visit.tiff", width = 4, height = 4, dpi = 500)

model.diversity.1 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + pnrs, data = .)
summary(model.diversity.1) # lesional week12: ppnrs: beta = 0.07422, p = 0.0368
model.diversity.1.1 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + pnrs, data = .)
summary(model.diversity.1.1) # lesional week 12, treatment: ppnrs: beta = 0.19780, p =  0.00176
model.diversity.1.2 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8" & trt01p == "Placebo") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + pnrs, data = .)
summary(model.diversity.1.2) # lesional week 12, placebo: ppnrs: beta = -0.008296, p =  0.895754
model.diversity.1.3 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V3" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + pnrs, data = .)
summary(model.diversity.1.3) # lesional baseline, treatment: ppnrs: beta = 0.16225, p =  0.1616   
model.diversity.1.4 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V3" & trt01p == "Placebo") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + pnrs, data = .)
summary(model.diversity.1.4) # lesional baseline, placebo: ppnrs: beta = -0.2434, p =  0.242  


## compare between baseline and week 12 in lesional

shannon_sample %>% 
  filter(lesional == "lesional") %>% 
  select(subjid, visit, shannon) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  mutate(diff = V8 - V3) %>% 
  select(subjid, diff) %>% 
  inner_join(., map) %>% 
  mutate(
    trt = ifelse(trt01p == "Placebo", "Placebo", "Nemolizumab"), 
    trt = factor(trt, levels = c("Placebo", "Nemolizumab"))
    ) %>% 
  ggplot(aes(x = -chg, y = diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = NULL,
       x = "Decrease of PP NRS level\nfrom baseline to week 12",
       y = "Change of Shannon diversity\nfrom baseline to week 12") +
  facet_wrap(~trt, scales = "free_x") +
  theme_bw() +
  theme(
    plot.margin = margin(l = 15, r = 15, t = 30),
    axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
    axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, family = "Times"),
    strip.text = element_text(size = 16, face = "bold", family = "Times"),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line()
    )
ggsave("manuscript figures/scatter.decrease.ppnrs.diffshannon_treatment.tiff", width = 4.2, height = 3)

model.diversity.2.1 = shannon_sample %>% 
  filter(lesional == "lesional") %>% 
  select(subjid, visit, shannon) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  mutate(diff = V8 - V3) %>% 
  select(subjid, diff) %>% 
  inner_join(., map) %>% 
  mutate(trt01p = factor(trt01p, levels = c("Placebo", "Nemolizumab 0.5mg/kg"))) %>% 
  filter(trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lm(diff ~ agegr1 + asex + adtype + acountry + I(-chg), data = .) # treatment: beta = -0.24045, p = 0.0174
summary(model.diversity.2.1)
model.diversity.2.2 = shannon_sample %>% 
  filter(lesional == "lesional") %>% 
  select(subjid, visit, shannon) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  mutate(diff = V8 - V3) %>% 
  select(subjid, diff) %>% 
  inner_join(., map) %>% 
  mutate(trt01p = factor(trt01p, levels = c("Placebo", "Nemolizumab 0.5mg/kg"))) %>% 
  filter(trt01p == "Placebo") %>% 
  lm(diff ~ agegr1 + asex + adtype + acountry + I(-chg), data = .) # treatment: beta = -0.1078, p = 0.3063 
summary(model.diversity.2.2)

re_abund = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  select(-count) %>% 
  pivot_wider(names_from = "otu", values_from = "re_abund")
re_abund_sample = inner_join(re_abund, shannon_sample, by = "sample")
cor.test(re_abund_sample$Phylo0001, re_abund_sample$shannon, method = "spearman")
re_abund_sample %>% 
  filter(lesional == "lesional") %>% 
  cor.test(formula = ~ `Phylo0001` + shannon, method = "spearman", data = .)
re_abund_sample %>% 
  filter(lesional == "non_lesional") %>% 
  cor.test(formula = ~ `Phylo0001` + shannon, method = "spearman", data = .)


## Association with IGA

iga = read.table(file = "dat/PN_IGA", sep = "\t", header = TRUE) %>% 
  rename_all(tolower) %>% 
  rename(
    "iga" = "aval", "iga_c" = "avalc", "iga_base" = "base", "iga_chg" = "chg"
    ) %>% 
  mutate(
    subjid = str_replace_all(subjid, pattern = "-", replacement = ""),
    iga_chg = replace_na(iga_chg, 0)
    ) %>% 
  filter(avisit != "")   
shannon_sample = iga %>% 
  filter(avisit == "Baseline" | avisit == "Week 12") %>% 
  mutate(visit = ifelse(avisit == "Baseline", "V3", "V8")) %>% 
  select(subjid, visit, iga, iga_c, iga_base, iga_chg) %>% 
  inner_join(shannon_sample, ., by = c("subjid", "visit")) 
shannon_sample %>% 
  filter(lesional == "lesional") %>% 
  mutate(
    trt01p = factor(ifelse(trt01p == "Placebo", "Placebo", "Nemolizumab"), levels = c("Placebo", "Nemolizumab")),
    visit = ifelse(visit == "V3", "Baseline", "Week 12")
    ) %>% 
  ggplot(aes(x = factor(iga), y = shannon)) +
  geom_boxplot() +
  geom_smooth(method = lm, se = FALSE, linetype = 1) +
  scale_y_continuous(limits = c(0, 4), 
                     expand = c(0, 0)) +
  facet_grid(visit ~ trt01p) +
  labs(title = NULL, 
       x = "IGA level", 
       y = "Shannon diversity") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
    axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
    axis.text = element_text(size = 12, family = "Times"),
    strip.text = element_text(size = 16, face = "bold", family = "Times"),
    plot.title = element_text(hjust = 0.5),
    panel.grid=element_blank(),
    strip.background = element_blank(),
    axis.line=element_line()
    )
ggsave("manuscript figures/scatter.shannon.iga_treatment_visit.tiff", width = 4.2, height = 3, dpi = 500)

mod_1.1 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.1)
mod_1.2 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8" & trt01p == "Placebo") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.2)
mod_1.3 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V3" & trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.3)
mod_1.4 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V3" & trt01p == "Placebo") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.4)
mod_1.5 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V8") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.5)
mod_1.6 = shannon_sample %>% 
  filter(lesional_visit == "lesional_V3") %>% 
  lm(shannon ~ agegr1 + asex + adtype + acountry + factor(iga), data = .)
summary(mod_1.6)










## Test between non-lesional and lesional, set subjid and batch as random effect

m_1 = shannon_sample %>%
  filter(visit == "V3") %>%
  select(subjid, lesional, shannon, visit) %>% 
  pivot_wider(names_from = lesional, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-c("subjid", "visit"), names_to = "lesional", values_to = "shannon") %>% 
  inner_join(., map_sample, by = c("subjid", "lesional", "visit")) %>% 
  lmer(shannon ~ agegr1 + asex + adtype + acountry + lesional + (1|subjid) + (1|batch), data = .)
anova(m_1)
summary(m_1)


## Test between baseline and week 12, set subjid and batch as random effect
m_2 = shannon_sample %>%
  filter(lesional == "lesional") %>%
  select(subjid, lesional, shannon, visit) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-c("subjid", "lesional"), names_to = "visit", values_to = "shannon") %>% 
  inner_join(., map_sample, by = c("subjid", "lesional", "visit")) %>% 
  lmer(shannon ~ agegr1 + asex + adtype + acountry + visit*trt01p + (1|subjid) + (1|batch), data = .)
anova(m_2)
summary(m_2)

m_2.1 = shannon_sample %>%
  filter(lesional_visit == "lesional_V8" & trt01p == "Nemolizumab 0.5mg/kg") %>%
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + (1|batch), data = .)
anova(m_2.1)
summary(m_2.1)

m_2.2 = shannon_sample %>%
  filter(lesional_visit == "lesional_V8" & trt01p == "Placebo") %>%
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + (1|batch), data = .)
anova(m_2.2)
summary(m_2.2)

m_2.3 = shannon_sample %>%
  filter(lesional_visit == "lesional_V3" & trt01p == "Nemolizumab 0.5mg/kg") %>%
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + (1|batch), data = .)
anova(m_2.3)
summary(m_2.3)

m_2.4 = shannon_sample %>%
  filter(lesional_visit == "lesional_V3" & trt01p == "Placebo") %>%
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + (1|batch), data = .)
anova(m_2.4)
summary(m_2.4)

m_3 = shannon_sample %>%
  filter(lesional == "lesional") %>%
  select(subjid, lesional, shannon, visit) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-c("subjid", "lesional"), names_to = "visit", values_to = "shannon") %>% 
  inner_join(., map_sample, by = c("subjid", "lesional", "visit")) %>% 
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + visit + (1|subjid) + (1|batch), data = .)
anova(m_3)
summary(m_3)

m_3.1 = shannon_sample %>%
  filter(lesional == "lesional") %>%
  select(subjid, lesional, shannon, visit) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-c("subjid", "lesional"), names_to = "visit", values_to = "shannon") %>% 
  inner_join(., map_sample, by = c("subjid", "lesional", "visit")) %>% 
  filter(trt01p == "Nemolizumab 0.5mg/kg") %>% 
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + visit + (1|subjid) + (1|batch), data = .)
anova(m_3.1)
summary(m_3.1)

m_3.2 = shannon_sample %>%
  filter(lesional == "lesional") %>%
  select(subjid, lesional, shannon, visit) %>% 
  pivot_wider(names_from = visit, values_from = shannon) %>% 
  drop_na() %>% 
  pivot_longer(-c("subjid", "lesional"), names_to = "visit", values_to = "shannon") %>% 
  inner_join(., map_sample, by = c("subjid", "lesional", "visit")) %>% 
  filter(trt01p == "Placebo") %>% 
  lmer(shannon ~ agegr1 + asex + adtype + acountry + pnrs + visit + (1|subjid) + (1|batch), data = .)
anova(m_3.2)
summary(m_3.2)
