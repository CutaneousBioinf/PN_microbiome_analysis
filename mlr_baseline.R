# This script performs machine learning methods to classify lesional/non-lesional samples at baseline

library(mikropml)
library(purrr)
library(broom)
library(ROCR)
library(mlr)
library(mlr3)
library(tidyverse)

include_genus = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  group_by(otu) %>% 
  summarize(mean = mean(re_abund)) %>% 
  filter(mean >= 0.001) %>% 
  pull(otu) 
include_genus_name = taxonomy %>% filter(otu %in% include_genus) %>% pull(genus)

demographics = map_sample %>% 
  select(sample, agegr1, asex, adtype, acountry, lesional, visit) %>% 
  mutate(lesional = as.character(lesional))
include_genus_re_abund = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  filter(otu %in% include_genus) %>%
  pivot_wider(names_from = "otu", values_from = "re_abund")

full_data = inner_join(demographics, include_genus_re_abund, by = "sample") %>% 
  as.data.frame()
sample = full_data$sample
full_data = full_data %>% select(-sample)

preprocessed_data = preprocess_data(full_data, outcome_colname = "lesional")$dat_transformed %>% 
  as.data.frame()
rownames(preprocessed_data) = sample
preprocessed_data = preprocessed_data %>% 
  dplyr::rename(c("adtype_IntrinsicAD" = "adtype_Intrinsic AD", 
                  "agegr1_over65years" = "agegr1_> 65 years")) %>% 
  as.data.frame()

baseline_data = preprocessed_data %>% filter(visit_V8 == 0) %>% select(-visit_V8)
week12_data = preprocessed_data %>% filter(visit_V8 == 1) %>% select(-visit_V8)

X_test = week12_data %>% select(-lesional)

## functions
iter_ml_results = function(data, 
                           learner, 
                           ParamSpace, 
                           randSearch, 
                           cvForTuning, 
                           training_fraction, 
                           seed){
  size = floor(training_fraction * nrow(data))
  set.seed(seed)
  index = sample(1:nrow(data), size = size)
  training = data[index,]
  testing = data[-index,]
  X_testing = testing %>% select(-lesional)
  y_testing = testing %>% select(lesional)
  Task = makeClassifTask(data = training, target = "lesional")
  set.seed(seed)
  tunedPars = tuneParams(learner, 
                         task = Task,
                         resampling = cvForTuning,
                         par.set = ParamSpace,
                         control = randSearch, 
                         measures = auc)
  tuned = setHyperPars(learner, par.vals = tunedPars$x)
  tunedModel = mlr::train(tuned, Task)
  prediction = predict(tunedModel, newdata = X_testing)$data
  cv_auc = tunedPars$y
  auc = measureAUC(probabilities = prediction$prob.lesional, 
                   truth = y_testing, 
                   positive = "lesional")
  res = data.frame(cv_auc = cv_auc, auc = auc, seed = seed)
  return(res)
}

### glmnet
learner = makeLearner("classif.glmnet", predict.type = "prob")
ParamSpace = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 0),
  makeNumericParam("lambda", lower = 0.001, upper = 0.5)
)
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

glm_results = lapply(1:50, function(x){
  iter_ml_results(baseline_data, 
                  learner, 
                  ParamSpace, 
                  randSearch,
                  cvForTuning,
                  training_fraction = 0.8,
                  seed = x)}
  )
combine_glm = glm_results %>% data.table::rbindlist(.)
write_csv(combine_glm, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_baseline_glm.csv")

### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 100)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

svm_results = lapply(1:50, function(x){
  iter_ml_results(baseline_data, 
                  learner, 
                  ParamSpace, 
                  randSearch,
                  cvForTuning,
                  training_fraction = 0.8,
                  seed = x)}
  )

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 70)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

rf_results = lapply(1:50, function(x){
  iter_ml_results(baseline_data, 
                  learner, 
                  ParamSpace, 
                  randSearch,
                  cvForTuning,
                  training_fraction = 0.8,
                  seed = x)}
)

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 100)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

xgb_results = lapply(1:50, function(x){
  iter_ml_results(baseline_data, 
                  learner, 
                  ParamSpace, 
                  randSearch,
                  cvForTuning,
                  training_fraction = 0.8,
                  seed = x)}
)


## compare different model performance

wd = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/"
setwd(wd)
glm_results = read_csv(file = "mlr_baseline_glm.csv") %>% 
  mutate(method = "glm")
svm_results = read.csv(file = "mlr_baseline_svm.csv") %>% 
  mutate(method = "svm")
rf_results = read.csv(file = "mlr_baseline_rf.csv") %>% 
  mutate(method = "rf")
xgb_results = read.csv(file = "mlr_baseline_xgb.csv") %>% 
  mutate(method = "xgb")
combine_results = data.table::rbindlist(list(glm_results, 
                                             svm_results, 
                                             rf_results, 
                                             xgb_results))

combine_results %>% 
  group_by(method) %>% 
  summarize(mean = mean(auc)) %>% View()
combine_results %>% 
  ggplot(aes(x = method, y = auc)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.5, linewidth = 1.2, lty = 2) +
    scale_y_continuous(breaks = 1:5/5, 
                       limits = c(0.2, 1.0)) +
    scale_x_discrete(breaks = c("glm", "rf", "svm", "xgb"),
                     labels = c("L2-logistic", "RF", "SVM", "XGB")) +
    labs(x = NULL, y = "Mean AUC") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
          axis.text.x = element_text(size = 16, 
                                     color = "black", face = "bold", family = "Times"),
          axis.text.y = element_text(size = 13, face = "bold", family = "Times"))
ggsave("/Users/zrayw/Desktop/PN_microbiome/analysis/manuscript figures/mlr_baseline1.tiff", width = 4, height = 3, dpi = 500)
      

## SVM to model baseine
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

Task = makeClassifTask(data = baseline_data, target = "lesional")
set.seed(10086)
tunedPars = tuneParams(learner, 
                       task = Task,
                       resampling = cvForTuning,
                       par.set = ParamSpace,
                       control = randSearch, 
                       measures = auc)
tuned = setHyperPars(learner, par.vals = tunedPars$x)
tunedModel = mlr::train(tuned, Task)
tunedPars$y

importance = generateFeatureImportanceData(Task, "permutation.importance",
                                           learner)
saveRDS(importance, file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_baseline_importance.Rds")
importance = readRDS(file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_baseline_importance.Rds")

importance$res %>% 
  t() %>% 
  as_tibble(rownames = "features") %>%
  mutate(names = c(include_genus_name, 
                   "Austria", "France", "Germany", "Poland", 
                   "Intrinsic AD", "age > 65", "Male")) %>% 
  arrange(desc(mmce)) %>% 
  filter(str_detect(features, pattern = "Phylo")) %>% 
  top_n(mmce, n = 10) %>% View()
  ggplot(aes(x = mmce, y = reorder(names, mmce))) +
  geom_bar(stat = "identity", fill = "dark grey") +
  scale_x_continuous(position = "top") +
  labs(x = "Feature Importance",
       y = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13.5, color = "black", face = "bold.italic", family = "Times"))
ggsave("/Users/zrayw/Desktop/PN_microbiome/analysis/manuscript figures/feature importance2.tiff", width = 6, height = 2.5, dpi = 500)


## predict week 12 status

predict_week12 = predict(tunedModel, newdata = X_test)
predict_week12$data %>% 
  as_tibble(rownames = "sample") %>% 
  inner_join(., map_sample) %>% 
  group_by(group1) %>% 
  summarize(sum = sum(prob.lesional <= 0.5),
            n = n())


## Iterative feature importance from svm
iter_imp_results = function(data, 
                            learner, 
                            training_fraction, 
                            seed){
  size = floor(training_fraction * nrow(data))
  set.seed(seed)
  index = sample(1:nrow(data), size = size)
  training = data[index,]
  testing = data[-index,]
  X_testing = testing %>% select(-lesional)
  y_testing = testing %>% select(lesional)
  Task = makeClassifTask(data = training, target = "lesional")
  set.seed(seed)
  importance = generateFeatureImportanceData(Task, "permutation.importance",
                                             learner)
  return(importance$res)
}

learner = makeLearner("classif.svm", predict.type = "prob")
svm_results = lapply(1:50, function(x){
  iter_imp_results(baseline_data, 
                   learner, 
                   training_fraction = 0.8,
                   seed = x)})

results = svm_results %>% 
  data.table::rbindlist() %>% 
  t() %>% 
  as_tibble(rownames = "taxon") %>% 
  pivot_longer(-taxon, names_to = "1", values_to = "value") 

svm_results %>% 
  data.table::rbindlist() %>% 
  t() %>% 
  as_tibble(rownames = "taxon") %>% 
  pivot_longer(-taxon, names_to = "1", values_to = "value") %>% 
  group_by(taxon) %>% 
  summarize(mean = mean(value)) %>% 
  arrange(desc(mean)) %>% 
  View()



