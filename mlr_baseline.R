# This script performs machine learning methods to classify lesional/non-lesional samples at baseline

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ROCR)
library(mikropml)
library(mlr)
library(mlr3)

setwd("/Users/zrayw/Desktop/Alex_Lab/PN_microbiome/analysis")


## generate training/testing data for machine learning

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
  column_to_rownames("sample")
sample_name = rownames(full_data)

### preprocess data using mikropml
preprocessed_data = preprocess_data(full_data, outcome_colname = "lesional")$dat_transformed %>% 
  as.data.frame()
rownames(preprocessed_data) = sample_name
preprocessed_data = preprocessed_data %>% 
  dplyr::rename(c("adtype_IntrinsicAD" = "adtype_Intrinsic AD", 
                  "agegr1_over65years" = "agegr1_> 65 years")) %>% 
  as.data.frame()
baseline_data = preprocessed_data %>% filter(visit_V8 == 0) %>% select(-visit_V8)


## machine learning iterating functions: cv training and prediction on testing

iter_ml_results = function(
    data, learner, ParamSpace, randSearch, cvForTuning, 
    training_fraction, seed){
  size = floor(training_fraction * nrow(data))
  set.seed(seed)
  index = sample(1:nrow(data), size = size)
  training = data[index,]
  testing = data[-index,]
  X_testing = testing %>% select(-lesional)
  y_testing = testing %>% select(lesional)
  Task = makeClassifTask(data = training, target = "lesional")
  set.seed(seed)
  tunedPars = tuneParams(
    learner, task = Task,
    resampling = cvForTuning, par.set = ParamSpace,
    control = randSearch, measures = auc
    )
  tuned = setHyperPars(learner, par.vals = tunedPars$x)
  tunedModel = mlr::train(tuned, Task)
  prediction = predict(tunedModel, newdata = X_testing)$data
  cv_auc = tunedPars$y
  auc = measureAUC(
    probabilities = prediction$prob.lesional, 
    truth = y_testing, positive = "lesional"
    )
  res = data.frame(cv_auc = cv_auc, auc = auc, seed = seed)
  return(res)
}


## train machine learning model on lesional/non-lesional baseline data

### glmnet
learner = makeLearner("classif.glmnet", predict.type = "prob")
ParamSpace = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 0),
  makeNumericParam("lambda", lower = 0.001, upper = 0.5)
)
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      baseline_data, learner, ParamSpace, 
      randSearch, cvForTuning, training_fraction = 0.8,
      seed = x
    )    
    }
)
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_baseline_glm.csv")

### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 100)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      baseline_data, learner, ParamSpace, 
      randSearch, cvForTuning, training_fraction = 0.8,
      seed = x
    )    
  }
)
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_baseline_svm.csv")

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 70)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      baseline_data, learner, ParamSpace, 
      randSearch, cvForTuning, training_fraction = 0.8,
      seed = x
    )    
  }
)
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_baseline_rf.csv")

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 100)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      baseline_data, learner, ParamSpace, 
      randSearch, cvForTuning, training_fraction = 0.8,
      seed = x
    )    
  }
)
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_baseline_xgb.csv")


## compare different model performance

glm_results = read_csv(file = "machine learning data/mlr_baseline_glm.csv") %>% 
  mutate(method = "glm")
svm_results = read.csv(file = "machine learning data/mlr_baseline_svm.csv") %>% 
  mutate(method = "svm")
rf_results = read.csv(file = "machine learning data/mlr_baseline_rf.csv") %>% 
  mutate(method = "rf")
xgb_results = read.csv(file = "machine learning data/mlr_baseline_xgb.csv") %>% 
  mutate(method = "xgb")
combined_results = data.table::rbindlist(
  list(glm_results, svm_results, rf_results, xgb_results)
  )

### wilcox test between different methods
glm = glm_results$auc
svm = svm_results$auc
xgb = xgb_results$auc
rf = rf_results$auc
wilcox.test(svm, glm, paired = TRUE)
wilcox.test(svm, rf, paired = TRUE)
wilcox.test(svm, xgb, paired = TRUE)

combined_results %>% 
  group_by(method) %>% 
  summarize(mean = mean(auc)) %>% View()
combined_results %>% 
  ggplot(aes(x = method, y = auc)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.5, linewidth = 1.2, lty = 2) +
    scale_y_continuous(breaks = 0:5/5, limits = c(0, 1.0)) +
    scale_x_discrete(
      breaks = c("glm", "rf", "svm", "xgb"),
      labels = c("L2-logistic", "RF", "SVM", "XGB")
      ) +
    labs(x = NULL, y = "AUC") +
    theme_classic() +
    theme(
      axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
      axis.text.x = element_text(size = 16, color = "black", face = "bold", family = "Times"),
      axis.text.y = element_text(size = 13, face = "bold", family = "Times")
      )
ggsave("manuscript figures/box.auc_mlr_baseline.tiff", width = 4, height = 3, dpi = 500)
      

## SVM to model baseine lesional/non-lesional

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
tunedPars = tuneParams(
  learner, task = Task,
  resampling = cvForTuning, par.set = ParamSpace,
  control = randSearch, measures = auc
  )
tuned = setHyperPars(learner, par.vals = tunedPars$x)
tunedModel = mlr::train(tuned, Task)
tunedPars$y # avg test auc: 0.7137263 
importance = generateFeatureImportanceData(Task, "permutation.importance", learner)
saveRDS(importance, file = "machine learning data/mlr_baseline_importance.RDS")

### plot top importance
importance = readRDS(file = "machine learning data/mlr_baseline_importance.RDS")
importance$res %>% 
  t() %>% 
  as_tibble(rownames = "features") %>% 
  filter(str_detect(features, "Phylo")) %>% 
  mutate(
    names = include_genus_name
    ) %>% 
  arrange(desc(mmce)) %>% 
  filter(str_detect(features, pattern = "Phylo")) %>% 
  top_n(mmce, n = 10) %>% 
  ggplot(aes(x = mmce, y = reorder(names, mmce))) +
  geom_bar(stat = "identity", fill = "dark grey") +
  scale_x_continuous(position = "top") +
  labs(x = "Feature Importance", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold", family = "Times"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold.italic", family = "Times"),
    axis.text.x = element_text(size = 12, color = "black", family = "Times")
    )
ggsave("manuscript figures/bar.importance_les.nonles.tiff", width = 6, height = 2.5, dpi = 500)
