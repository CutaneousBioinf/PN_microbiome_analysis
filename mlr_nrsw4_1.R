# This script performs machine learning methods on baseline lesional samples in treatment groups to predict week 4 responding status

library(mikropml)
library(purrr)
library(broom)
library(ROCR)
library(mlr)
library(mlr3)
library(tidyverse)


demographics = map_sample %>% 
  filter(trt01p == "Nemolizumab 0.5mg/kg" & lesional_visit == "lesional_V3") %>%
  select(sample, agegr1, asex, adtype, acountry, pnrs, nrs_w4)

genus_re_abund = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  filter(otu %in% abund_genus) %>% 
  pivot_wider(names_from = "otu", values_from = "re_abund")

full_data = inner_join(demographics, genus_re_abund, by = "sample") %>% select(-sample)
preprocessed_data = preprocess_data(full_data, outcome_colname = "nrs_w4")$dat_transformed


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
  X_testing = testing %>% select(-nrs_w4)
  y_testing = testing %>% select(nrs_w4)
  Task = makeClassifTask(data = training, target = "nrs_w4")
  set.seed(seed)
  tunedPars = tuneParams(learner, 
                         task = Task,
                         resampling = cvForTuning,
                         par.set = ParamSpace,
                         control = randSearch, 
                         measures = acc)
  tuned = setHyperPars(learner, par.vals = tunedPars$x)
  tunedModel = mlr::train(tuned, Task)
  prediction = predict(tunedModel, newdata = X_testing)$data
  cv_acc = tunedPars$y
  auc = measureAUC(probabilities = prediction$prob.Y, 
                   truth = y_testing, 
                   positive = "Y")
  res = data.frame(cv_acc = cv_acc, auc = auc, seed = seed)
  return(res)
}


## microbiome + demographics
data = preprocessed_data %>% select(-pnrs) %>% 
  dplyr::rename(c("adtype_IntrinsicAD" = "adtype_Intrinsic AD", 
                  "agegr1_over65years" = "agegr1_> 65 years")) %>% 
  as.data.frame()

### glmnet
learner = makeLearner("classif.glmnet", predict.type = "prob")
ParamSpace = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 0),
  makeNumericParam("lambda", lower = 0.001, upper = 1)
)
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_full_glm.csv")


### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combined_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_full_svm.csv")

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 20)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_full_rf.csv")

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_full_xgb.csv")



## modeling only demographics results
data = preprocessed_data %>% 
  rename(c("adtype_IntrinsicAD" = "adtype_Intrinsic AD", 
           "agegr1_over65years" = "agegr1_> 65 years")) %>% 
  select("acountry_Austria", "acountry_France", "acountry_Germany", "acountry_Poland", 
         "adtype_IntrinsicAD", "agegr1_over65years", "asex_Male", "nrs_w4")

### glmnet
learner = makeLearner("classif.glmnet", predict.type = "prob")
ParamSpace = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 0),
  makeNumericParam("lambda", lower = 0.001, upper = 1)
)
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_base_glm.csv")


### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_base_svm.csv")

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 20)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_base_rf.csv")

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)

results = lapply(1:50, function(x){
  iter_ml_results(data,
                  learner, 
                  ParamSpace, 
                  randSearch, 
                  cvForTuning, 
                  training_fraction = 0.8,
                  seed = x)
})

combined_results = results %>% data.table::rbindlist(.)
write_csv(combine_results, 
          file = "/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/mlr_nrsw4_base_xgb.csv")

setwd("/Users/zrayw/Desktop/PN_microbiome/analysis/machine learning data/")
base_glm = read_csv(file = "mlr_nrsw4_base_glm.csv") %>% 
  mutate(method = "glm",
         covariate = "Base")
base_svm = read.csv(file = "mlr_nrsw4_base_svm.csv") %>% 
  mutate(method = "svm",
         covariate = "Base")
base_rf = read.csv(file = "mlr_nrsw4_base_rf.csv") %>% 
  mutate(method = "rf",
         covariate = "Base")
base_xgb = read.csv(file = "mlr_nrsw4_base_xgb.csv") %>% 
  mutate(method = "xgb",
         covariate = "Base")

full_glm = read_csv(file = "mlr_nrsw4_full_glm.csv") %>% 
  mutate(method = "glm",
         covariate = "Full")
full_svm = read.csv(file = "mlr_nrsw4_full_svm.csv") %>% 
  mutate(method = "svm",
         covariate = "Full")
full_rf = read.csv(file = "mlr_nrsw4_full_rf.csv") %>% 
  mutate(method = "rf",
         covariate = "Full")
full_xgb = read.csv(file = "mlr_nrsw4_full_xgb.csv") %>% 
  mutate(method = "xgb",
         covariate = "Full")
combine_results = data.table::rbindlist(list(base_glm, 
                                             base_svm, 
                                             base_rf, 
                                             base_xgb,
                                             full_glm, 
                                             full_svm, 
                                             full_rf, 
                                             full_xgb))
combine_results %>% 
  group_by(method, covariate) %>% 
  summarize(auc = mean(auc), .groups = "drop") %>% 
  mutate(method1 = rep(1:4, each = 2)) %>% 
  ggplot(aes(x = method1, y = auc, color = covariate)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.5, linewidth = 1.2, lty = 2) +
  scale_y_continuous(breaks = 1:5/5, 
                     limits = c(0.2, 1.0)) +
  scale_x_continuous(breaks = 1:4,
                   labels = c("L2-logistic", "RF", "SVM", "XGB")) +
  scale_color_manual(name = NULL,
                     breaks = c("Base", "Full"),
                     values = c("red", "blue")) +
  labs(x = NULL, y = "Mean AUC") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2), 
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 16, 
                                   color = "black", face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"))
ggsave("/Users/zrayw/Desktop/PN_microbiome/analysis/manuscript figures/mlr_nrsw41.tiff", width = 4, height = 3.5, dpi = 500)

