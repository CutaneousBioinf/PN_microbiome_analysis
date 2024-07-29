# This script performs machine learning methods on baseline lesional samples in treatment groups to predict week 4 responding status

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

demographics = map_sample %>% 
  filter(trt01p == "Nemolizumab 0.5mg/kg" & lesional_visit == "lesional_V3") %>%
  select(sample, agegr1, asex, adtype, acountry, pnrs, nrs_w4)
genus_re_abund = shared %>% 
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>% 
  group_by(sample) %>% 
  mutate(re_abund = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  filter(otu %in% common_genus) %>% 
  pivot_wider(names_from = "otu", values_from = "re_abund")
full_data = inner_join(demographics, genus_re_abund, by = "sample") %>% select(-sample)

### preprocess data
preprocessed_data = preprocess_data(full_data, outcome_colname = "nrs_w4")$dat_transformed


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


## using microbiome + demographics to predict reponse status at week 4 (Full model)

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
res = lapply(
  1:50, 
  function(x){
  iter_ml_results(
    data, learner, ParamSpace, randSearch,
    cvForTuning, training_fraction = 0.8, seed = x)
})
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_full_glm.csv")

### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_full_svm.csv")

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 20)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_full_rf.csv")

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_full_xgb.csv")


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
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_base_glm.csv")

### svm radial model
learner = makeLearner("classif.svm", predict.type = "prob")
ParamSpace = makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 10),
  makeNumericParam("gamma", lower = 0.1, 10))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_base_svm.csv")

### random forest
learner = makeLearner("classif.randomForest", predict.type = "prob")
ParamSpace = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 20)
)
randSearch = makeTuneControlRandom(maxit = 20)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_base_rf.csv")

### xgboost
learner = makeLearner("classif.xgboost", predict.type = "prob")
ParamSpace = makeParamSet(makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                          makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                          makeNumericParam("subsample",lower = 0.5,upper = 1), 
                          makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
res = lapply(
  1:50, 
  function(x){
    iter_ml_results(
      data, learner, ParamSpace, randSearch,
      cvForTuning, training_fraction = 0.8, seed = x)
  })
res = res %>% data.table::rbindlist(.)
write_csv(res, file = "machine learning data/mlr_nrsw4_base_xgb.csv")


## plot auc for different methods and base/full model

base_glm = read_csv(file = "machine learning data/mlr_nrsw4_base_glm.csv") %>% 
  mutate(method = "glm", covariate = "Base")
base_svm = read.csv(file = "machine learning data/mlr_nrsw4_base_svm.csv") %>% 
  mutate(method = "svm", covariate = "Base")
base_rf = read.csv(file = "machine learning data/mlr_nrsw4_base_rf.csv") %>% 
  mutate(method = "rf", covariate = "Base")
base_xgb = read.csv(file = "machine learning data/mlr_nrsw4_base_xgb.csv") %>% 
  mutate(method = "xgb", covariate = "Base")
full_glm = read_csv(file = "machine learning data/mlr_nrsw4_full_glm.csv") %>% 
  mutate(method = "glm", covariate = "Full")
full_svm = read.csv(file = "machine learning data/mlr_nrsw4_full_svm.csv") %>% 
  mutate(method = "svm", covariate = "Full")
full_rf = read.csv(file = "machine learning data/mlr_nrsw4_full_rf.csv") %>% 
  mutate(method = "rf", covariate = "Full")
full_xgb = read.csv(file = "machine learning data/mlr_nrsw4_full_xgb.csv") %>% 
  mutate(method = "xgb", covariate = "Full")
combined_results = data.table::rbindlist(
  list(base_glm, base_svm, base_rf, base_xgb,
       full_glm, full_svm, full_rf, full_xgb
       )
  )

### wilcox test on full model auc between different methods
glm = full_glm$auc
svm = full_svm$auc
xgb = full_xgb$auc
rf = full_rf$auc
wilcox.test(glm, rf, paired = TRUE)
wilcox.test(glm, svm, paired = TRUE)
wilcox.test(glm, xgb, paired = TRUE)

### wilcox test on different methods auc between full/base model
glm1 = full_glm$auc
glm2 = base_glm$auc
wilcox.test(glm1, glm2, paired = TRUE)
svm1 = full_svm$auc
svm2 = base_svm$auc
wilcox.test(svm1, svm2, paired = TRUE)
rf1 = full_rf$auc
rf2 = base_rf$auc
wilcox.test(rf1, rf2, paired = TRUE)
xgb1 = full_xgb$auc
xgb2 = base_xgb$auc
wilcox.test(xgb1, xgb2, paired = TRUE)

combined_results %>% 
  group_by(method, covariate) %>% 
  summarize(auc = mean(auc), .groups = "drop") %>% 
  mutate(method1 = rep(1:4, each = 2)) %>% 
  ggplot(aes(x = method1, y = auc, color = covariate)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.5, linewidth = 1.2, lty = 2) +
  scale_y_continuous(breaks = 0:5/5, limits = c(0, 1.0)) +
  scale_x_continuous(breaks = 1:4, labels = c("L2-logistic", "RF", "SVM", "XGB")) +
  scale_color_manual(
    name = NULL, breaks = c("Base", "Full"),
    values = c("red", "blue")
    ) +
  labs(x = NULL, y = "Mean AUC") +
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.2), 
    legend.text = element_text(size = 16, face = "bold", family = "Times"),
    axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold", family = "Times"),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times")
    )
ggsave("manuscript figures/line.mean.auc_mlr_nrsw4.tiff", width = 4, height = 3.5, dpi = 500)

combined_results %>% 
  ggplot(aes(x = method, y = auc, color = covariate)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linewidth = 1.2, lty = 2) +
  scale_y_continuous(breaks = 0:5/5, limits = c(0, 1.0)) +
  scale_x_discrete(
    breaks = c("glm", "rf", "svm", "xgb"),
    labels = c("L2-logistic", "RF", "SVM", "XGB")
    ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Base", "Full"),
    values = c("red", "blue")
    ) +
  labs(x = NULL, y = "AUC") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 16, face = "bold", family = "Times"),
    axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold", family = "Times"),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times")
    )
ggsave("manuscript figures/box.mean.auc_mlr_nrsw4.tiff", width = 5, height = 3.5, dpi = 500)


## Use ridge logistic to fit the final model

learner = makeLearner("classif.glmnet", predict.type = "prob")
ParamSpace = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 0),
  makeNumericParam("lambda", lower = 0.001, upper = 1)
)
randSearch = makeTuneControlRandom(maxit = 200)
cvForTuning = makeResampleDesc("RepCV", folds = 5, reps = 50)
Task = makeClassifTask(data = data, target = "nrs_w4")
set.seed(1)
tunedPars = tuneParams(
  learner, task = Task, resampling = cvForTuning,
  par.set = ParamSpace, control = randSearch, 
  measures = acc
  )
tuned = setHyperPars(learner, par.vals = tunedPars$x)
tunedModel = mlr::train(tuned, Task)
importance = generateFeatureImportanceData(
  Task, "permutation.importance",learner
  )
importance$res %>% 
  t() %>% 
  as_tibble(rownames = "features") %>%
  filter(str_detect(features, "Phylo")) %>% 
  mutate(
    abs = abs(mmce),
    names = common_genus_name
    ) %>% 
  top_n(abs, n = 5) %>% 
  ggplot(aes(x = mmce, y = reorder(names, abs))) +
  geom_bar(stat = "identity", fill = "dark grey") +
  scale_x_continuous(position = "top") +
  labs(x = "Feature Importance",
       y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold", family = "Times"),
    axis.text.y = element_text(size = 12, face = "bold", family = "Times")
    )
ggsave("manuscript figures/bar.importance_nrsw4.tiff", width = 4.5, height = 2.5)





