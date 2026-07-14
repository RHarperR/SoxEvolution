library(dplyr)
library(randomForest)
library(caret)
ecology <- read.table('21_ecology_factor.txt', header = TRUE, row.names = NULL)
geography <- read.table('21_geography_factor.txt', header = TRUE, row.names = NULL)
phy_spe <- read.table('21_phylo.spe_factor.txt', header = TRUE, row.names = NULL)
strategy <- read.table('21_metabolism_rDsr.txt', header = TRUE, row.names = NULL)
combine <- merge(ecology,geography,by = 'accession')
combine <- merge(combine,phy_spe,by = 'accession')
combine <- merge(combine,strategy,by = 'accession')
combine <- combine[,-1]
combine[, c(1,11)] <- lapply(combine[, c(1,11)], as.factor)
head(combine,n=6)
summary(combine$strategy)

trainIndex <- createDataPartition(combine$strategy, p = .7, list = FALSE)
df_train <- combine[trainIndex,]
df_test <- combine[-trainIndex,]


final_rf <- randomForest(
  strategy ~ .,
  data = df_train,
  mtry = 3,
  ntree = 1400,
  importance = TRUE,
  #classwt = class_weights
)
print(final_rf)
train_pred <- predict(final_rf, df_train)
train_conf_matrix <- caret::confusionMatrix(train_pred, df_train$strategy)
test_pred <- predict(final_rf, df_test)
test_conf_matrix <- caret::confusionMatrix(test_pred, df_test$strategy)

library(rfPermute)
factor_rfP <- rfPermute(
  strategy ~ .,
  mtry = 3,
  ntree = 1400,
  data = combine,
  importance = TRUE,
  num.cores = 10
)

factor_rfP
importance_factor.scale <- data.frame(
  importance(factor_rfP, scale = TRUE),
  check.names = FALSE
)
importance_factor.scale
p_values <- importance_factor.scale[,12]
print(p_values)
adj_p_values <- p.adjust(p_values, method = "BH")
print(adj_p_values)
