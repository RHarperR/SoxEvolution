Random Forest and rfPermute
================
2026-07-15

``` r
library(dplyr)
```

    ## 
    ## 载入程辑包：'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(randomForest)
```

    ## randomForest 4.7-1.1

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## 载入程辑包：'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(caret)
```

    ## 载入需要的程辑包：ggplot2

    ## 
    ## 载入程辑包：'ggplot2'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     margin

    ## 载入需要的程辑包：lattice

``` r
library(rfPermute)
```

    ## Welcome to rfPermute v2.5.1
    ## See rfPermuteTutorial() for a guide to the package.

    ## 
    ## 载入程辑包：'rfPermute'

    ## The following object is masked from 'package:caret':
    ## 
    ##     confusionMatrix

``` r
ecology   <- read.table('21_ecology_factor.txt',    header = TRUE, row.names = NULL)
geography <- read.table('21_geography_factor.txt',   header = TRUE, row.names = NULL)
phy_spe   <- read.table('21_phylo.spe_factor.txt',   header = TRUE, row.names = NULL)
strategy  <- read.table('21_metabolism_rDsr.txt',    header = TRUE, row.names = NULL)
```

``` r
combine <- merge(ecology, geography, by = 'accession')
combine <- merge(combine, phy_spe,  by = 'accession')
combine <- merge(combine, strategy, by = 'accession')
combine <- combine[, -1]
combine[, c(1, 11)] <- lapply(combine[, c(1, 11)], as.factor)
```

``` r
set.seed(1)
trainIndex <- createDataPartition(combine$strategy, p = 0.7, list = FALSE)
df_train <- combine[trainIndex, ]
df_test  <- combine[-trainIndex, ]
```

``` r
set.seed(12346)
tune_grid <- expand.grid(
  .mtry = seq(2, ncol(df_train) - 1, by = 1) # 
)

# Set initial ntree value
train_control <- trainControl(
  method = "cv",
  number = 10,
  search = "grid"
)

# 
rf_model <- train(
  strategy ~ .,
  data = df_train,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "Accuracy",
  ntree = 500, 
  importance = TRUE,
)

rf_model$bestTune
```

    ##   mtry
    ## 2    3

``` r
best_mtry <- as.numeric(rf_model$bestTune)
ntree_values <- seq(1000, 1500, by = 100)
set.seed(12346)
results <- data.frame(ntree = integer(), Accuracy = double())
for (nt in ntree_values) {
  model <- randomForest(
    strategy ~ .,
    data = df_train,
    mtry = best_mtry,  # Use the best mtry from Step 1
    ntree = nt,
    importance = TRUE,
  )
  results <- rbind(results, data.frame(ntree = nt, Accuracy = mean(model$err.rate[, "OOB"])))
}

# Find best ntree
best_ntree <- results[which.min(results$Accuracy), "ntree"]
best_ntree
```

    ## [1] 1400

# Random forest

``` r
set.seed(123)
final_rf <- randomForest(
  strategy ~ .,
  data       = df_train,
  mtry       = 3,
  ntree      = 1400,
  importance = TRUE
)

print(final_rf)
```

    ## 
    ## Call:
    ##  randomForest(formula = strategy ~ ., data = df_train, mtry = 3,      ntree = 1400, importance = TRUE) 
    ##                Type of random forest: classification
    ##                      Number of trees: 1400
    ## No. of variables tried at each split: 3
    ## 
    ##         OOB estimate of  error rate: 19.72%
    ## Confusion matrix:
    ##        both others rDsr soxCD class.error
    ## both      5      1    4    65  0.93333333
    ## others    1    274   46   155  0.42436975
    ## rDsr      4     19  234    72  0.28875380
    ## soxCD     3     57   32  1355  0.06357982

# Model evaluation

``` r
train_pred <- predict(final_rf, df_train)
train_cm   <- caret::confusionMatrix(train_pred, df_train$strategy)

test_pred  <- predict(final_rf, df_test)
test_cm    <- caret::confusionMatrix(test_pred, df_test$strategy)
```

## Training and testing

``` r
train_cm
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction both others rDsr soxCD
    ##     both     75      0    0     0
    ##     others    0    476    0     0
    ##     rDsr      0      0  329     0
    ##     soxCD     0      0    0  1447
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9984, 1)
    ##     No Information Rate : 0.6218     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: both Class: others Class: rDsr Class: soxCD
    ## Sensitivity              1.00000        1.0000      1.0000       1.0000
    ## Specificity              1.00000        1.0000      1.0000       1.0000
    ## Pos Pred Value           1.00000        1.0000      1.0000       1.0000
    ## Neg Pred Value           1.00000        1.0000      1.0000       1.0000
    ## Prevalence               0.03223        0.2046      0.1414       0.6218
    ## Detection Rate           0.03223        0.2046      0.1414       0.6218
    ## Detection Prevalence     0.03223        0.2046      0.1414       0.6218
    ## Balanced Accuracy        1.00000        1.0000      1.0000       1.0000

## Testing

``` r
test_cm
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction both others rDsr soxCD
    ##     both      0      0    3     4
    ##     others    1    114   10    27
    ##     rDsr      2     23   95    18
    ##     soxCD    28     66   33   570
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7837          
    ##                  95% CI : (0.7568, 0.8089)
    ##     No Information Rate : 0.6227          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.5777          
    ##                                           
    ##  Mcnemar's Test P-Value : 4.496e-08       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: both Class: others Class: rDsr Class: soxCD
    ## Sensitivity             0.000000        0.5616     0.67376       0.9208
    ## Specificity             0.992731        0.9520     0.94959       0.6613
    ## Pos Pred Value          0.000000        0.7500     0.68841       0.8178
    ## Neg Pred Value          0.968592        0.8943     0.94626       0.8350
    ## Prevalence              0.031187        0.2042     0.14185       0.6227
    ## Detection Rate          0.000000        0.1147     0.09557       0.5734
    ## Detection Prevalence    0.007042        0.1529     0.13883       0.7012
    ## Balanced Accuracy       0.496366        0.7568     0.81167       0.7911

# rfPermute

``` r
set.seed(123)
factor_rfP <- rfPermute(
  strategy ~ .,
  data       = combine,
  mtry       = 3,
  ntree      = 1400,
  importance = TRUE,
  num.cores  = 10,
  num.rep = 400
)

factor_rfP
```

    ## An rfPermute model
    ## 
    ##                Type of random forest: classification 
    ##                      Number of trees: 1400 
    ## No. of variables tried at each split: 3 
    ##        No. of permutation replicates: 400 
    ##                           Start time: 2026-07-15 11:47:20 
    ##                             End time: 2026-07-15 12:07:08 
    ##                             Run time: 19.8 mins 
    ## 
    ##         both others rDsr soxCD pct.correct LCI_0.95 UCI_0.95
    ## both       2      0   10    94        1.89    0.229     6.65
    ## others     1    382   68   228       56.26   52.434    60.03
    ## rDsr       4     33  324   109       68.94   64.537    73.10
    ## soxCD      5     66   48  1947       94.24   93.147    95.21
    ## Overall   NA     NA   NA    NA       79.95   78.543    81.30

## Variable importance and significance

``` r
importance_factor_scale <- data.frame(
  importance(factor_rfP, scale = TRUE),
  check.names = FALSE
)

p_values <- importance_factor_scale[, 12]
adj_p_values <- p.adjust(p_values, method = "bonferroni")

df.gini <- data.frame(
  MeanDecreaseGini = importance_factor_scale[, 11],
  Adj_P_Value = adj_p_values,
  row.names = rownames(importance_factor_scale)
)

print(df.gini)
```

    ##             MeanDecreaseGini Adj_P_Value
    ## pH                 236.72045  0.17456359
    ## spe_PC3            291.37179  0.02493766
    ## spe_PC2            307.17302  0.02493766
    ## Temperature        177.55261  1.00000000
    ## spe_PC1            254.23205  0.02493766
    ## Oxygen              58.47345  0.02493766
    ## Salinity           112.04108  1.00000000
    ## geo_PC3            131.66291  1.00000000
    ## geo_PC2            131.74477  1.00000000
    ## geo_PC1            125.80138  1.00000000

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 26200)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Chinese (Simplified)_China.936 
    ## [2] LC_CTYPE=Chinese (Simplified)_China.936   
    ## [3] LC_MONETARY=Chinese (Simplified)_China.936
    ## [4] LC_NUMERIC=C                              
    ## [5] LC_TIME=Chinese (Simplified)_China.936    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] rfPermute_2.5.1      caret_6.0-93         lattice_0.22-7      
    ## [4] ggplot2_4.0.3        randomForest_4.7-1.1 dplyr_1.1.2         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10          lubridate_1.9.2      listenv_0.9.0       
    ##  [4] class_7.3-21         digest_0.6.31        ipred_0.9-13        
    ##  [7] foreach_1.5.2        parallelly_1.34.0    R6_2.6.1            
    ## [10] plyr_1.8.8           hardhat_1.2.0        stats4_4.1.3        
    ## [13] e1071_1.7-13         evaluate_0.20        pillar_1.11.1       
    ## [16] rlang_1.1.0          rstudioapi_0.14      data.table_1.14.8   
    ## [19] rpart_4.1.19         Matrix_1.5-3         rmarkdown_2.29      
    ## [22] splines_4.1.3        gower_1.0.1          stringr_1.6.0       
    ## [25] S7_0.2.2             proxy_0.4-27         compiler_4.1.3      
    ## [28] xfun_0.49            pkgconfig_2.0.3      globals_0.16.2      
    ## [31] htmltools_0.5.4      nnet_7.3-18          tidyselect_1.2.1    
    ## [34] tibble_3.2.1         prodlim_2019.11.13   codetools_0.2-19    
    ## [37] future_1.32.0        withr_3.0.3          MASS_7.3-58.3       
    ## [40] recipes_1.0.5        ModelMetrics_1.2.2.2 grid_4.1.3          
    ## [43] nlme_3.1-162         gtable_0.3.6         lifecycle_1.0.5     
    ## [46] magrittr_2.0.3       pROC_1.18.0          scales_1.4.0        
    ## [49] future.apply_1.10.0  cli_3.6.1            stringi_1.7.12      
    ## [52] farver_2.1.1         reshape2_1.4.4       timeDate_4022.108   
    ## [55] generics_0.1.4       vctrs_0.6.1          lava_1.7.2.1        
    ## [58] RColorBrewer_1.1-3   iterators_1.0.14     tools_4.1.3         
    ## [61] dichromat_2.0-0.1    glue_1.6.2           purrr_1.0.1         
    ## [64] parallel_4.1.3       fastmap_1.1.1        survival_3.5-3      
    ## [67] yaml_2.3.7           timechange_0.2.0     knitr_1.49          
    ## [70] swfscMisc_1.6
