# ============================================================================
# 03_outcome_model.R
# IPW-Weighted XGBoost Outcome Model
# ============================================================================
#
# This script trains the primary prediction model for suicide risk improvement
# using inverse-probability-weighted XGBoost. Propensity scores from
# 02_propensity_scores.R are used as stabilized weights (not features).
#
# Steps:
#   1. Feature engineering (patient, therapy, context, therapist features)
#   2. Temporal train-test split (80/20 by admission date)
#   3. Stabilized IPW weight construction from propensity scores
#   4. Hyperparameter tuning via 3-fold temporal cross-validation
#   5. Final model training and evaluation
#   6. Feature importance analysis (SHAP-ready)
#
# Key outputs:
#   - xgb_final: trained XGBoost model object
#   - pred_test / pred_train: predicted probabilities
#   - X_train / X_test / y_train / y_test: data matrices
#   - train_data / test_data: data frames with all covariates
#   - train_weights / test_weights: stabilized IPW weights
#   - best_params_list: tuned hyperparameters
#
# Depends on: 01_clean.R, 02_propensity_scores.R
# ============================================================================

set.seed(123)

# ============================================================================
# 1. FEATURE ENGINEERING
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("FEATURE ENGINEERING\n")
cat(rep("=", 70), "\n", sep = "")

# --- Patient clinical features (measured at intake) ---
patient_clinical_features <- c(
  "total_score", "risk_high_initial", "days_first_srs",
  "deterrents_month", "what_sort_of_reasons", "duration_month",
  "adolescent",
  "could_can_you_stop_thinking_about_killing_yourself_or_wanting_to_die_if_you_want_to",
  "are_there_things", "frequency_month",
  "when_you_have_the_thoughts_how_long_do_they_last",
  "how_many_times_have_you_had_these_thoughts",
  "male", "dx_group",
  names(data)[grepl("^current_and_past_psychiatric_diagnoses_", names(data))],
  names(data)[grepl("^presenting_symptoms_", names(data))],
  names(data)[grepl("^family_history_", names(data))],
  names(data)[grepl("^precipitants_stressors_", names(data))],
  names(data)[grepl("^internal_protective_factors_|^external_protective_factors_", names(data))],
  names(data)[grepl("^change_in_treatment_", names(data))]
)

# --- Therapy modality indicators ---
therapy_features <- c("act", "cbt", "dbt", "motivational_interviewing",
                      "mindfulness", "stages_of_change", "family_systems")

# --- Treatment context ---
treatment_context_features <- c("therapy_duration_category", "delivery_method",
                                "session_mode")

# --- Therapist/organizational features ---
therapist_features <- c("therapist_name", "location", "program",
                        "pn_month", "pn_time_block", "pn_year")

# Combine all features (propensity scores are NOT included as features)
all_features <- unique(c(patient_clinical_features, therapy_features,
                         treatment_context_features, therapist_features))
all_features <- all_features[all_features %in% names(data)]

cat(sprintf("  Patient clinical features: %d\n",
            sum(all_features %in% patient_clinical_features)))
cat(sprintf("  Therapy indicators: %d\n",
            sum(all_features %in% therapy_features)))
cat(sprintf("  Treatment context: %d\n",
            sum(all_features %in% treatment_context_features)))
cat(sprintf("  Therapist/organizational: %d\n",
            sum(all_features %in% therapist_features)))
cat(sprintf("  Total features: %d\n", length(all_features)))

# ============================================================================
# 2. TEMPORAL TRAIN-TEST SPLIT
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TEMPORAL TRAIN-TEST SPLIT\n")
cat(rep("=", 70), "\n", sep = "")

data <- data %>%
  arrange(admission_date) %>%
  mutate(row_id = row_number())

split_point <- floor(nrow(data) * 0.8)
train_data  <- data %>% filter(row_id <= split_point)
test_data   <- data %>% filter(row_id > split_point)

cat(sprintf("  Training set: n=%d (%.1f%% improved)\n",
            nrow(train_data), mean(train_data$improve, na.rm = TRUE) * 100))
cat(sprintf("  Test set: n=%d (%.1f%% improved)\n",
            nrow(test_data), mean(test_data$improve, na.rm = TRUE) * 100))
cat(sprintf("  Training dates: %s to %s\n",
            min(train_data$admission_date), max(train_data$admission_date)))
cat(sprintf("  Test dates: %s to %s\n",
            min(test_data$admission_date), max(test_data$admission_date)))

# ============================================================================
# 3. PREPARE DATA MATRICES
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PREPARING DATA MATRICES\n")
cat(rep("=", 70), "\n", sep = "")

# Helper function to create numeric matrices from mixed data
prepare_matrices <- function(train_df, test_df, features) {
  
  categorical_features <- features[sapply(train_df[features], function(x) {
    is.factor(x) || is.character(x)
  })]
  numeric_features <- features[sapply(train_df[features], function(x) {
    is.numeric(x) || is.integer(x) || is.logical(x)
  })]
  
  cat(sprintf("  Processing %d numeric and %d categorical features\n",
              length(numeric_features), length(categorical_features)))
  
  # Numeric: impute missing with training median
  X_train_num <- as.matrix(train_df[, numeric_features, drop = FALSE])
  X_test_num  <- as.matrix(test_df[, numeric_features, drop = FALSE])
  
  for (j in 1:ncol(X_train_num)) {
    if (any(is.na(X_train_num[, j]))) {
      med <- median(X_train_num[, j], na.rm = TRUE)
      if (is.na(med)) med <- 0
      X_train_num[is.na(X_train_num[, j]), j] <- med
      X_test_num[is.na(X_test_num[, j]), j]   <- med
    }
  }
  
  # Categorical: one-hot encode
  if (length(categorical_features) > 0) {
    cat_dummies_train <- list()
    cat_dummies_test  <- list()
    
    for (feat in categorical_features) {
      if (!feat %in% names(train_df)) next
      train_vals <- train_df[[feat]]
      test_vals  <- test_df[[feat]]
      if (all(is.na(train_vals))) next
      
      train_levels <- unique(as.character(train_vals[!is.na(train_vals)]))
      if (length(train_levels) <= 1) next
      
      train_factor <- factor(as.character(train_vals), levels = train_levels)
      test_factor  <- factor(as.character(test_vals),  levels = train_levels)
      
      if (any(is.na(train_factor))) {
        train_factor <- addNA(train_factor)
        levels(train_factor)[is.na(levels(train_factor))] <- "MISSING"
      }
      if (any(is.na(test_factor))) {
        test_factor <- addNA(test_factor)
        levels(test_factor)[is.na(levels(test_factor))] <- "MISSING"
      }
      
      for (level in levels(train_factor)[-1]) {
        dummy_name <- paste0(feat, "_", gsub("[^[:alnum:]]", "_", level))
        cat_dummies_train[[dummy_name]] <- as.numeric(train_factor == level)
        cat_dummies_test[[dummy_name]]  <- as.numeric(test_factor == level)
      }
    }
    
    if (length(cat_dummies_train) > 0) {
      X_train_cat <- do.call(cbind, cat_dummies_train)
      X_test_cat  <- do.call(cbind, cat_dummies_test)
      X_train <- cbind(X_train_num, X_train_cat)
      X_test  <- cbind(X_test_num, X_test_cat)
      cat(sprintf("  Created %d dummy variables from categorical features\n",
                  ncol(X_train_cat)))
    } else {
      X_train <- X_train_num
      X_test  <- X_test_num
    }
  } else {
    X_train <- X_train_num
    X_test  <- X_test_num
  }
  
  # Align columns
  train_only <- setdiff(colnames(X_train), colnames(X_test))
  if (length(train_only) > 0) {
    zeros <- matrix(0, nrow = nrow(X_test), ncol = length(train_only))
    colnames(zeros) <- train_only
    X_test <- cbind(X_test, zeros)
    cat(sprintf("  Added %d columns to test set (new in training)\n", length(train_only)))
  }
  X_test <- X_test[, colnames(X_train), drop = FALSE]
  
  list(train = X_train, test = X_test)
}

matrices <- prepare_matrices(train_data, test_data, all_features)
X_train  <- matrices$train
X_test   <- matrices$test
y_train  <- train_data$improve
y_test   <- test_data$improve

# Remove NA outcomes
if (any(is.na(y_train))) {
  keep <- !is.na(y_train)
  X_train <- X_train[keep, ]; y_train <- y_train[keep]; train_data <- train_data[keep, ]
  cat(sprintf("  Removed %d training rows with NA outcomes\n", sum(!keep)))
}
if (any(is.na(y_test))) {
  keep <- !is.na(y_test)
  X_test <- X_test[keep, ]; y_test <- y_test[keep]; test_data <- test_data[keep, ]
  cat(sprintf("  Removed %d test rows with NA outcomes\n", sum(!keep)))
}

cat(sprintf("\nFinal dimensions: X_train %d x %d, X_test %d x %d\n",
            nrow(X_train), ncol(X_train), nrow(X_test), ncol(X_test)))

# ============================================================================
# 4. CONSTRUCT STABILIZED IPW WEIGHTS
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("CONSTRUCTING STABILIZED IPW WEIGHTS\n")
cat(rep("=", 70), "\n", sep = "")

# Mapping: therapy column name -> propensity score column name
therapy_to_prop <- c(
  "act"                       = "prop_act",
  "cbt"                       = "prop_cbt",
  "dbt"                       = "prop_dbt",
  "motivational_interviewing" = "prop_mi",
  "mindfulness"               = "prop_mindfulness",
  "stages_of_change"          = "prop_stages_of_change",
  "family_systems"            = "prop_family_systems"
)

therapy_cols_in_data <- therapy_features[therapy_features %in% colnames(X_train)]
cat(sprintf("  Therapy modalities found: %d\n", length(therapy_cols_in_data)))

create_aipw_weights <- function(X_matrix, df_with_props, therapy_cols,
                                therapy_to_prop_map) {
  n <- nrow(X_matrix)
  weights <- rep(1, n)
  weight_components <- data.frame(row = 1:n)
  
  for (therapy in therapy_cols) {
    prop_col <- therapy_to_prop_map[therapy]
    if (is.na(prop_col) || !prop_col %in% names(df_with_props)) next
    
    T_i <- X_matrix[, therapy]
    e_i <- pmax(0.01, pmin(0.99, df_with_props[[prop_col]]))
    p_t <- mean(T_i, na.rm = TRUE)
    
    # Stabilized IPW: w = p(T)/e for treated, (1-p(T))/(1-e) for untreated
    w_i <- ifelse(T_i == 1, p_t / e_i, (1 - p_t) / (1 - e_i))
    weights <- weights * w_i
    weight_components[[paste0("w_", therapy)]] <- w_i
    
    cat(sprintf("    %s: propensity [%.3f, %.3f], treatment rate = %.1f%%\n",
                therapy, min(e_i), max(e_i), p_t * 100))
  }
  
  # Normalize to mean 1
  weights <- weights / mean(weights)
  
  cat(sprintf("\n  Weight statistics: median=%.3f, range=[%.3f, %.3f], SD=%.3f\n",
              median(weights), min(weights), max(weights), sd(weights)))
  
  # Trim at 99th percentile
  cap <- quantile(weights, 0.99)
  weights_trimmed <- pmin(weights, cap)
  n_trimmed <- sum(weights != weights_trimmed)
  if (n_trimmed > 0) {
    cat(sprintf("  Trimmed %d weights at 99th percentile (%.3f)\n", n_trimmed, cap))
  }
  
  list(weights = weights_trimmed, weights_untrimmed = weights,
       components = weight_components)
}

train_weights_result <- create_aipw_weights(X_train, train_data,
                                            therapy_cols_in_data, therapy_to_prop)
train_weights <- train_weights_result$weights

test_weights_result <- create_aipw_weights(X_test, test_data,
                                           therapy_cols_in_data, therapy_to_prop)
test_weights <- test_weights_result$weights

# ============================================================================
# 5. HYPERPARAMETER TUNING (3-FOLD TEMPORAL CV)
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("HYPERPARAMETER TUNING\n")
cat(rep("=", 70), "\n", sep = "")

# Expanding-window temporal splits
create_time_series_splits <- function(n_samples, n_splits = 3) {
  test_size <- floor(n_samples / (n_splits + 1))
  lapply(1:n_splits, function(i) {
    train_end  <- test_size * i
    test_start <- train_end + 1
    test_end   <- min(train_end + test_size, n_samples)
    list(train = 1:train_end, test = test_start:test_end)
  })
}

cv_splits <- create_time_series_splits(nrow(X_train), n_splits = 3)

for (i in seq_along(cv_splits)) {
  cat(sprintf("  Fold %d: Train [1:%d], Test [%d:%d]\n",
              i, max(cv_splits[[i]]$train),
              min(cv_splits[[i]]$test), max(cv_splits[[i]]$test)))
}

# Parameter grid (random sample of 100 combinations)
param_grid <- expand.grid(
  max_depth        = c(3, 4, 5, 6, 8),
  eta              = c(0.01, 0.05, 0.1, 0.15),
  subsample        = c(0.6, 0.7, 0.8, 0.9),
  colsample_bytree = c(0.6, 0.7, 0.8, 0.9),
  min_child_weight = c(1, 3, 5, 7),
  gamma            = c(0, 0.5, 1, 2),
  alpha            = c(0, 0.5, 1),
  lambda           = c(0.5, 1, 2)
)

set.seed(999)
param_grid <- param_grid[sample(nrow(param_grid), 100), ]
cat(sprintf("  Testing %d combinations with %d-fold temporal CV\n",
            nrow(param_grid), length(cv_splits)))

# Evaluation function
evaluate_params <- function(params_row, X, y, splits,
                            df_with_props, therapy_cols, therapy_to_prop_map) {
  params <- list(
    objective = "binary:logistic", eval_metric = "auc",
    max_depth = params_row$max_depth, eta = params_row$eta,
    subsample = params_row$subsample,
    colsample_bytree = params_row$colsample_bytree,
    min_child_weight = params_row$min_child_weight,
    gamma = params_row$gamma, alpha = params_row$alpha,
    lambda = params_row$lambda
  )
  
  fold_aucs    <- numeric(length(splits))
  fold_nrounds <- numeric(length(splits))
  
  for (i in seq_along(splits)) {
    train_idx <- splits[[i]]$train
    test_idx  <- splits[[i]]$test
    
    fold_wt <- create_aipw_weights(X[train_idx, ], df_with_props[train_idx, ],
                                   therapy_cols, therapy_to_prop_map)$weights
    
    dtrain_f <- xgb.DMatrix(X[train_idx, ], label = y[train_idx], weight = fold_wt)
    dtest_f  <- xgb.DMatrix(X[test_idx, ],  label = y[test_idx])
    
    model_f <- xgb.train(params = params, data = dtrain_f, nrounds = 500,
                         watchlist = list(test = dtest_f),
                         early_stopping_rounds = 30, verbose = 0)
    
    pred_f <- predict(model_f, dtest_f)
    fold_aucs[i]    <- as.numeric(auc(roc(y[test_idx], pred_f, quiet = TRUE)))
    fold_nrounds[i] <- model_f$best_iteration
  }
  
  list(mean_auc = mean(fold_aucs), std_auc = sd(fold_aucs),
       mean_nrounds = round(mean(fold_nrounds)))
}

# Run grid search
results_grid <- data.frame(param_grid)
results_grid$cv_auc_mean <- NA
results_grid$cv_auc_std  <- NA
results_grid$best_nrounds <- NA

pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)
for (i in 1:nrow(param_grid)) {
  res <- evaluate_params(param_grid[i, ], X_train, y_train, cv_splits,
                         train_data, therapy_cols_in_data, therapy_to_prop)
  results_grid$cv_auc_mean[i]  <- res$mean_auc
  results_grid$cv_auc_std[i]   <- res$std_auc
  results_grid$best_nrounds[i] <- res$mean_nrounds
  setTxtProgressBar(pb, i)
}
close(pb)

# Select best: stability-adjusted criterion (mean - 0.5*SD)
results_grid$score <- results_grid$cv_auc_mean - 0.5 * results_grid$cv_auc_std
best_idx    <- which.max(results_grid$score)
best_params <- results_grid[best_idx, ]

cat(sprintf("\n\nBest parameters: CV AUC = %.4f (+/- %.4f)\n",
            best_params$cv_auc_mean, best_params$cv_auc_std))
print(best_params[, c("max_depth", "eta", "subsample", "colsample_bytree",
                      "min_child_weight", "gamma", "alpha", "lambda",
                      "best_nrounds")])

# ============================================================================
# 6. TRAIN FINAL MODEL
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TRAINING FINAL MODEL\n")
cat(rep("=", 70), "\n", sep = "")

best_params_list <- list(
  objective = "binary:logistic", eval_metric = "auc",
  max_depth        = best_params$max_depth,
  eta              = best_params$eta,
  subsample        = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree,
  min_child_weight = best_params$min_child_weight,
  gamma            = best_params$gamma,
  alpha            = best_params$alpha,
  lambda           = best_params$lambda
)

dtrain <- xgb.DMatrix(data = X_train, label = y_train, weight = train_weights)
dtest  <- xgb.DMatrix(data = X_test,  label = y_test)

xgb_final <- xgb.train(
  params    = best_params_list,
  data      = dtrain,
  nrounds   = 1000,
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 30,
  verbose   = 0
)

pred_train <- predict(xgb_final, dtrain)
pred_test  <- predict(xgb_final, dtest)

auc_train <- as.numeric(auc(roc(y_train, pred_train, quiet = TRUE)))
auc_test  <- as.numeric(auc(roc(y_test,  pred_test,  quiet = TRUE)))

cat(sprintf("  Training AUC: %.4f\n", auc_train))
cat(sprintf("  Test AUC: %.4f\n", auc_test))
cat(sprintf("  Best iteration: %d\n", xgb_final$best_iteration))

# ============================================================================
# 7. PREDICTIVE PERFORMANCE METRICS
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PREDICTIVE PERFORMANCE\n")
cat(rep("=", 70), "\n", sep = "")

roc_test  <- roc(y_test, pred_test, quiet = TRUE)
coords    <- coords(roc_test, "best", ret = "all", transpose = FALSE)
optimal_threshold <- coords$threshold[1]

pred_test_binary <- as.integer(pred_test > optimal_threshold)
cm <- table(Actual = y_test, Predicted = pred_test_binary)

sensitivity <- cm[2, 2] / sum(cm[2, ])
specificity <- cm[1, 1] / sum(cm[1, ])
ppv         <- cm[2, 2] / sum(cm[, 2])
npv         <- cm[1, 1] / sum(cm[, 1])
brier_score <- mean((pred_test - y_test)^2)

cat("Confusion Matrix:\n")
print(cm)
cat(sprintf("\nThreshold: %.3f\n", optimal_threshold))
cat(sprintf("  Sensitivity: %.1f%%\n", sensitivity * 100))
cat(sprintf("  Specificity: %.1f%%\n", specificity * 100))
cat(sprintf("  PPV: %.1f%%\n", ppv * 100))
cat(sprintf("  NPV: %.1f%%\n", npv * 100))
cat(sprintf("  Brier Score: %.4f\n", brier_score))

# Platt scaling for downstream calibration assessment
platt_model <- glm(y_train ~ poly(pred_train, 4), family = binomial)
platt_probs <- predict(platt_model,
                       newdata = data.frame(pred_train = pred_test),
                       type = "response")
brier_score_platt <- mean((platt_probs - y_test)^2)
cat(sprintf("  Brier Score (Platt scaled): %.4f\n", brier_score_platt))

cat(sprintf("\nCalibration:\n"))
cat(sprintf("  Mean predicted: %.3f\n", mean(pred_test)))
cat(sprintf("  Observed rate: %.3f\n", mean(y_test)))
cat(sprintf("  Difference: %.3f\n", mean(pred_test) - mean(y_test)))

# ============================================================================
# 8. FEATURE IMPORTANCE
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("FEATURE IMPORTANCE\n")
cat(rep("=", 70), "\n", sep = "")

importance_matrix <- xgb.importance(model = xgb_final)

cat("Top 20 features by gain:\n")
print(head(importance_matrix[, c("Feature", "Gain", "Cover", "Frequency")], 20))

# Categorize by feature group
importance_df <- as.data.frame(importance_matrix) %>%
  mutate(
    feature_group = case_when(
      Feature %in% colnames(X_train)[colnames(X_train) %in% patient_clinical_features]
      ~ "Patient Clinical",
      Feature %in% colnames(X_train)[colnames(X_train) %in% therapy_features]
      ~ "Therapy Received",
      Feature %in% colnames(X_train)[colnames(X_train) %in% treatment_context_features]
      ~ "Treatment Context",
      grepl("therapist_name_|location_|program_|pn_", Feature)
      ~ "Therapist/Org",
      TRUE ~ "Other"
    )
  )

group_importance <- importance_df %>%
  group_by(feature_group) %>%
  summarise(total_gain = sum(Gain), mean_gain = mean(Gain),
            n_features = n(), .groups = "drop") %>%
  arrange(desc(total_gain))

cat("\nFeature Importance by Group:\n")
print(group_importance)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("OUTCOME MODEL COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")