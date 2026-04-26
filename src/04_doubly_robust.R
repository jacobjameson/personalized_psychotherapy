# ============================================================================
# 04_doubly_robust.R
# Doubly Robust Estimation of Therapy-Specific Associations
# ============================================================================
#
# This script estimates the average treatment effect (ATE) of each therapy
# modality on suicide risk improvement using augmented inverse propensity
# weighting (AIPW) with 5-fold cross-fitting.
#
# Each modality is treated as a binary indicator (received vs. not received).
# ATEs represent the bias-adjusted association between receiving each
# modality and improvement, relative to the counterfactual of not receiving
# it under the observed confounded allocation.
#
# Key outputs:
#   - dr_results: list of ATE estimates, SEs, CIs, p-values per modality
#
# Depends on: 03_outcome_model.R (provides X_train, y_train, train_data,
#             therapy_to_prop, best_params_list, therapy_features)
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("DOUBLY ROBUST ESTIMATION OF THERAPY EFFECTS\n")
cat(rep("=", 70), "\n", sep = "")

# Therapy columns present in the feature matrix
therapy_cols <- therapy_features[therapy_features %in% colnames(X_train)]
cat(sprintf("  Modalities to evaluate: %d\n", length(therapy_cols)))
cat("  Method: AIPW with 5-fold cross-fitting\n\n")

dr_results <- list()

# Create cross-fitting folds (random, not temporal — within training set)
set.seed(123)
n_folds  <- 5
fold_ids <- sample(rep(1:n_folds, length.out = nrow(X_train)))

# Features excluding all therapy indicators (to avoid conditioning on treatment)
X_without_therapy <- X_train[, !colnames(X_train) %in% therapy_cols]

for (therapy in therapy_cols) {
  
  prop_col <- therapy_to_prop[therapy]
  if (is.na(prop_col) || !prop_col %in% colnames(train_data)) next
  
  T_i <- X_train[, therapy]
  e_i <- pmax(0.01, pmin(0.99, train_data[[prop_col]]))
  
  n_treated <- sum(T_i == 1)
  n_control <- sum(T_i == 0)
  
  if (n_treated < 20 || n_control < 20) {
    cat(sprintf("  %s: SKIPPED (n_treated=%d, n_control=%d)\n",
                therapy, n_treated, n_control))
    next
  }
  
  # --- Cross-fitted outcome predictions ---
  mu_1_cf <- numeric(nrow(X_train))
  mu_0_cf <- numeric(nrow(X_train))
  
  for (fold in 1:n_folds) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    n_treated_fold <- sum(T_i[train_idx] == 1)
    n_control_fold <- sum(T_i[train_idx] == 0)
    
    # Fall back to marginal means if fold is too small
    if (n_treated_fold < 5 || n_control_fold < 5) {
      mu_1_cf[test_idx] <- mean(y_train[train_idx[T_i[train_idx] == 1]])
      mu_0_cf[test_idx] <- mean(y_train[train_idx[T_i[train_idx] == 0]])
      next
    }
    
    # Outcome model among treated
    treated_idx <- train_idx[T_i[train_idx] == 1]
    dtrain_t    <- xgb.DMatrix(X_without_therapy[treated_idx, ],
                               label = y_train[treated_idx])
    model_t     <- xgb.train(params = best_params_list, data = dtrain_t,
                             nrounds = 200, verbose = 0)
    
    # Outcome model among untreated
    control_idx <- train_idx[T_i[train_idx] == 0]
    dtrain_c    <- xgb.DMatrix(X_without_therapy[control_idx, ],
                               label = y_train[control_idx])
    model_c     <- xgb.train(params = best_params_list, data = dtrain_c,
                             nrounds = 200, verbose = 0)
    
    # Out-of-sample predictions for held-out fold
    dtest_fold <- xgb.DMatrix(X_without_therapy[test_idx, ])
    mu_1_cf[test_idx] <- predict(model_t, dtest_fold)
    mu_0_cf[test_idx] <- predict(model_c, dtest_fold)
  }
  
  # --- AIPW estimator ---
  tau_i <- mu_1_cf - mu_0_cf +
    T_i * (y_train - mu_1_cf) / e_i -
    (1 - T_i) * (y_train - mu_0_cf) / (1 - e_i)
  
  ate <- mean(tau_i)
  se  <- sd(tau_i) / sqrt(length(tau_i))
  p_value <- 2 * (1 - pnorm(abs(ate / se)))
  
  dr_results[[therapy]] <- list(
    ate       = ate,
    se        = se,
    ci_lower  = ate - 1.96 * se,
    ci_upper  = ate + 1.96 * se,
    n_treated = n_treated,
    n_control = n_control,
    p_value   = p_value
  )
  
  sig <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else
    if (p_value < 0.05) "*" else ""
  cat(sprintf("  %s: ATE = %.3f (95%% CI: %.3f to %.3f), p=%.3f %s [n=%d]\n",
              therapy, ate, ate - 1.96 * se, ate + 1.96 * se,
              p_value, sig, n_treated))
}

# ============================================================================
# SUMMARY
# ============================================================================

if (length(dr_results) > 0) {
  p_values <- sapply(dr_results, function(x) x$p_value)
  ates     <- sapply(dr_results, function(x) x$ate)
  
  cat(sprintf("\n  Therapies evaluated: %d\n", length(dr_results)))
  cat(sprintf("  Significant at p<0.05: %d\n", sum(p_values < 0.05)))
  cat(sprintf("  Significant at p<0.01: %d\n", sum(p_values < 0.01)))
  cat(sprintf("  Largest effect: %s (%.1f pp, p=%.3f)\n",
              names(which.max(ates)), max(ates) * 100,
              dr_results[[names(which.max(ates))]]$p_value))
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("DOUBLY ROBUST ESTIMATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")