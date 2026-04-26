# ============================================================================
# 05_counterfactual.R
# Counterfactual Analysis for Personalized Therapy Selection
# ============================================================================
#
# For each patient in the validation set, this script:
#   1. Computes predicted improvement under each of the 15 most common
#      therapy combinations (holding patient features constant)
#   2. Identifies the model-recommended combination (highest predicted prob)
#   3. Calculates the personalization gain (recommended - observed)
#   4. Summarizes gains by baseline risk quartile, initial risk level, etc.
#   5. Exports all results for manuscript
#
# Key outputs:
#   - personalization_gains: vector of per-patient predicted gains
#   - optimal_combos: matrix of recommended therapy indicators
#   - personalization_viz_data: data frame for Figure 4
#   - manuscript_results: comprehensive results list (saved as .rds)
#
# Depends on: 03_outcome_model.R, 04_doubly_robust.R
# ============================================================================

library(dplyr)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("COUNTERFACTUAL ANALYSIS\n")
cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# 1. IDENTIFY TOP THERAPY COMBINATIONS
# ============================================================================

therapy_cols <- therapy_features[therapy_features %in% colnames(X_test)]

combo_matrix  <- X_test[, therapy_cols]
combo_strings <- apply(combo_matrix, 1, paste, collapse = "-")
combo_counts  <- table(combo_strings)
top_combos    <- names(sort(combo_counts, decreasing = TRUE))[1:min(15, length(combo_counts))]

cat(sprintf("  Evaluating %d therapy combinations\n", length(top_combos)))
cat(sprintf("  Coverage: %.1f%% of test episodes\n",
            100 * sum(combo_strings %in% top_combos) / length(combo_strings)))

# ============================================================================
# 2. COMPUTE PERSONALIZATION GAINS
# ============================================================================

personalization_gains <- numeric(nrow(X_test))
optimal_combos <- matrix(0, nrow = nrow(X_test), ncol = length(therapy_cols))
colnames(optimal_combos) <- therapy_cols

pb <- txtProgressBar(min = 0, max = nrow(X_test), style = 3)

for (i in 1:nrow(X_test)) {
  patient_features <- X_test[i, ]
  baseline_prob    <- pred_test[i]
  
  best_prob  <- baseline_prob
  best_combo <- patient_features[therapy_cols]
  
  for (combo_str in top_combos) {
    combo_vals <- as.numeric(strsplit(combo_str, "-")[[1]])
    
    cf_features <- patient_features
    cf_features[therapy_cols] <- combo_vals
    cf_matrix <- xgb.DMatrix(matrix(cf_features, nrow = 1))
    cf_prob   <- predict(xgb_final, cf_matrix)
    
    if (cf_prob > best_prob) {
      best_prob  <- cf_prob
      best_combo <- combo_vals
    }
  }
  
  personalization_gains[i] <- best_prob - baseline_prob
  optimal_combos[i, ]      <- best_combo
  setTxtProgressBar(pb, i)
}
close(pb)

# ============================================================================
# 3. SUMMARIZE PERSONALIZATION GAINS
# ============================================================================

mean_gain          <- mean(personalization_gains)
median_gain        <- median(personalization_gains)
pct_benefit        <- mean(personalization_gains > 0.01) * 100
pct_large_benefit  <- mean(personalization_gains > 0.05) * 100

cat(sprintf("\nPersonalization Gain Summary:\n"))
cat(sprintf("  Mean gain: %.3f (%.1f%% relative improvement)\n",
            mean_gain, mean_gain / mean(pred_test) * 100))
cat(sprintf("  Median gain: %.3f\n", median_gain))
cat(sprintf("  Patients with gain >1pp: %.1f%%\n", pct_benefit))
cat(sprintf("  Patients with gain >5pp: %.1f%%\n", pct_large_benefit))
cat(sprintf("  Maximum gain: %.3f\n", max(personalization_gains)))
cat(sprintf("  NNT: %.0f\n", ifelse(mean_gain > 0, 1 / mean_gain, Inf)))

if (sum(personalization_gains > 0) > 0) {
  mean_gain_ben <- mean(personalization_gains[personalization_gains > 0])
  cat(sprintf("  Mean gain among benefiters: %.3f (%.1f%% relative)\n",
              mean_gain_ben,
              mean_gain_ben / mean(pred_test[personalization_gains > 0]) * 100))
}

# --- Visualization data for Figure 4 ---
personalization_viz_data <- data.frame(
  baseline_prob    = pred_test,
  gain             = personalization_gains * 100,
  initial_risk     = test_data$risk_level_initial,
  observed_outcome = y_test,
  gain_category    = cut(personalization_gains * 100,
                         breaks = c(-Inf, 0, 1, 5, Inf),
                         labels = c("None", "Minimal (0-1pp)",
                                    "Moderate (1-5pp)", "Large (>5pp)"),
                         include.lowest = TRUE)
)

# --- Gains by baseline probability quartile ---
summary_by_quartile <- personalization_viz_data %>%
  mutate(prob_quartile = cut(
    baseline_prob,
    breaks = quantile(baseline_prob, probs = c(0, 0.25, 0.5, 0.75, 1)),
    labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"),
    include.lowest = TRUE
  )) %>%
  group_by(prob_quartile) %>%
  summarise(
    n                = n(),
    mean_baseline    = mean(baseline_prob),
    mean_gain        = mean(gain),
    median_gain      = median(gain),
    pct_benefit_1pp  = mean(gain > 1) * 100,
    pct_benefit_5pp  = mean(gain > 5) * 100,
    max_gain         = max(gain),
    .groups = "drop"
  )

cat("\nGains by Baseline Probability Quartile:\n")
print(summary_by_quartile)

# --- Clinical impact by gain threshold ---
cat("\nClinical Impact by Gain Threshold:\n")
cat(sprintf("%-20s  %10s  %10s  %12s  %10s\n",
            "Min gain", "N switched", "% cohort", "Mean gain pp", "Obs improve"))

for (min_gain in c(0.00, 0.01, 0.02, 0.03, 0.05, 0.10)) {
  idx <- personalization_gains >= min_gain
  n_sw <- sum(idx)
  if (n_sw > 0) {
    cat(sprintf(">= %2.0f pp %15d  %10.1f  %12.1f  %10.1f\n",
                min_gain * 100, n_sw,
                100 * n_sw / length(y_test),
                mean(personalization_gains[idx]) * 100,
                100 * mean(y_test[idx])))
  }
}

# ============================================================================
# 4. EXPORT RESULTS
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("EXPORTING RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

manuscript_results <- list(
  model_type = "AIPW-weighted XGBoost",
  
  model_performance = data.frame(
    metric = c("AUC_train", "AUC_test", "Sensitivity", "Specificity",
               "PPV", "NPV", "Brier", "Brier_Platt"),
    value = c(auc_train, auc_test, sensitivity, specificity, ppv, npv,
              brier_score, brier_score_platt)
  ),
  
  feature_importance = group_importance,
  
  therapy_associations = if (length(dr_results) > 0) {
    data.frame(
      therapy  = names(dr_results),
      ate      = sapply(dr_results, function(x) x$ate),
      se       = sapply(dr_results, function(x) x$se),
      ci_lower = sapply(dr_results, function(x) x$ci_lower),
      ci_upper = sapply(dr_results, function(x) x$ci_upper),
      p_value  = sapply(dr_results, function(x) x$p_value),
      n_treated = sapply(dr_results, function(x) x$n_treated)
    ) %>% arrange(desc(ate))
  } else NULL,
  
  personalization_summary = data.frame(
    metric = c("mean_gain", "median_gain", "pct_benefit", "pct_large_benefit",
               "max_gain", "nnt"),
    value = c(mean_gain, median_gain, pct_benefit / 100, pct_large_benefit / 100,
              max(personalization_gains), ifelse(mean_gain > 0, 1 / mean_gain, NA))
  ),
  
  personalization_by_quartile = summary_by_quartile,
  
  xgboost_model = xgb_final,
  
  weight_summary = data.frame(
    dataset = c("Training", "Test"),
    mean   = c(mean(train_weights),   mean(test_weights)),
    median = c(median(train_weights), median(test_weights)),
    sd     = c(sd(train_weights),     sd(test_weights)),
    min    = c(min(train_weights),    min(test_weights)),
    max    = c(max(train_weights),    max(test_weights))
  ),
  
  predictions = data.frame(
    predicted_prob        = pred_test,
    personalization_gain  = personalization_gains,
    optimal_therapy       = apply(optimal_combos, 1, function(x)
      paste(therapy_cols[x == 1], collapse = "+"))
  )
)

saveRDS(manuscript_results, "outputs/models/aipw_xgboost_results.rds")
cat("  Saved: outputs/models/aipw_xgboost_results.rds\n")

# Counterfactual table
sink("outputs/tables/counterfactual_analysis.txt")
cat(rep("=", 70), "\n", sep = "")
cat("PERSONALIZATION GAINS BY BASELINE PROBABILITY QUARTILE\n")
cat(rep("=", 70), "\n", sep = "")
print(summary_by_quartile, width = Inf)
sink()
cat("  Saved: outputs/tables/counterfactual_analysis.txt\n")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("COUNTERFACTUAL ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")