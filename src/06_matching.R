# ============================================================================
# 06_matching.R
# Matched analysis: model-concordant vs model-discordant care
# Run after: 01_clean.R through 05_counterfactual.R
# ============================================================================
# 
# Estimand: Compare observed outcomes between patients who happened to receive
# the model-recommended therapy combination (concordant) vs. those who received
# a different combination (discordant), matched on baseline characteristics.
#
# Note: This is an exploratory, descriptive analysis. Concordance is incidental
# — therapists did not consult the model. This analysis is retained for
# directional evidence only and should not be interpreted as a causal estimate
# of deployment benefit.
# ============================================================================

set.seed(123)

cat("================================================================\n")
cat("MATCHED ANALYSIS: CONCORDANT VS DISCORDANT CARE\n")
cat("================================================================\n\n")

# ============================================================================
# 1. IDENTIFY MODEL-CONCORDANT EPISODES IN TEST SET
# ============================================================================

cat("--- 1. Identifying concordant vs discordant episodes ---\n\n")

therapy_cols <- c("act", "cbt", "dbt", "motivational_interviewing",
                  "mindfulness", "stages_of_change", "family_systems")

# Observed combination string for each test patient
observed_combos <- apply(
  X_test[, therapy_cols], 1,
  function(x) paste(x, collapse = "-")
)

# Model-recommended combination string for each test patient
# optimal_combos was saved in 05_counterfactual.R
recommended_combos <- apply(
  optimal_combos, 1,
  function(x) paste(as.integer(x), collapse = "-")
)

# Concordant = therapist chose what model would recommend
concordant <- as.integer(observed_combos == recommended_combos)

cat(sprintf("Test set: N = %d\n", length(concordant)))
cat(sprintf("Model-concordant episodes: %d (%.1f%%)\n",
            sum(concordant), 100 * mean(concordant)))
cat(sprintf("Model-discordant episodes: %d (%.1f%%)\n",
            sum(1 - concordant), 100 * mean(1 - concordant)))

# ============================================================================
# 2. BUILD MATCHING DATASET
# ============================================================================

cat("\n--- 2. Building matching dataset ---\n\n")

# Patient features for matching — use numeric patient clinical features only
# (no therapist/facility — those would be tautologically related to concordance)
patient_numeric_features <- c(
  "total_score", "risk_high_initial", "adolescent", "male", "days_first_srs"
)
patient_numeric_features <- patient_numeric_features[
  patient_numeric_features %in% names(test_data)
]

# Add diagnosis group as factor
match_df <- test_data %>%
  mutate(
    concordant   = concordant,
    outcome      = y_test,
    pred_prob    = pred_test,
    gain         = personalization_gains
  ) %>%
  select(all_of(c(patient_numeric_features, "dx_group", "program",
                  "concordant", "outcome", "pred_prob", "gain"))) %>%
  mutate(
    dx_group = as.factor(dx_group),
    program  = as.factor(program)
  )

# Drop any rows with missing values in matching variables
match_vars <- c(patient_numeric_features, "dx_group", "program")
match_df_complete <- match_df %>%
  filter(complete.cases(select(., all_of(match_vars))))

cat(sprintf("Complete cases for matching: %d\n", nrow(match_df_complete)))
cat(sprintf("Concordant: %d | Discordant: %d\n",
            sum(match_df_complete$concordant),
            sum(1 - match_df_complete$concordant)))

# ============================================================================
# 3. CHECK BASELINE BALANCE BEFORE MATCHING
# ============================================================================

cat("\n--- 3. Baseline balance before matching ---\n\n")

balance_before <- match_df_complete %>%
  group_by(concordant) %>%
  summarise(
    n                = n(),
    pct_high_risk    = round(100 * mean(risk_high_initial, na.rm = TRUE), 1),
    mean_cssrs       = round(mean(total_score, na.rm = TRUE), 1),
    pct_adolescent   = round(100 * mean(adolescent, na.rm = TRUE), 1),
    pct_male         = round(100 * mean(male, na.rm = TRUE), 1),
    pct_rtc          = round(100 * mean(program == "RTC"), 1),
    mean_pred_prob   = round(mean(pred_prob), 3),
    obs_improve_rate = round(100 * mean(outcome), 1),
    .groups = "drop"
  ) %>%
  mutate(group = ifelse(concordant == 1, "Concordant", "Discordant")) %>%
  select(group, everything(), -concordant)

cat("Balance before matching:\n")
print(balance_before, row.names = FALSE)

# ============================================================================
# 4. PROPENSITY SCORE MATCHING (1:1 NEAREST NEIGHBOR)
# ============================================================================

cat("\n--- 4. Propensity score matching (1:1 nearest neighbor) ---\n\n")

# Build matching formula from available patient features
match_formula_vars <- match_vars[match_vars %in% names(match_df_complete)]
match_formula <- as.formula(
  paste("concordant ~", paste(match_formula_vars, collapse = " + "))
)
cat("Matching formula:\n")
print(match_formula)

# Run matching
match_out <- tryCatch(
  matchit(
    formula   = match_formula,
    data      = match_df_complete,
    method    = "nearest",
    distance  = "logit",
    ratio     = 1,
    replace   = FALSE,
    caliper   = 0.2,   # 0.2 SD of propensity score
    std.caliper = TRUE
  ),
  error = function(e) {
    cat(sprintf("Matching with caliper failed: %s\nRetrying without caliper...\n",
                e$message))
    matchit(
      formula  = match_formula,
      data     = match_df_complete,
      method   = "nearest",
      distance = "logit",
      ratio    = 1,
      replace  = FALSE
    )
  }
)

cat("\nMatching summary:\n")
print(summary(match_out, un = FALSE))

# Extract matched data
matched_data <- match.data(match_out)

n_concordant  <- sum(matched_data$concordant == 1)
n_discordant  <- sum(matched_data$concordant == 0)

cat(sprintf("\nMatched sample: %d concordant, %d discordant\n",
            n_concordant, n_discordant))

# ============================================================================
# 5. BALANCE AFTER MATCHING
# ============================================================================

cat("\n--- 5. Balance after matching ---\n\n")

balance_after <- matched_data %>%
  group_by(concordant) %>%
  summarise(
    n                = n(),
    pct_high_risk    = round(100 * mean(risk_high_initial, na.rm = TRUE), 1),
    mean_cssrs       = round(mean(total_score, na.rm = TRUE), 1),
    pct_adolescent   = round(100 * mean(adolescent, na.rm = TRUE), 1),
    pct_male         = round(100 * mean(male, na.rm = TRUE), 1),
    pct_rtc          = round(100 * mean(program == "RTC"), 1),
    mean_pred_prob   = round(mean(pred_prob), 3),
    obs_improve_rate = round(100 * mean(outcome), 1),
    .groups = "drop"
  ) %>%
  mutate(group = ifelse(concordant == 1, "Concordant", "Discordant")) %>%
  select(group, everything(), -concordant)

cat("Balance after matching:\n")
print(balance_after, row.names = FALSE)

# ============================================================================
# 6. PRIMARY OUTCOME COMPARISON IN MATCHED SAMPLE
# ============================================================================

cat("\n--- 6. Outcome comparison in matched sample ---\n\n")

concordant_outcomes   <- matched_data$outcome[matched_data$concordant == 1]
discordant_outcomes   <- matched_data$outcome[matched_data$concordant == 0]

rate_concordant  <- mean(concordant_outcomes)
rate_discordant  <- mean(discordant_outcomes)
risk_diff        <- rate_concordant - rate_discordant
risk_diff_pp     <- risk_diff * 100

# Bootstrap CI for risk difference
set.seed(123)
n_boot <- 1000
boot_diffs <- numeric(n_boot)

for (b in 1:n_boot) {
  boot_idx <- sample(seq_len(n_concordant), replace = TRUE)
  boot_diffs[b] <- mean(concordant_outcomes[boot_idx]) -
    mean(discordant_outcomes[boot_idx])
}

ci_lower_pp <- quantile(boot_diffs, 0.025) * 100
ci_upper_pp <- quantile(boot_diffs, 0.975) * 100

# Chi-squared test
ct <- table(matched_data$concordant, matched_data$outcome)
chi_test <- chisq.test(ct)

cat(sprintf("Concordant improvement rate:   %.1f%% (n = %d)\n",
            rate_concordant * 100, n_concordant))
cat(sprintf("Discordant improvement rate:   %.1f%% (n = %d)\n",
            rate_discordant * 100, n_discordant))
cat(sprintf("Risk difference:               %.1f pp\n", risk_diff_pp))
cat(sprintf("95%% CI (bootstrap):            %.1f to %.1f pp\n",
            ci_lower_pp, ci_upper_pp))
cat(sprintf("Chi-squared p-value:           %.3f\n", chi_test$p.value))
cat(sprintf("\nNNT: %.0f\n", 1 / risk_diff))


# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

cat("\n--- 8. Saving results ---\n\n")

matching_results <- list(
  concordance_rate   = mean(concordant),
  n_concordant_full  = sum(concordant),
  n_discordant_full  = sum(1 - concordant),
  n_matched_pairs    = n_concordant,
  rate_concordant    = rate_concordant,
  rate_discordant    = rate_discordant,
  risk_diff_pp       = risk_diff_pp,
  ci_lower_pp        = ci_lower_pp,
  ci_upper_pp        = ci_upper_pp,
  chi_p              = chi_test$p.value,
  nnt                = 1 / risk_diff,
  balance_before     = balance_before,
  balance_after      = balance_after,
  matched_data       = matched_data
)

saveRDS(matching_results, "outputs/models/matching_results.rds")

# Summary table
results_table <- data.frame(
  Metric  = c("N matched pairs",
              "Concordant improvement rate",
              "Discordant improvement rate",
              "Risk difference (pp)",
              "95% CI lower (pp)",
              "95% CI upper (pp)",
              "p-value (chi-squared)",
              "NNT"),
  Value   = c(n_concordant,
              paste0(round(rate_concordant * 100, 1), "%"),
              paste0(round(rate_discordant * 100, 1), "%"),
              round(risk_diff_pp, 1),
              round(ci_lower_pp, 1),
              round(ci_upper_pp, 1),
              round(chi_test$p.value, 3),
              round(1 / risk_diff))
)

results_table
