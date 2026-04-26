# ============================================================================
# 08_sensitivity.R
# Sensitivity and supplemental analyses for revision response
# Run after: 01_clean.R, 02_propensity_scores.R, 03_outcome_model.R,
#            04_doubly_robust.R, 05_counterfactual.R
# ============================================================================

set.seed(123)

therapy_cols <- c("act", "cbt", "dbt", "motivational_interviewing",
                  "mindfulness", "stages_of_change", "family_systems")

therapy_labels <- c(
  "act"                       = "ACT",
  "cbt"                       = "CBT",
  "dbt"                       = "DBT",
  "motivational_interviewing"  = "MI",
  "mindfulness"               = "Mindfulness",
  "stages_of_change"          = "Stages of Change",
  "family_systems"            = "Family Systems"
)

therapy_to_prop <- c(
  "act"                       = "prop_act",
  "cbt"                       = "prop_cbt",
  "dbt"                       = "prop_dbt",
  "motivational_interviewing"  = "prop_mi",
  "mindfulness"               = "prop_mindfulness",
  "stages_of_change"          = "prop_stages_of_change",
  "family_systems"            = "prop_family_systems"
)

cat("================================================================\n")
cat("SENSITIVITY ANALYSES\n")
cat("================================================================\n\n")

# ============================================================================
# 1. IPW WEIGHT DISTRIBUTIONS 
# ============================================================================

cat("--- 1. IPW Weight Distributions ---\n\n")

weight_table <- data.frame()

for (col in names(therapy_labels)) {
  prop_col <- therapy_to_prop[col]
  if (!prop_col %in% names(train_data)) next
  
  T_i  <- X_train[, col]
  e_i  <- pmax(0.01, pmin(0.99, train_data[[prop_col]]))
  p_t  <- mean(T_i, na.rm = TRUE)
  w_i  <- ifelse(T_i == 1, p_t / e_i, (1 - p_t) / (1 - e_i))
  w_trim <- pmin(w_i, quantile(w_i, 0.99))
  
  weight_table <- rbind(weight_table, data.frame(
    Modality     = therapy_labels[col],
    Min          = round(min(w_trim), 2),
    P25          = round(quantile(w_trim, 0.25), 2),
    Median       = round(median(w_trim), 2),
    Mean         = round(mean(w_trim), 2),
    P75          = round(quantile(w_trim, 0.75), 2),
    P95          = round(quantile(w_trim, 0.95), 2),
    Max          = round(max(w_trim), 2),
    Pct_over_10  = round(100 * mean(w_i > 10), 1),
    Prop_lt_05   = round(100 * mean(train_data[[prop_col]] < 0.05, na.rm = TRUE), 1),
    Prop_gt_95   = round(100 * mean(train_data[[prop_col]] > 0.95, na.rm = TRUE), 1)
  ))
}

cat("IPW Weight Distribution Table:\n")
print(weight_table, row.names = FALSE)
write.csv(weight_table,
          "outputs/tables/sensitivity_ipw_weights.csv",
          row.names = FALSE)

# Density plot
weight_long <- data.frame()
for (col in names(therapy_labels)) {
  prop_col <- therapy_to_prop[col]
  if (!prop_col %in% names(train_data)) next
  T_i   <- X_train[, col]
  e_i   <- pmax(0.01, pmin(0.99, train_data[[prop_col]]))
  p_t   <- mean(T_i, na.rm = TRUE)
  w_i   <- ifelse(T_i == 1, p_t / e_i, (1 - p_t) / (1 - e_i))
  w_trim <- pmin(w_i, quantile(w_i, 0.99))
  weight_long <- rbind(weight_long, data.frame(
    modality = therapy_labels[col],
    weight   = w_trim
  ))
}
weight_long$modality <- factor(weight_long$modality,
                               levels = therapy_labels[c("stages_of_change","motivational_interviewing",
                                                         "cbt","dbt","mindfulness","family_systems","act")])

p_weights <- ggplot(weight_long, aes(x = weight, color = modality, fill = modality)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ modality, nrow = 2, scales = "free_y") +
  scale_x_continuous(limits = c(0, 7), breaks = c(0, 1, 2, 3, 5, 7)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Stabilized IPW value", y = "Density") + 
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(size = 8, color = "grey50", hjust = 0))

ggsave("outputs/figures/efigure_ipw_weights.png",
       p_weights, width = 10, height = 5, dpi = 300, bg = "white")
ggsave("outputs/figures/efigure_ipw_weights.pdf",
       p_weights, width = 10, height = 5, bg = "white")
cat("Saved: outputs/figures/efigure_ipw_weights.png\n\n")

# ============================================================================
# 2. XGBOOST PROPENSITY SENSITIVITY 
# ============================================================================

cat("--- 2. XGBoost Propensity Sensitivity ---\n\n")

encode_for_xgb <- function(df, cols) {
  num_cols <- cols[sapply(df[cols], is.numeric)]
  cat_cols  <- cols[sapply(df[cols], function(x) is.character(x) | is.factor(x))]
  X_num <- as.matrix(df[, num_cols, drop = FALSE])
  X_num[is.na(X_num)] <- 0
  if (length(cat_cols) > 0) {
    dummies <- model.matrix(
      as.formula(paste("~", paste(cat_cols, collapse = " + "), "- 1")),
      data = df[, cat_cols, drop = FALSE]
    )
    return(cbind(X_num, dummies))
  }
  return(X_num)
}

combined_cols <- c(patient_features,
                   "therapist_name", "location", "program")
combined_cols <- combined_cols[combined_cols %in% names(model_df)]

xgb_prop_results <- list()

for (col in names(therapy_labels)) {
  if (!col %in% names(model_df)) next
  
  y_train_t <- ifelse(is.na(train_df[[col]]), 0, train_df[[col]])
  y_valid_t <- ifelse(is.na(valid_df[[col]]), 0, valid_df[[col]])
  if (length(unique(y_train_t)) < 2) next
  
  X_tr  <- encode_for_xgb(train_df, combined_cols)
  X_val <- encode_for_xgb(valid_df, combined_cols)
  
  # Align columns
  missing_val <- setdiff(colnames(X_tr), colnames(X_val))
  if (length(missing_val) > 0) {
    z <- matrix(0, nrow = nrow(X_val), ncol = length(missing_val))
    colnames(z) <- missing_val
    X_val <- cbind(X_val, z)
  }
  X_val <- X_val[, colnames(X_tr), drop = FALSE]
  
  class_counts <- table(y_train_t)
  w <- ifelse(y_train_t == 1,
              length(y_train_t) / (2 * class_counts["1"]),
              length(y_train_t) / (2 * class_counts["0"]))
  
  xgb_mod <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc",
                  max_depth = 4, eta = 0.05, subsample = 0.8,
                  colsample_bytree = 0.8, min_child_weight = 5),
    data    = xgb.DMatrix(X_tr, label = y_train_t, weight = w),
    nrounds = 200, verbose = 0
  )
  
  pred_val   <- predict(xgb_mod, xgb.DMatrix(X_val))
  auc_xgb    <- as.numeric(auc(roc(y_valid_t, pred_val, quiet = TRUE)))
  ridge_auc  <- summary_stats$Combined_AUC[summary_stats$Therapy == therapy_labels[col]]
  
  cat(sprintf("  %s: Ridge AUC = %.3f | XGBoost AUC = %.3f\n",
              therapy_labels[col], ridge_auc, auc_xgb))
  
  xgb_prop_results[[col]] <- list(
    label     = therapy_labels[col],
    auc_ridge = ridge_auc,
    auc_xgb   = auc_xgb,
    ridge_adv = ridge_auc - auc_xgb
  )
}

propensity_comparison <- do.call(rbind, lapply(xgb_prop_results, function(x) {
  data.frame(Modality = x$label, Ridge_Combined = round(x$auc_ridge, 3),
             XGBoost_Combined = round(x$auc_xgb, 3),
             Ridge_Advantage  = round(x$ridge_adv, 3))
}))

# Add patient and therapist AUCs from summary_stats
propensity_comparison <- propensity_comparison %>%
  left_join(
    summary_stats %>%
      select(Therapy, Patient_AUC, Therapist_AUC, DeLong_p) %>%
      rename(Modality = Therapy,
             Patient_Model   = Patient_AUC,
             Therapist_Model = Therapist_AUC,
             p_value         = DeLong_p) %>%
      mutate(across(c(Patient_Model, Therapist_Model), ~round(.x, 3))),
    by = "Modality"
  ) %>%
  select(Modality, Patient_Model, Therapist_Model,
         Ridge_Combined, XGBoost_Combined, Ridge_Advantage, p_value) %>%
  arrange(desc(Ridge_Combined - Patient_Model))

cat("\nPropensity Model Comparison Table:\n")
print(propensity_comparison, row.names = FALSE)
write.csv(propensity_comparison,
          "outputs/tables/sensitivity_propensity_comparison.csv",
          row.names = FALSE)
cat("Saved: outputs/tables/sensitivity_propensity_comparison.csv\n\n")

# ============================================================================
# 3. E-VALUES 
# ============================================================================

cat("--- 3. E-values for Unmeasured Confounding ---\n\n")

baseline_rate <- mean(y_train, na.rm = TRUE)
cat(sprintf("Baseline improvement rate (training): %.3f\n\n", baseline_rate))

evalue_table <- data.frame()

for (col in names(therapy_labels)) {
  if (!col %in% names(dr_results)) next
  res    <- dr_results[[col]]
  ate    <- res$ate
  ci_low <- res$ci_lower
  label  <- therapy_labels[col]
  
  p0     <- baseline_rate
  RR     <- (p0 + ate) / p0
  RR_ci  <- (p0 + ci_low) / p0
  
  evalue    <- ifelse(RR > 1, RR + sqrt(RR * (RR - 1)), NA)
  evalue_ci <- ifelse(RR_ci > 1, RR_ci + sqrt(RR_ci * (RR_ci - 1)), 1)
  
  cat(sprintf("  %-20s  ATE = %.1fpp  CI lower = %.1fpp  E-value (CI) = %.2f\n",
              label, ate * 100, ci_low * 100, evalue_ci))
  
  evalue_table <- rbind(evalue_table, data.frame(
    Modality    = label,
    ATE_pp      = round(ate * 100, 1),
    CI_lower_pp = round(ci_low * 100, 1),
    p_value     = round(res$p_value, 3),
    E_value_CI  = round(evalue_ci, 2)
  ))
}

evalue_table <- evalue_table %>% arrange(desc(ATE_pp))
cat("\nE-value Table:\n")
print(evalue_table, row.names = FALSE)
write.csv(evalue_table,
          "outputs/tables/sensitivity_evalues.csv",
          row.names = FALSE)

# E-value forest plot
evalue_plot_data <- evalue_table %>%
  mutate(Modality = factor(Modality, levels = Modality[order(ATE_pp)]))

p_eval <- ggplot(evalue_plot_data, aes(y = Modality)) +
  geom_segment(aes(x = CI_lower_pp, xend = ATE_pp, y = Modality, yend = Modality),
               color = "#1D9E75", linewidth = 0.8) +
  geom_point(aes(x = ATE_pp), color = "#1D9E75", size = 3) +
  geom_point(aes(x = CI_lower_pp), color = "#1D9E75", size = 2, shape = 124) +
  geom_text(aes(x = ATE_pp + 0.8, label = sprintf("E = %.2f", E_value_CI)),
            hjust = 0, size = 3.2, color = "#5F5E5A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_x_continuous(name = "Average Treatment Effect (percentage points)",
                     breaks = seq(0, 16, 2), limits = c(-1, 17)) +
  labs(y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(size = 11.5, face = "bold", hjust = 0,
                                      margin = margin(b = 8)),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9, color = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.2),
    plot.margin        = margin(8, 8, 4, 8)
  )

ggsave("outputs/figures/efigure_evalues.png",
       p_eval, width = 7, height = 4, dpi = 300, bg = "white")
ggsave("outputs/figures/efigure_evalues.pdf",
       p_eval, width = 7, height = 4, bg = "white")
cat("Saved: outputs/figures/efigure_evalues.png\n\n")

# ============================================================================
# 4. MONTHLY CALIBRATION
# ============================================================================

cat("--- 4. Monthly Calibration and Drift ---\n\n")

test_data_cal <- test_data %>%
  mutate(pred_prob  = pred_test,
         observed   = y_test,
         month_year = floor_date(as.Date(admission_date), "month"))

monthly_perf <- test_data_cal %>%
  group_by(month_year) %>%
  summarise(
    n         = n(),
    obs_rate  = round(100 * mean(observed), 1),
    mean_pred = round(100 * mean(pred_prob), 1),
    brier     = round(mean((pred_prob - observed)^2), 4),
    auroc     = ifelse(n >= 20 & length(unique(observed)) > 1,
                       round(as.numeric(auc(roc(observed, pred_prob, quiet = TRUE))), 3),
                       NA),
    .groups = "drop"
  ) %>%
  arrange(month_year)

cat("Monthly performance:\n")
print(monthly_perf)

# Calibration intercept and slope
cal_intercept_model <- glm(observed ~ offset(qlogis(pred_prob)),
                           data = test_data_cal, family = binomial)
cal_intercept <- round(coef(cal_intercept_model)[1], 3)

cal_slope_model <- glm(observed ~ qlogis(pred_prob),
                       data = test_data_cal, family = binomial)
cal_slope <- round(coef(cal_slope_model)[2], 3)

cat(sprintf("\nCalibration intercept: %.3f\n", cal_intercept))
cat(sprintf("Calibration slope:     %.3f\n", cal_slope))
cat(sprintf("Mean predicted:        %.3f\n", mean(test_data_cal$pred_prob)))
cat(sprintf("Observed rate:         %.3f\n", mean(test_data_cal$observed)))

# Decile calibration
cal_deciles <- test_data_cal %>%
  mutate(decile = ntile(pred_prob, 10)) %>%
  group_by(decile) %>%
  summarise(mean_pred = mean(pred_prob),
            obs_rate  = mean(observed),
            n         = n(), .groups = "drop")

cat("\nDecile calibration:\n")
print(cal_deciles)

write.csv(monthly_perf,
          "outputs/tables/sensitivity_monthly_calibration.csv",
          row.names = FALSE)
write.csv(cal_deciles,
          "outputs/tables/sensitivity_decile_calibration.csv",
          row.names = FALSE)
cat("Saved calibration tables.\n\n")

# ============================================================================
# 5. SINGLE-MODALITY SENSITIVITY (R3 M2)
# ============================================================================

cat("--- 5. Single-Modality Sensitivity ---\n\n")

data <- data %>%
  mutate(n_modalities = rowSums(select(., all_of(therapy_cols)), na.rm = TRUE))

cat("Modality count distribution:\n")
data %>%
  count(n_modalities) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

single_mod_test <- test_data %>%
  mutate(n_modalities = rowSums(select(., all_of(therapy_cols)), na.rm = TRUE)) %>%
  filter(n_modalities == 1)

idx_single <- which(test_data$row_id %in% single_mod_test$row_id)
# If row_id not available, use row position matching
idx_single <- which(rowSums(
  test_data %>% select(all_of(therapy_cols)), na.rm = TRUE) == 1)

cat(sprintf("\nSingle-modality test episodes: %d (%.1f%% of test set)\n",
            length(idx_single),
            100 * length(idx_single) / nrow(test_data)))

if (length(idx_single) >= 30) {
  gains_single  <- personalization_gains[idx_single]
  y_single      <- y_test[idx_single]
  pred_single   <- pred_test[idx_single]
  
  auroc_single  <- as.numeric(auc(roc(y_single, pred_single, quiet = TRUE)))
  
  single_summary <- data.frame(
    Metric        = c("N (test set)", "Observed improvement rate", "AUROC",
                      "Mean gain (pp)", "% benefiting (>1pp)", "NNT"),
    Full_Sample   = c(length(y_test),
                      round(100 * mean(y_test), 1),
                      round(auc_test, 3),
                      round(mean(personalization_gains) * 100, 1),
                      round(100 * mean(personalization_gains > 0.01), 1),
                      round(1 / mean(personalization_gains))),
    Single_Mod    = c(length(y_single),
                      round(100 * mean(y_single), 1),
                      round(auroc_single, 3),
                      round(mean(gains_single) * 100, 1),
                      round(100 * mean(gains_single > 0.01), 1),
                      round(1 / mean(gains_single)))
  )
  
  cat("\nFull sample vs. single-modality comparison:\n")
  print(single_summary, row.names = FALSE)
  write.csv(single_summary,
            "outputs/tables/sensitivity_single_modality.csv",
            row.names = FALSE)
}
cat("\n")

# ============================================================================
# 6. DEMOGRAPHIC SUBGROUP GAINS
# ============================================================================

cat("--- 6. Demographic Subgroup Personalization Gains ---\n\n")

subgroups <- list(
  list(name = "Adolescent", var = "adolescent", val = 1),
  list(name = "Adult",      var = "adolescent", val = 0),
  list(name = "Female",     var = "male",       val = 0),
  list(name = "Male",       var = "male",       val = 1),
  list(name = "RTC",        var = "program",    val = "RTC"),
  list(name = "PHP",        var = "program",    val = "PHP"),
  list(name = "OP",         var = "program",    val = "OP")
)

subgroup_results <- data.frame()

for (sg in subgroups) {
  idx <- which(test_data[[sg$var]] == sg$val)
  if (length(idx) < 30) next
  
  gains_sg <- personalization_gains[idx]
  y_sg     <- y_test[idx]
  pred_sg  <- pred_test[idx]
  
  auroc_sg <- tryCatch(
    as.numeric(auc(roc(y_sg, pred_sg, quiet = TRUE))),
    error = function(e) NA)
  
  subgroup_results <- rbind(subgroup_results, data.frame(
    Subgroup     = sg$name,
    N            = length(idx),
    Obs_rate     = round(100 * mean(y_sg), 1),
    AUROC        = round(auroc_sg, 3),
    Mean_gain_pp = round(mean(gains_sg) * 100, 1),
    Pct_over_1pp = round(100 * mean(gains_sg > 0.01), 1),
    NNT          = round(1 / mean(gains_sg))
  ))
}

cat("Demographic subgroup results:\n")
print(subgroup_results, row.names = FALSE)
write.csv(subgroup_results,
          "outputs/tables/sensitivity_subgroup_gains.csv",
          row.names = FALSE)
cat("Saved: outputs/tables/sensitivity_subgroup_gains.csv\n\n")

# ============================================================================
# 7. HTE / CLINICAL SUBGROUP GAINS
# ============================================================================

cat("--- 7. HTE by Clinical Subgroup ---\n\n")

# 7a. By initial risk level
cat("By initial risk level:\n")
hte_risk <- data.frame(
  risk  = test_data$risk_level_initial,
  gain  = personalization_gains,
  obs   = y_test
) %>%
  group_by(risk) %>%
  summarise(n = n(),
            mean_gain_pp = round(mean(gain) * 100, 1),
            pct_over_1pp = round(100 * mean(gain > 0.01), 1),
            NNT          = round(1 / mean(gain)),
            .groups = "drop")
print(hte_risk)

# 7b. By C-SSRS quartile
cat("\nBy baseline C-SSRS quartile:\n")
hte_cssrs <- data.frame(
  cssrs = test_data$total_score,
  gain  = personalization_gains,
  obs   = y_test
) %>%
  mutate(q = ntile(cssrs, 4)) %>%
  group_by(q) %>%
  summarise(n            = n(),
            mean_cssrs   = round(mean(cssrs), 1),
            mean_gain_pp = round(mean(gain) * 100, 1),
            pct_over_1pp = round(100 * mean(gain > 0.01), 1),
            NNT          = round(1 / mean(gain)),
            .groups = "drop")
print(hte_cssrs)

# 7c. By predicted probability quartile
cat("\nBy baseline predicted probability quartile:\n")
hte_prob <- data.frame(
  pred = pred_test,
  gain = personalization_gains,
  obs  = y_test
) %>%
  mutate(q = ntile(pred, 4)) %>%
  group_by(q) %>%
  summarise(n            = n(),
            mean_pred    = round(mean(pred) * 100, 1),
            mean_gain_pp = round(mean(gain) * 100, 1),
            pct_over_1pp = round(100 * mean(gain > 0.01), 1),
            NNT          = round(1 / mean(gain)),
            .groups = "drop")
print(hte_prob)

# 7d. By primary diagnosis
cat("\nBy primary diagnosis (n >= 30):\n")
hte_dx <- data.frame(
  dx   = test_data$dx_group,
  gain = personalization_gains,
  obs  = y_test
) %>%
  group_by(dx) %>%
  summarise(n            = n(),
            mean_gain_pp = round(mean(gain) * 100, 1),
            pct_over_1pp = round(100 * mean(gain > 0.01), 1),
            NNT          = round(1 / mean(gain)),
            .groups = "drop") %>%
  filter(n >= 30) %>%
  arrange(desc(mean_gain_pp))
print(hte_dx)

write.csv(hte_risk,  "outputs/tables/sensitivity_hte_risk.csv",  row.names = FALSE)
write.csv(hte_cssrs, "outputs/tables/sensitivity_hte_cssrs.csv", row.names = FALSE)
write.csv(hte_prob,  "outputs/tables/sensitivity_hte_prob.csv",  row.names = FALSE)
write.csv(hte_dx,    "outputs/tables/sensitivity_hte_dx.csv",    row.names = FALSE)
cat("Saved HTE tables.\n\n")

# ============================================================================
# 8. TOP 30 COMBINATIONS SENSITIVITY
# ============================================================================

cat("--- 8. Top 30 Combinations Sensitivity ---\n\n")

combo_strings  <- apply(X_test[, therapy_cols], 1, paste, collapse = "-")
combo_counts   <- table(combo_strings)
top15_combos   <- names(sort(combo_counts, decreasing = TRUE)[1:15])
top30_combos   <- names(sort(combo_counts, decreasing = TRUE)[1:min(30, length(combo_counts))])

pct_covered_15 <- round(100 * mean(combo_strings %in% top15_combos), 1)
pct_covered_30 <- round(100 * mean(combo_strings %in% top30_combos), 1)

cat(sprintf("Coverage — top 15: %.1f%% | top 30: %.1f%%\n\n",
            pct_covered_15, pct_covered_30))

cat("Running top 30 personalization analysis...\n")
gains_top30 <- numeric(nrow(X_test))

pb <- txtProgressBar(min = 0, max = nrow(X_test), style = 3)
for (i in 1:nrow(X_test)) {
  pf       <- X_test[i, ]
  best_p   <- pred_test[i]
  for (combo_str in top30_combos) {
    combo_vals <- as.numeric(strsplit(combo_str, "-")[[1]])
    cf         <- pf
    cf[therapy_cols] <- combo_vals
    cf_p <- predict(xgb_final, xgb.DMatrix(matrix(cf, nrow = 1)))
    if (cf_p > best_p) best_p <- cf_p
  }
  gains_top30[i] <- best_p - pred_test[i]
  setTxtProgressBar(pb, i)
}
close(pb)

top30_summary <- data.frame(
  Metric = c("Coverage of validation episodes",
             "Mean personalization gain",
             "Median personalization gain",
             "% with gain > 1pp",
             "% with gain > 5pp",
             "NNT",
             "Maximum gain"),
  Top_15 = c(paste0(pct_covered_15, "%"),
             paste0(round(mean(personalization_gains) * 100, 1), " pp"),
             paste0(round(median(personalization_gains) * 100, 1), " pp"),
             paste0(round(100 * mean(personalization_gains > 0.01), 1), "%"),
             paste0(round(100 * mean(personalization_gains > 0.05), 1), "%"),
             as.character(round(1 / mean(personalization_gains))),
             paste0(round(max(personalization_gains) * 100, 1), " pp")),
  Top_30 = c(paste0(pct_covered_30, "%"),
             paste0(round(mean(gains_top30) * 100, 1), " pp"),
             paste0(round(median(gains_top30) * 100, 1), " pp"),
             paste0(round(100 * mean(gains_top30 > 0.01), 1), "%"),
             paste0(round(100 * mean(gains_top30 > 0.05), 1), "%"),
             as.character(round(1 / mean(gains_top30))),
             paste0(round(max(gains_top30) * 100, 1), " pp"))
)

cat("\nTop 15 vs Top 30 comparison:\n")
print(top30_summary, row.names = FALSE)
write.csv(top30_summary,
          "outputs/tables/sensitivity_top30_combinations.csv",
          row.names = FALSE)
cat("Saved: outputs/tables/sensitivity_top30_combinations.csv\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("================================================================\n")
cat("SENSITIVITY ANALYSES COMPLETE\n")
cat("================================================================\n\n")

cat("Tables saved to outputs/tables/:\n")
cat("  sensitivity_ipw_weights.csv\n")
cat("  sensitivity_propensity_comparison.csv\n")
cat("  sensitivity_evalues.csv\n")
cat("  sensitivity_monthly_calibration.csv\n")
cat("  sensitivity_decile_calibration.csv\n")
cat("  sensitivity_single_modality.csv\n")
cat("  sensitivity_subgroup_gains.csv\n")
cat("  sensitivity_hte_risk.csv\n")
cat("  sensitivity_hte_cssrs.csv\n")
cat("  sensitivity_hte_prob.csv\n")
cat("  sensitivity_hte_dx.csv\n")
cat("  sensitivity_top30_combinations.csv\n\n")

cat("Figures saved to outputs/figures/:\n")
cat("  efigure_ipw_weights.png/.pdf\n")
cat("  efigure_evalues.png/.pdf\n")

