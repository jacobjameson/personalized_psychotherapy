# ============================================================================
# 02_propensity_scores.R
# Propensity Score Estimation for Psychotherapy Modality Assignment
# ============================================================================
#
# This script estimates the probability of receiving each of seven
# psychotherapy modalities using L2-regularized (ridge) logistic regression.
# Three models are trained per modality:
#   1. Patient clinical features only
#   2. Therapist/organizational features only
#   3. Combined (all features)
#
# Key outputs:
#   - Propensity scores (from combined model) merged into `data`
#   - ROC curves and DeLong tests comparing model types
#   - Table 2 summary statistics
#
# Depends on: 01_clean.R (provides `data`)
# ============================================================================

set.seed(42)

# ============================================================================
# 1. DEFINE FEATURE GROUPS
# ============================================================================

# Therapy modality indicators
therapy_cols <- c("act", "cbt", "dbt", "motivational_interviewing",
                  "mindfulness", "stages_of_change", "family_systems",
                  "trauma_informed")

# Patient clinical features (extracted from structured intake fields)
dx_cols        <- names(data)[grepl("^current_and_past_psychiatric_diagnoses_", names(data))]
symp_cols      <- names(data)[grepl("^presenting_symptoms_", names(data))]
fam_cols       <- names(data)[grepl("^family_history_", names(data))]
precip_cols    <- names(data)[grepl("^precipitants_stressors_", names(data))]
prot_cols      <- names(data)[grepl("^internal_protective_factors_|^external_protective_factors_", names(data))]
treatment_cols <- names(data)[grepl("^change_in_treatment_", names(data))]

# Derived variables
if ("risk_level_initial" %in% names(data)) {
  data$risk_high_initial <- as.integer(data$risk_level_initial == "High")
}
data$adolescent <- as.integer(data$age_group == "Adolescent")
data$male       <- as.integer(data$sex_fs == "male")

# All patient features for predicting therapy assignment
patient_features <- unique(c(
  "total_score", "risk_high_initial", "adolescent", "male",
  dx_cols, symp_cols, fam_cols, precip_cols, prot_cols, treatment_cols
))

# Therapist/organizational features
therapist_feature_cols <- c("therapist_name", "location", "program")

cat("Feature counts:\n")
cat(sprintf("  Diagnosis features: %d\n", length(dx_cols)))
cat(sprintf("  Symptom features: %d\n", length(symp_cols)))
cat(sprintf("  Family history features: %d\n", length(fam_cols)))
cat(sprintf("  Precipitant/stressor features: %d\n", length(precip_cols)))
cat(sprintf("  Protective factor features: %d\n", length(prot_cols)))
cat(sprintf("  Total patient features: %d\n", length(patient_features)))

# ============================================================================
# 2. CREATE MODELING DATASET AND TEMPORAL SPLIT
# ============================================================================

# Metadata columns needed for modeling
metadata_cols <- c("therapist_name", "location", "program",
                   "pn_month", "pn_time_block", "pn_year", "admission_date")

# Assemble modeling features
model_features <- unique(c(
  "risk_high_initial", "total_score", "adolescent", "male",
  therapy_cols, dx_cols, symp_cols, fam_cols, precip_cols, prot_cols, treatment_cols
))

model_df <- data[, unique(c(model_features, metadata_cols))]
model_df <- model_df[, !duplicated(names(model_df))]
model_df$row_idx <- 1:nrow(model_df)

# Filter patient features to those present in model_df
patient_features <- patient_features[patient_features %in% names(model_df)]

# Temporal split: earliest 80% for training, latest 20% for validation
model_df <- model_df %>%
  arrange(admission_date) %>%
  mutate(split_row_id = row_number())

split_idx <- floor(nrow(model_df) * 0.8)
train_df  <- model_df %>% filter(split_row_id <= split_idx)
valid_df  <- model_df %>% filter(split_row_id > split_idx)

cat(sprintf("\nTemporal split:\n"))
cat(sprintf("  Training set: %d observations (80%%)\n", nrow(train_df)))
cat(sprintf("  Validation set: %d observations (20%%)\n", nrow(valid_df)))

# ============================================================================
# 3. MODEL TRAINING FUNCTION
# ============================================================================
# For each therapy modality, trains three L2-regularized logistic regression
# models (patient-only, therapist-only, combined), evaluates discrimination
# via AUROC in the validation set, and generates propensity scores for the
# full dataset.
# ============================================================================

get_curves_for_therapy <- function(train_df, valid_df, full_df, therapy_col,
                                   label, patient_features,
                                   therapist_cols = c("therapist_name", "location", "program")) {
  
  # Extract outcome vectors
  y_train <- ifelse(is.na(train_df[[therapy_col]]), 0, train_df[[therapy_col]])
  y_test  <- ifelse(is.na(valid_df[[therapy_col]]), 0, valid_df[[therapy_col]])
  y_full  <- ifelse(is.na(full_df[[therapy_col]]),  0, full_df[[therapy_col]])
  
  # Require sufficient variance and class balance
  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) {
    cat(sprintf("  Skipping %s: insufficient outcome variance\n", label))
    return(NULL)
  }
  if (sum(y_test == 1) < 5 || sum(y_test == 0) < 5) {
    cat(sprintf("  Skipping %s: insufficient test set size\n", label))
    return(NULL)
  }
  
  available_therapist_cols <- therapist_cols[therapist_cols %in% names(train_df)]
  
  # --- Helper: prepare design matrices for glmnet ---
  prepare_data_for_glmnet <- function(train_data, test_data, feature_cols,
                                      include_categorical = FALSE) {
    num_cols <- sapply(train_data[, feature_cols, drop = FALSE], is.numeric)
    numeric_features    <- feature_cols[num_cols]
    categorical_features <- feature_cols[!num_cols]
    
    # Numeric features: impute missing with training median
    if (length(numeric_features) > 0) {
      X_train_num <- as.matrix(train_data[, numeric_features, drop = FALSE])
      X_test_num  <- as.matrix(test_data[, numeric_features, drop = FALSE])
      for (i in 1:ncol(X_train_num)) {
        median_val <- median(X_train_num[, i], na.rm = TRUE)
        if (is.na(median_val)) median_val <- 0
        X_train_num[is.na(X_train_num[, i]), i] <- median_val
        X_test_num[is.na(X_test_num[, i]), i]   <- median_val
      }
    } else {
      X_train_num <- matrix(nrow = nrow(train_data), ncol = 0)
      X_test_num  <- matrix(nrow = nrow(test_data),  ncol = 0)
    }
    
    # Categorical features: one-hot encode, align columns
    if (include_categorical && length(categorical_features) > 0) {
      formula_obj <- as.formula(paste("~", paste(categorical_features, collapse = " + "), "- 1"))
      X_train_cat <- model.matrix(formula_obj, data = train_data[, categorical_features, drop = FALSE])
      X_test_cat  <- model.matrix(formula_obj, data = test_data[, categorical_features, drop = FALSE])
      
      # Align columns between train and test
      missing_cols <- setdiff(colnames(X_train_cat), colnames(X_test_cat))
      if (length(missing_cols) > 0) {
        zeros <- matrix(0, nrow = nrow(X_test_cat), ncol = length(missing_cols))
        colnames(zeros) <- missing_cols
        X_test_cat <- cbind(X_test_cat, zeros)
      }
      extra_cols <- setdiff(colnames(X_test_cat), colnames(X_train_cat))
      if (length(extra_cols) > 0) {
        X_test_cat <- X_test_cat[, !colnames(X_test_cat) %in% extra_cols, drop = FALSE]
      }
      X_test_cat <- X_test_cat[, colnames(X_train_cat), drop = FALSE]
    } else {
      X_train_cat <- matrix(nrow = nrow(train_data), ncol = 0)
      X_test_cat  <- matrix(nrow = nrow(test_data),  ncol = 0)
    }
    
    list(train = cbind(X_train_num, X_train_cat),
         test  = cbind(X_test_num, X_test_cat))
  }
  
  tryCatch({
    # Prepare design matrices for validation and full dataset
    patient_data   <- prepare_data_for_glmnet(train_df, valid_df, patient_features)
    therapist_data <- prepare_data_for_glmnet(train_df, valid_df, available_therapist_cols,
                                              include_categorical = TRUE)
    combined_data  <- prepare_data_for_glmnet(train_df, valid_df,
                                              c(patient_features, available_therapist_cols),
                                              include_categorical = TRUE)
    
    patient_data_full   <- prepare_data_for_glmnet(train_df, full_df, patient_features)
    therapist_data_full <- prepare_data_for_glmnet(train_df, full_df, available_therapist_cols,
                                                   include_categorical = TRUE)
    combined_data_full  <- prepare_data_for_glmnet(train_df, full_df,
                                                   c(patient_features, available_therapist_cols),
                                                   include_categorical = TRUE)
    
    if (ncol(therapist_data$train) == 0) {
      cat(sprintf("  Warning: No therapist features for %s\n", label))
      return(NULL)
    }
    
    # Class weights to address imbalanced therapy frequencies
    class_counts <- table(y_train)
    if (length(class_counts) < 2) return(NULL)
    weights <- ifelse(y_train == 1,
                      length(y_train) / (2 * class_counts["1"]),
                      length(y_train) / (2 * class_counts["0"]))
    
    # Train L2-regularized logistic regression (ridge) with 5-fold CV
    cv_patient   <- cv.glmnet(x = patient_data$train,   y = y_train, family = "binomial",
                              alpha = 0, weights = weights, type.measure = "auc", nfolds = 5)
    cv_therapist <- cv.glmnet(x = therapist_data$train, y = y_train, family = "binomial",
                              alpha = 0, weights = weights, type.measure = "auc", nfolds = 5)
    cv_combined  <- cv.glmnet(x = combined_data$train,  y = y_train, family = "binomial",
                              alpha = 0, weights = weights, type.measure = "auc", nfolds = 5)
    
    # Predictions for validation set (evaluation)
    y_pred_patient   <- predict(cv_patient,   newx = patient_data$test,   s = "lambda.min", type = "response")[, 1]
    y_pred_therapist <- predict(cv_therapist, newx = therapist_data$test, s = "lambda.min", type = "response")[, 1]
    y_pred_combined  <- predict(cv_combined,  newx = combined_data$test,  s = "lambda.min", type = "response")[, 1]
    
    # Propensity scores for full dataset (used downstream in AIPW)
    prop_patient_full   <- predict(cv_patient,   newx = patient_data_full$test,   s = "lambda.min", type = "response")[, 1]
    prop_therapist_full <- predict(cv_therapist, newx = therapist_data_full$test, s = "lambda.min", type = "response")[, 1]
    prop_combined_full  <- predict(cv_combined,  newx = combined_data_full$test,  s = "lambda.min", type = "response")[, 1]
    
    # ROC curves and AUC
    roc_patient   <- roc(y_test, y_pred_patient,   quiet = TRUE, direction = "<")
    roc_therapist <- roc(y_test, y_pred_therapist, quiet = TRUE, direction = "<")
    roc_combined  <- roc(y_test, y_pred_combined,  quiet = TRUE, direction = "<")
    
    auc_p <- as.numeric(auc(roc_patient))
    auc_t <- as.numeric(auc(roc_therapist))
    auc_c <- as.numeric(auc(roc_combined))
    
    # DeLong tests comparing correlated AUCs
    delong_t_vs_p <- tryCatch({
      res <- roc.test(roc_therapist, roc_patient, method = "delong")
      list(statistic = res$statistic, p_value = res$p.value)
    }, error = function(e) list(statistic = NA, p_value = NA))
    
    delong_c_vs_p <- tryCatch({
      res <- roc.test(roc_combined, roc_patient, method = "delong")
      list(statistic = res$statistic, p_value = res$p.value)
    }, error = function(e) list(statistic = NA, p_value = NA))
    
    # Brier scores
    brier_patient   <- mean((y_test - y_pred_patient)^2)
    brier_therapist <- mean((y_test - y_pred_therapist)^2)
    brier_combined  <- mean((y_test - y_pred_combined)^2)
    
    cat(sprintf("  %s: Patient AUC=%.3f, Therapist AUC=%.3f, Combined AUC=%.3f (n=%d/%d, p=%.4f)\n",
                label, auc_p, auc_t, auc_c, sum(y_test == 1), length(y_test),
                delong_t_vs_p$p_value))
    
    list(
      therapy       = label,
      baseline_rate = mean(y_test),
      n_positive    = sum(y_test == 1),
      n_total       = length(y_test),
      patient   = list(fpr = 1 - roc_patient$specificities,
                       tpr = roc_patient$sensitivities, auc = auc_p),
      therapist = list(fpr = 1 - roc_therapist$specificities,
                       tpr = roc_therapist$sensitivities, auc = auc_t),
      combined  = list(fpr = 1 - roc_combined$specificities,
                       tpr = roc_combined$sensitivities, auc = auc_c),
      delong_tests = list(t_vs_p = delong_t_vs_p, c_vs_p = delong_c_vs_p),
      brier = list(patient = brier_patient, therapist = brier_therapist,
                   combined = brier_combined),
      propensity_scores = list(patient = prop_patient_full,
                               therapist = prop_therapist_full,
                               combined = prop_combined_full),
      models = list(patient = cv_patient, therapist = cv_therapist,
                    combined = cv_combined)
    )
    
  }, error = function(e) {
    cat(sprintf("  Error processing %s: %s\n", label, e$message))
    return(NULL)
  })
}

# ============================================================================
# 4. RUN MODELS FOR ALL THERAPY MODALITIES
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PROPENSITY SCORE ESTIMATION\n")
cat(rep("=", 70), "\n", sep = "")

therapy_labels <- list(
  "act"                       = "ACT",
  "cbt"                       = "CBT",
  "dbt"                       = "DBT",
  "motivational_interviewing" = "MI",
  "mindfulness"               = "Mindfulness",
  "stages_of_change"          = "Stages of Change",
  "family_systems"            = "Family Systems"
)

results <- list()

for (col in names(therapy_labels)) {
  label <- therapy_labels[[col]]
  if (col %in% names(model_df)) {
    res <- get_curves_for_therapy(
      train_df, valid_df, model_df, col, label,
      patient_features,
      therapist_cols = therapist_feature_cols
    )
    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  } else {
    cat(sprintf("  Warning: Therapy column '%s' not found\n", col))
  }
}

# ============================================================================
# 5. MERGE PROPENSITY SCORES INTO DATA
# ============================================================================

propensity_df <- data.frame(row_idx = model_df$row_idx)

for (i in seq_along(results)) {
  therapy_col_clean <- tolower(gsub(" ", "_", results[[i]]$therapy))
  propensity_df[[paste0("prop_", therapy_col_clean)]] <-
    results[[i]]$propensity_scores$combined
}

data <- data %>%
  mutate(row_idx = row_number()) %>%
  left_join(propensity_df, by = "row_idx") %>%
  select(-row_idx)

cat(sprintf("\nPropensity scores added: %d columns\n", length(results)))

# ============================================================================
# 6. SUMMARY STATISTICS AND TABLE 2
# ============================================================================

summary_stats <- data.frame(
  Therapy       = sapply(results, function(x) x$therapy),
  Baseline_Rate = sapply(results, function(x) x$baseline_rate),
  N_Positive    = sapply(results, function(x) x$n_positive),
  N_Total       = sapply(results, function(x) x$n_total),
  Patient_AUC   = sapply(results, function(x) x$patient$auc),
  Therapist_AUC = sapply(results, function(x) x$therapist$auc),
  Combined_AUC  = sapply(results, function(x) x$combined$auc),
  DeLong_Z      = sapply(results, function(x) x$delong_tests$t_vs_p$statistic),
  DeLong_p      = sapply(results, function(x) x$delong_tests$t_vs_p$p_value),
  Brier_Patient   = sapply(results, function(x) x$brier$patient),
  Brier_Therapist = sapply(results, function(x) x$brier$therapist),
  Brier_Combined  = sapply(results, function(x) x$brier$combined)
) %>%
  mutate(
    Delta_AUC_T_vs_P = Therapist_AUC - Patient_AUC,
    Delta_AUC_C_vs_P = Combined_AUC  - Patient_AUC,
    Delta_AUC_C_vs_T = Combined_AUC  - Therapist_AUC
  )

# Paired t-test: therapist vs patient AUC
t_test_result <- t.test(summary_stats$Therapist_AUC, summary_stats$Patient_AUC, paired = TRUE)
n_sig <- sum(summary_stats$DeLong_p < 0.05 & summary_stats$Delta_AUC_T_vs_P > 0, na.rm = TRUE)

# --- Print to console ---
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("RESULTS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

cat(sprintf("\nMean AUC — Patient: %.3f (SD=%.3f), Therapist: %.3f (SD=%.3f), Combined: %.3f (SD=%.3f)\n",
            mean(summary_stats$Patient_AUC),   sd(summary_stats$Patient_AUC),
            mean(summary_stats$Therapist_AUC), sd(summary_stats$Therapist_AUC),
            mean(summary_stats$Combined_AUC),  sd(summary_stats$Combined_AUC)))
cat(sprintf("Mean difference (Therapist - Patient): %.3f (95%% CI: %.3f to %.3f)\n",
            mean(summary_stats$Delta_AUC_T_vs_P),
            quantile(summary_stats$Delta_AUC_T_vs_P, 0.025),
            quantile(summary_stats$Delta_AUC_T_vs_P, 0.975)))
cat(sprintf("Paired t-test: t(%d) = %.2f, p = %.4f\n",
            t_test_result$parameter, t_test_result$statistic, t_test_result$p.value))
cat(sprintf("Significant therapist advantage (p<.05): %d/%d modalities\n", n_sig, nrow(summary_stats)))
cat(sprintf("Adding patient to therapist: mean delta AUC = %.3f\n", mean(summary_stats$Delta_AUC_C_vs_T)))
cat(sprintf("Adding therapist to patient: mean delta AUC = %.3f\n", mean(summary_stats$Delta_AUC_C_vs_P)))

# --- Save Table 2 ---
sink("outputs/tables/table2_propensity.txt")
cat(rep("=", 90), "\n", sep = "")
cat("TABLE 2: THERAPY ASSIGNMENT PROPENSITY MODELS\n")
cat(rep("=", 90), "\n", sep = "")
print(summary_stats %>%
        select(Therapy, Patient_AUC, Therapist_AUC, Combined_AUC,
               Delta_AUC_T_vs_P, Delta_AUC_C_vs_T, DeLong_p) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        arrange(desc(Delta_AUC_T_vs_P)))
cat("\n")
cat(sprintf("Paired t-test (Therapist vs Patient): t(%d) = %.2f, p = %.4f\n",
            t_test_result$parameter, t_test_result$statistic, t_test_result$p.value))
sink()

cat("\nTable 2 saved to: outputs/tables/table2_propensity.txt\n")

# ============================================================================
# 7. SAVE RESULTS OBJECTS
# ============================================================================

# Store results and summary_stats for use by downstream scripts
propensity_results <- results
propensity_summary <- summary_stats

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PROPENSITY SCORE ESTIMATION COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("  Propensity scores for %d modalities merged into `data`\n", length(results)))
cat(sprintf("  Results stored in `propensity_results` and `propensity_summary`\n"))
