# ============================================================================
# run_all.R
# Master script to reproduce all analyses and outputs
# ============================================================================

rm(list = ls())
gc()

cat("\n==============================\n")
cat("RUNNING FULL PIPELINE\n")
cat("==============================\n\n")

# ----------------------------------------------------------------------------
# Load required packages
# ----------------------------------------------------------------------------
packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "lubridate",
  "glmnet",
  "xgboost",
  "pROC",
  "MatchIt",
  "cowplot",
  "gridExtra",
  "ComplexHeatmap",
  "grid"
)

invisible(lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

# ----------------------------------------------------------------------------
# Run scripts in order
# ----------------------------------------------------------------------------
scripts <- c(
  "src/01_clean.R",
  "src/02_propensity_scores.R",
  "src/03_outcome_model.R",
  "src/04_doubly_robust.R",
  "src/05_counterfactual.R",
  "src/06_matching.R",
  "src/07_figures.R",
  "src/08_tables.R",
  "src/09_sensitivity.R"
)

for (s in scripts) {
  cat("\n---------------------------------\n")
  cat("Running:", s, "\n")
  cat("---------------------------------\n")
  source(s)
}

cat("\n==============================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================\n")