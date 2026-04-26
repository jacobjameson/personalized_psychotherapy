# ============================================================================
# 07_figures.R
# Publication Figures for Manuscript
# ============================================================================
#
# Generates all main-text and supplemental figures with a unified theme.
#
# Main figures:
#   Figure 1 — Therapy modality combination UpSet plot
#   Figure 2 — ROC curves for propensity models (7 panels)
#   Figure 3 — SHAP values (A) and doubly robust ATEs (B)
#   Figure 4 — Personalization gains by baseline risk (A) and level (B)
#
#
# Depends on: 01–05 (all analysis objects in environment)
# ============================================================================

# Core color palette
pal <- list(
  red       = "#C1121F",
  red_light = "#E76F51",
  teal      = "#10B981",
  amber     = "#F59E0B",
  indigo    = "#6366F1",
  dark      = "#1a365d",
  grey_bg   = "grey95",
  grey_pt   = "grey80",
  black     = "#333333"
)


# ============================================================================
# FIGURE 1: THERAPY MODALITY COMBINATIONS (UPSET PLOT)
# ============================================================================

cat("Creating Figure 1: UpSet plot...\n")

# Colors (consistent with all figures)
col_red       <- "#C1121F"
col_red_light <- "#E76F51"

therapy_labels_fig <- c(
  "cbt"                       = "CBT",
  "dbt"                       = "DBT",
  "act"                       = "ACT",
  "motivational_interviewing" = "Motivational Interviewing",
  "mindfulness"               = "Mindfulness",
  "stages_of_change"          = "Stages of Change",
  "family_systems"            = "Family Systems"
)

upset_data <- data %>%
  select(all_of(names(therapy_labels_fig))) %>%
  mutate(across(everything(), ~ . == 1))
colnames(upset_data) <- therapy_labels_fig
upset_data <- upset_data[rowSums(upset_data) > 0, ]

m <- make_comb_mat(as.matrix(upset_data))
m <- m[comb_size(m) >= 10]

pdf("outputs/figures/figure1_modality_combinations.pdf", width = 10, height = 4)
UpSet(
  m,
  comb_order = order(-comb_size(m)),
  top_annotation = HeatmapAnnotation(
    "Patients per\nCombination" = anno_barplot(
      comb_size(m),
      border = FALSE,
      gp = gpar(fill = col_red),
      height = unit(4, "cm"),
      add_numbers = FALSE
    ),
    annotation_name_side = "left",
    annotation_name_rot  = 0,
    annotation_name_gp   = gpar(fontsize = 11)
  ),
  right_annotation = rowAnnotation(
    "Set Size" = anno_barplot(
      set_size(m),
      border = FALSE,
      gp = gpar(fill = col_red_light),
      width = unit(3, "cm"),
      add_numbers = TRUE,
      numbers_gp  = gpar(fontsize = 9)
    ),
    annotation_name_side = "bottom",
    annotation_name_rot  = 0,
    annotation_name_gp   = gpar(fontsize = 11)
  ),
  pt_size      = unit(4, "mm"),
  lwd          = 2,
  comb_col     = "black",
  bg_col       = "grey96",
  bg_pt_col    = "grey80",
  row_names_gp = gpar(fontsize = 11),
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

cat("  Saved: outputs/figures/figure1_modality_combinations.pdf\n")

# ============================================================================
# FIGURE 2: THERAPY ASSIGNMENT PROPENSITY MODELS
#   Panel A — ROC curves for seven psychotherapy modalities
#   Panel B — AUROC comparison: patient features vs. therapist factors
# ============================================================================


col_patient   <- "#10B981"
col_therapist <- "#F59E0B"
col_combined  <- "#6366F1"

# ============================================================================
# PANEL A: FACETED ROC CURVES
# ============================================================================

roc_long <- do.call(rbind, lapply(seq_along(propensity_results), function(i) {
  res <- propensity_results[[i]]
  rbind(
    data.frame(therapy = res$therapy, fpr = res$patient$fpr,
               tpr = res$patient$tpr, model = "Patient Features"),
    data.frame(therapy = res$therapy, fpr = res$therapist$fpr,
               tpr = res$therapist$tpr, model = "Therapist Factors"),
    data.frame(therapy = res$therapy, fpr = res$combined$fpr,
               tpr = res$combined$tpr, model = "Combined")
  )
}))

roc_long$model <- factor(roc_long$model,
                         levels = c("Patient Features", "Therapist Factors", "Combined"))

therapy_order <- propensity_summary %>%
  mutate(gap = Therapist_AUC - Patient_AUC) %>%
  arrange(desc(gap)) %>%
  pull(Therapy)

roc_long$therapy <- factor(roc_long$therapy, levels = therapy_order)

# AUC annotations (black, factual)
auc_annotations <- do.call(rbind, lapply(seq_along(propensity_results), function(i) {
  res <- propensity_results[[i]]
  data.frame(
    therapy = res$therapy,
    label   = c(sprintf("Patient AUC = %.2f", res$patient$auc),
                sprintf("Therapist AUC = %.2f", res$therapist$auc),
                sprintf("Combined AUC = %.2f", res$combined$auc)),
    x       = 0.97,
    y       = c(0.22, 0.12, 0.02)
  )
}))
auc_annotations$therapy <- factor(auc_annotations$therapy, levels = therapy_order)

panel_a <- ggplot(roc_long, aes(x = fpr, y = tpr, color = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_text(data = auc_annotations,
            aes(x = x, y = y, label = label),
            hjust = 1, size = 2.3, color = "black",
            show.legend = FALSE, inherit.aes = FALSE) +
  scale_color_manual(
    name   = NULL,
    values = c("Patient Features"  = col_patient,
               "Therapist Factors" = col_therapist,
               "Combined"          = col_combined)
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
  facet_wrap(~ therapy, nrow = 2) +
  labs(
    title = "Receiver operating characteristic curves for therapy assignment",
    x     = "False positive rate",
    y     = "True positive rate"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title         = element_text(size = 11.5, face = "bold", hjust = 0,
                                      margin = margin(b = 8)),
    strip.text         = element_text(size = 10, face = "bold"),
    strip.background   = element_rect(fill = "grey96", color = "grey80"),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 8, color = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "#F0F0F0", linewidth = 0.2),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    legend.key.width   = unit(1.2, "cm"),
    legend.margin      = margin(t = 0, b = 0),
    panel.spacing      = unit(0.6, "lines"),
    aspect.ratio       = 1,
    plot.margin        = margin(8, 8, 4, 8)
  )

# ============================================================================
# PANEL B: PAIRED AUROC FOREST PLOT WITH DELONG BRACKETS
# ============================================================================

fig2b_data <- do.call(rbind, lapply(seq_along(propensity_results), function(i) {
  res <- propensity_results[[i]]
  auc_p <- res$patient$auc
  auc_t <- res$therapist$auc
  n_pos <- res$n_positive
  n_neg <- res$n_total - res$n_positive
  
  se_hanley <- function(auc, np, nn) {
    q1 <- auc / (2 - auc)
    q2 <- 2 * auc^2 / (1 + auc)
    sqrt((auc * (1 - auc) + (np - 1) * (q1 - auc^2) + (nn - 1) * (q2 - auc^2)) / (np * nn))
  }
  
  se_p <- se_hanley(auc_p, n_pos, n_neg)
  se_t <- se_hanley(auc_t, n_pos, n_neg)
  
  rbind(
    data.frame(Therapy = res$therapy, model = "Patient Features",
               auc = auc_p, ci_lo = auc_p - 1.96 * se_p, ci_hi = auc_p + 1.96 * se_p),
    data.frame(Therapy = res$therapy, model = "Therapist Factors",
               auc = auc_t, ci_lo = auc_t - 1.96 * se_t, ci_hi = auc_t + 1.96 * se_t)
  )
}))

fig2b_data$model <- factor(fig2b_data$model,
                           levels = c("Patient Features", "Therapist Factors"))

therapy_gap <- fig2b_data %>%
  group_by(Therapy) %>%
  summarise(gap = auc[model == "Therapist Factors"] - auc[model == "Patient Features"],
            .groups = "drop") %>%
  arrange(gap)

fig2b_data$Therapy <- factor(fig2b_data$Therapy, levels = therapy_gap$Therapy)
fig2b_data <- fig2b_data %>%
  mutate(y_num  = as.numeric(Therapy),
         y_plot = ifelse(model == "Patient Features", y_num - 0.15, y_num + 0.15))

# DeLong annotations
delong_data <- propensity_summary %>%
  select(Therapy, Patient_AUC, Therapist_AUC, DeLong_p) %>%
  mutate(
    sig_label = ifelse(DeLong_p < 0.001, "p < .001", sprintf("p = %.3f", DeLong_p)),
    gap_label = sprintf("Delta == %.3f", Therapist_AUC - Patient_AUC)
  )
delong_data$Therapy <- factor(delong_data$Therapy, levels = therapy_gap$Therapy)
delong_data <- delong_data %>% mutate(y_num = as.numeric(Therapy))

bracket_x <- 0.80

panel_b <- ggplot(fig2b_data, aes(y = y_plot)) +
  
  # Chance line
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  
  # CIs (thick, colored)
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = model),
                 height = 0, linewidth = 0.9) +
  
  # Point estimates (filled circle with white border)
  geom_point(aes(x = auc, fill = model), 
             size = 3.5, shape = 21, color = "white", stroke = 0.5) +
  
  # Bracket: vertical
  geom_segment(data = delong_data,
               aes(x = bracket_x, xend = bracket_x,
                   y = y_num - 0.15, yend = y_num + 0.15),
               color = "grey60", linewidth = 0.35) +
  
  # Bracket: ticks
  geom_segment(data = delong_data,
               aes(x = bracket_x - 0.005, xend = bracket_x,
                   y = y_num + 0.15, yend = y_num + 0.15),
               color = "grey60", linewidth = 0.35) +
  geom_segment(data = delong_data,
               aes(x = bracket_x - 0.005, xend = bracket_x,
                   y = y_num - 0.15, yend = y_num - 0.15),
               color = "grey60", linewidth = 0.35) +
  
  # Delta
  geom_text(
    data = delong_data,
    aes(x = bracket_x + 0.012, y = y_num + 0.2, label = gap_label),
    hjust = 0,
    size = 2.7,
    color = "black",
    fontface = "bold",
    parse = TRUE
  ) +
  
  # P-value
  geom_text(data = delong_data,
            aes(x = bracket_x + 0.012, y = y_num - 0.2, label = sig_label),
            hjust = 0, size = 2.3, color = "grey50") +
  
  scale_color_manual(
    name   = NULL,
    values = c("Patient Features" = col_patient, "Therapist Factors" = col_therapist)
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Patient Features" = col_patient, "Therapist Factors" = col_therapist)
  ) +
  scale_y_continuous(
    breaks = 1:nrow(therapy_gap),
    labels = levels(fig2b_data$Therapy),
    expand = expansion(add = 0.5)
  ) +
  scale_x_continuous(
    limits = c(0.42, 0.92),
    breaks = seq(0.45, 0.80, 0.05),
    name   = "AUROC (95% CI)"
  ) +
  labs(
    title = "AUROC comparison: patient features vs. therapist factors",
    y     = NULL
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 21, size = 3.5, stroke = 0.5)),
    fill  = "none"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(size = 11.5, face = "bold", hjust = 0,
                                      margin = margin(b = 8)),
    axis.text.y        = element_text(size = 10, color = "black"),
    axis.text.x        = element_text(size = 9, color = "black"),
    axis.title.x       = element_text(size = 10, margin = margin(t = 6)),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#F0F0F0", linewidth = 0.2),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    plot.margin        = margin(8, 10, 4, 8)
  )

# ============================================================================
# COMBINE
# ============================================================================

fig2_combined <- cowplot::plot_grid(
  panel_a,
  panel_b,
  ncol = 1,
  rel_heights = c(0.6, 0.4),
  labels = c("A", "B"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0.0,
  label_y = 1.0
)

ggsave("outputs/figures/figure2_propensity_combined.pdf",
       fig2_combined, width = 7, height = 8, dpi = 300, bg = "white")
ggsave("outputs/figures/figure2_propensity_combined.png",
       fig2_combined, width = 8, height = 8, dpi = 300, bg = "white")

cat("Saved: outputs/figures/figure2_propensity_combined.pdf\n")
# ============================================================================
# FIGURE 3A: SHAP VALUES
# ============================================================================

cat("Creating Figure 3A: SHAP values...\n")

# Compute SHAP values
set.seed(123)
shap_n     <- min(10000, nrow(X_test))
shap_idx   <- sample(1:nrow(X_test), shap_n)
shap_data  <- xgb.DMatrix(X_test[shap_idx, ])
shap_raw   <- predict(xgb_final, shap_data, predcontrib = TRUE, approxcontrib = FALSE)
shap_vals  <- shap_raw[, -ncol(shap_raw)]
feat_names <- colnames(X_test)

# Filter to patient-level features (exclude provider/location fixed effects)
include_patterns <- c(
  "^total_score$", "^risk_high_initial$", "^deterrents_month$",
  "^what_sort_of_reasons$", "^duration_month$", "^adolescent$", "^male$",
  "^frequency_month$", "^are_there_things$",
  "^how_many_times_have_you_had_these_thoughts$",
  "^when_you_have_the_thoughts_how_long_do_they_last$",
  "^could_can_you_stop_thinking_about_killing_yourself",
  "^current_and_past_psychiatric_diagnoses_",
  "^presenting_symptoms_", "^family_history_",
  "^precipitants_stressors_",
  "^internal_protective_factors_", "^external_protective_factors_",
  "^change_in_treatment_", "^dx_group_",
  "^act$", "^cbt$", "^dbt$", "^motivational_interviewing$",
  "^mindfulness$", "^stages_of_change$", "^family_systems$",
  "^therapy_duration_category_", "^delivery_method_", "^session_mode_"
)

exclude_patterns <- c(
  "^therapist_name_", "^location_", "^program_",
  "^pn_month", "^pn_time_block", "^pn_year",
  "^intake_to_pn$", "^days_first_srs$"
)

keep <- sapply(feat_names, function(f) {
  inc <- any(sapply(include_patterns, function(p) grepl(p, f)))
  exc <- any(sapply(exclude_patterns, function(p) grepl(p, f)))
  inc & !exc
})

filtered_shap <- shap_vals[, keep]
mean_abs      <- colMeans(abs(filtered_shap))
names(mean_abs) <- feat_names[keep]
top20         <- names(sort(mean_abs, decreasing = TRUE)[1:20])

# Build plot data
shap_plot_data <- do.call(rbind, lapply(top20, function(feat) {
  fi <- which(feat_names == feat)
  data.frame(feature = feat,
             feature_value = X_test[shap_idx, fi],
             shap_value    = shap_vals[, fi])
}))

shap_plot_data$feature_clean <- case_when(
  shap_plot_data$feature == "risk_high_initial"
  ~ "High suicide risk at intake",
  shap_plot_data$feature == "total_score"
  ~ "C-SSRS total score (0-25)",
  shap_plot_data$feature == "how_many_times_have_you_had_these_thoughts"
  ~ "SI frequency (lifetime)",
  shap_plot_data$feature == "when_you_have_the_thoughts_how_long_do_they_last"
  ~ "SI duration (episode length)",
  shap_plot_data$feature == "frequency_month"
  ~ "SI frequency (past month)",
  shap_plot_data$feature == "precipitants_stressors_chronic_physical_pain_or_other_acute_medical_problem_e_g_cns_disorders"
  ~ "Stressor: Chronic pain / acute medical problem",
  shap_plot_data$feature == "internal_protective_factors_able_to_access_care_willing_to_reach_out"
  ~ "PF (internal): Access to care",
  shap_plot_data$feature == "internal_protective_factors_religious_beliefs"
  ~ "PF (internal): Religious beliefs",
  shap_plot_data$feature == "are_there_things"
  ~ "Protective factors present",
  shap_plot_data$feature == "internal_protective_factors_identifies_reasons_for_living"
  ~ "PF (internal): Reasons for living",
  shap_plot_data$feature == "external_protective_factors_positive_therapeutic_relationships"
  ~ "PF (external): Positive therapeutic relationships",
  shap_plot_data$feature == "could_can_you_stop_thinking_about_killing_yourself_or_wanting_to_die_if_you_want_to"
  ~ "Controllability of suicidal thoughts",
  shap_plot_data$feature == "external_protective_factors_high_academic_achievement"
  ~ "PF (external): High academic achievement",
  shap_plot_data$feature == "male"
  ~ "Male sex",
  shap_plot_data$feature == "precipitants_stressors_social_isolation"
  ~ "Stressor: Social isolation",
  shap_plot_data$feature == "cbt"
  ~ "CBT received",
  shap_plot_data$feature == "dx_group_Trauma_Related_Disorder"
  ~ "Dx: Trauma-related disorder",
  shap_plot_data$feature == "dx_group_Depressive_Disorder"
  ~ "Dx: Depressive disorder",
  shap_plot_data$feature == "adolescent"
  ~ "Adolescent patient",
  shap_plot_data$feature == "mindfulness"
  ~ "Mindfulness received",
  TRUE ~ gsub("_", " ", shap_plot_data$feature)
)
# Normalize feature values for color scale
shap_plot_data <- shap_plot_data %>%
  group_by(feature) %>%
  mutate(feature_norm = (feature_value - min(feature_value, na.rm = TRUE)) /
           (max(feature_value, na.rm = TRUE) - min(feature_value, na.rm = TRUE) + 1e-10)) %>%
  ungroup()

# Order by mean |SHAP|
feat_order <- shap_plot_data %>%
  group_by(feature_clean) %>%
  summarise(mean_abs = mean(abs(shap_value)), .groups = "drop") %>%
  arrange(desc(mean_abs))

shap_plot_data$feature_clean <- factor(shap_plot_data$feature_clean,
                                       levels = rev(feat_order$feature_clean))

panel_a <- ggplot(shap_plot_data, aes(x = shap_value, y = feature_clean)) +
  geom_jitter(aes(color = feature_norm),
              height = 0.15, width = 0, size = 0.8, alpha = 0.7) +
  scale_color_gradient2(
    low = "#6366F1", mid = "grey95", high = "#C1121F",
    midpoint = 0.5, name = "Feature\nvalue\n",
    breaks = c(0, 0.999), labels = c("Low", "High")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = pal$red, linewidth = 0.3) +
  labs(
    title    = "Feature importance for improvement prediction",
    x = "SHAP value (impact on prediction)",
    y = NULL
  ) +
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

# ============================================================================
# FIGURE 3B: DOUBLY ROBUST ATEs
# ============================================================================

cat("Creating Figure 3B: Therapy ATEs...\n")

dr_df <- data.frame(
  therapy   = names(dr_results),
  ate       = sapply(dr_results, function(x) x$ate),
  ci_lower  = sapply(dr_results, function(x) x$ci_lower),
  ci_upper  = sapply(dr_results, function(x) x$ci_upper)
)

dr_df$therapy_clean <- case_when(
  dr_df$therapy == "cbt"                       ~ "CBT",
  dr_df$therapy == "dbt"                       ~ "DBT",
  dr_df$therapy == "act"                       ~ "ACT",
  dr_df$therapy == "motivational_interviewing" ~ "Motivational Interviewing",
  dr_df$therapy == "mindfulness"               ~ "Mindfulness",
  dr_df$therapy == "stages_of_change"          ~ "Stages of Change",
  dr_df$therapy == "family_systems"            ~ "Family Systems",
  TRUE ~ dr_df$therapy
)

dr_df <- dr_df %>% arrange(ate)
dr_df$therapy_clean <- factor(dr_df$therapy_clean, levels = dr_df$therapy_clean)

panel_b <- ggplot(dr_df, aes(x = therapy_clean, y = ate * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = pal$red, linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100),
                width = 0, linewidth = 0.9, color = pal$red) +
  geom_point(size = 3.5, fill = pal$red, color = "white", shape = 21, stroke = 0.5) +
  coord_flip() +
  labs(
    title    = "Bias-adjusted therapy-specific associations with improvement",
    x = NULL,
    y = "Difference in improvement probability (pp)"
  ) +
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

# Combine panels
fig3 <- cowplot::plot_grid(
  panel_a, panel_b + theme(legend.position = "none"),
  ncol = 1, rel_heights = c(1.2, 0.8),
  labels = c("A", "B"), label_size = 14, align = "v"
)
fig3

ggsave("outputs/figures/figure3_shap_and_ate.pdf",
       fig3, width = 10, height = 9, dpi = 300, bg = "white")

ggsave("outputs/figures/figure3_shap_and_ate.png",
       fig3, width = 8.5, height = 7, dpi = 300, bg = "white")

cat("  Saved: outputs/figures/figure3_shap_and_ate.pdf\n")

# ============================================================================
# FIGURE 4: PERSONALIZATION GAINS BY BASELINE RISK
#   Panel A — Mean gain by quartile with individual patient jitter
#   Panel B — Individual gains across baseline probability (scatter + trend)
# ============================================================================



cat("Creating Figure 4: Personalization gains...\n")

# Colors
col_red       <- "#C1121F"
col_red_light <- "#E76F51"
col_amber     <- "#F59E0B"
col_indigo    <- "#6366F1"

# ============================================================================
# PANEL A: MEAN GAIN BY QUARTILE WITH JITTERED INDIVIDUAL POINTS
# ============================================================================

# Add quartile to individual-level data
plot_data_a <- personalization_viz_data %>%
  mutate(
    prob_quartile = cut(
      baseline_prob,
      breaks = quantile(baseline_prob, probs = c(0, 0.25, 0.5, 0.75, 1)),
      labels = c("Q1\n(Lowest)", "Q2", "Q3", "Q4\n(Highest)"),
      include.lowest = TRUE
    )
  )

# Quartile summaries
quartile_summary <- plot_data_a %>%
  group_by(prob_quartile) %>%
  summarise(
    n         = n(),
    mean_gain = mean(gain),
    se_gain   = sd(gain) / sqrt(n()),
    ci_lo     = mean_gain - 1.96 * se_gain,
    ci_hi     = mean_gain + 1.96 * se_gain,
    .groups   = "drop"
  )

panel_a <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  
  # Violin (faint, shows distribution shape)
  geom_violin(data = plot_data_a,
              aes(x = prob_quartile, y = gain),
              fill = "#E76F51", color = NA, alpha = 0.2,
              width = 0.7, scale = "width") +
  geom_jitter(data = plot_data_a,
              aes(x = prob_quartile, y = gain),
              width = 0.25, height = 0, size = 0.8, alpha = 0.07,
              color = 'grey') +
  
  # CI bars
  geom_errorbar(data = quartile_summary,
                aes(x = prob_quartile, ymin = ci_lo, ymax = ci_hi),
                width = 0, linewidth = 0.9, color = col_red) +
  
  # Mean points
  geom_point(data = quartile_summary,
             aes(x = prob_quartile, y = mean_gain),
             size = 4, shape = 21, fill = col_red, color = "white", stroke = 0.5) +
  
  # Value labels
  geom_text(data = quartile_summary,
            aes(x = prob_quartile, y = mean_gain,
                label = sprintf("%.1f pp", mean_gain)),
            vjust = -1.8, size = 3.2, color = "black", fontface = "bold") +
  
  # Sample size
  geom_text(data = quartile_summary,
            aes(x = prob_quartile, y = -1.5,
                label = sprintf("n = %d", n)),
            size = 2.5, color = "grey50") +
  
  scale_y_continuous(
    limits = c(-2.5, 30),
    breaks = seq(0, 30, 5),
    name   = "Predicted personalization gain (pp)"
  ) +
  labs(
    title = "Personalization gains by baseline improvement probability",
    x     = "Baseline predicted probability quartile"
  ) +
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

# ============================================================================
# PANEL B: GAINS BY INITIAL RISK LEVEL
# ============================================================================

risk_summary <- personalization_viz_data %>%
  group_by(initial_risk) %>%
  summarise(
    mean_gain   = mean(gain),
    median_gain = median(gain),
    n           = n(),
    .groups     = "drop"
  )

panel_b <- ggplot(personalization_viz_data,
                  aes(x = initial_risk, y = gain)) +
  geom_violin(aes(fill = initial_risk), alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, outlier.alpha = 0.3,
               color = "#333333", fill = NA) +
  
  # Mean annotation
  geom_text(data = risk_summary,
            aes(x = initial_risk, y = mean_gain,
                label = sprintf("mean = %.1f pp", mean_gain)),
            hjust = -0.3, size = 3, color = "black", fontface='bold') +
  
  scale_fill_manual(values = c("Moderate" = col_amber, "High" = col_red)) +
  labs(
    title = "Predicted personalization gains by initial suicide risk level",
    x     = "Initial risk level",
    y     = "Predicted personalization gain (pp)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(size = 11.5, face = "bold", hjust = 0,
                                      margin = margin(b = 8)),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9, color = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.2),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.margin        = margin(8, 8, 4, 8)
  )

# ============================================================================
# COMBINE
# ============================================================================

fig4 <- cowplot::plot_grid(
  panel_a, panel_b,
  ncol = 1, rel_heights = c(0.5, 0.5),
  labels = c("A", "B"),
  label_size = 15, label_fontface = "bold",
  label_x = 0.0, label_y = 1.0
)
fig4
ggsave("outputs/figures/figure4_personalization_gains.pdf",
       fig4, width = 8, height = 10, dpi = 300, bg = "white")
ggsave("outputs/figures/figure4_personalization_gains.png",
       fig4, width = 6, height = 6, dpi = 300, bg = "white")


cat("  Saved: outputs/figures/figure4_personalization_gains.pdf\n")
# ============================================================================
# SUPPLEMENTAL: eFIGURE 2 — TEMPORAL TRENDS
# ============================================================================

cat("Creating eFigure 2: Temporal trends...\n")

temporal_data <- data %>%
  mutate(month_year = floor_date(as.Date(admission_date), "month")) %>%
  group_by(month_year) %>%
  summarise(
    n             = n(),
    improve_rate  = mean(improve, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Training vs test split date
split_date <- max(train_data$admission_date)

p_temporal <- ggplot(temporal_data, aes(x = month_year, y = improve_rate)) +
  geom_vline(xintercept = as.numeric(as.Date(split_date)),
             linetype = "dashed", color = pal$red, linewidth = 0.5) +
  geom_point(color = pal$teal, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = pal$indigo,
              fill = pal$indigo, alpha = 0.15, linewidth = 0.8) +
  geom_hline(yintercept = mean(train_data$improve, na.rm = TRUE) * 100,
             linetype = "dotted", color = "grey50", linewidth = 0.4) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "3 months") +
  annotate("text", x = as.Date(split_date) - 60, y = 95,
           label = "Training", hjust = 1, size = 3.5, color = "grey40") +
  annotate("text", x = as.Date(split_date) + 60, y = 95,
           label = "Test", hjust = 0, size = 3.5, color = "grey40") +
  labs(
    x     = "Admission Date",
    y     = "Proportion Achieving Risk Reduction (%)"
  ) +
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
p_temporal
ggsave("outputs/figures/efigure2_temporal_trends.png",
       p_temporal, width = 9, height = 5, dpi = 300, bg = "white")

cat("  Saved: outputs/figures/efigure2_temporal_trends.png\n")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("ALL FIGURES COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")