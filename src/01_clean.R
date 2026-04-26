# ============================================================================
# 01_clean.R
# Data Preparation for Personalized Psychotherapy Analysis
# ============================================================================
#
# This script reads and merges three data sources from Discovery Behavioral
# Health (DBH) electronic health records:
#   1. MH SA  — Suicide assessments (SAFE-T/C-SSRS) at intake
#   2. MH SRS — Suicide risk re-screens during treatment
#   3. MH PN  — Progress notes documenting psychotherapy modalities
#   4. Demo   — Patient demographics, program, and diagnosis
#
# Output: analytic dataset `data` containing treatment episodes for patients
#         at moderate-to-high suicide risk with documented psychotherapy
#         exposure and follow-up risk assessment.
# ============================================================================

# ============================================================================
# 1. SUICIDE ASSESSMENTS AT INTAKE (MH SA)
# ============================================================================

mh_sa <- read_excel("~/Sue Goldie Dropbox/Jacob Jameson/DBH data/MH SA.xlsx")
mh_sa <- clean_names(mh_sa)

# Rename unwieldy risk-level columns
mh_sa <- mh_sa %>%
  rename(
    risk_final    = final_b_risk_level_b_including_any_change_based_on_clinical_judgment_if_applicable_u,
    risk_high_click = please_click_here_for_high_risk_level_if_applicable,
    risk_mod_click  = please_click_here_for_moderate_risk_level_if_applicable,
    risk_low_click  = please_click_here_for_low_risk_level_if_applicable
  )

# Derive initial risk level from structured click fields, falling back to
# the free-text risk_level field when click fields are absent
mh_sa <- mh_sa %>%
  rowwise() %>%
  mutate(
    risk_level_initial = case_when(
      !is.na(risk_low_click)  & !str_trim(tolower(as.character(risk_low_click)))  %in% c("", "na", "n/a") ~ "Low",
      !is.na(risk_mod_click)  & !str_trim(tolower(as.character(risk_mod_click)))  %in% c("", "na", "n/a") ~ "Moderate",
      !is.na(risk_high_click) & !str_trim(tolower(as.character(risk_high_click))) %in% c("", "na", "n/a") ~ "High",
      !is.na(risk_level) & str_detect(tolower(risk_level), "low")      ~ "Low",
      !is.na(risk_level) & str_detect(tolower(risk_level), "moderate") ~ "Moderate",
      !is.na(risk_level) & str_detect(tolower(risk_level), "high")     ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    risk_level_initial = factor(risk_level_initial,
                                levels = c("Low", "Moderate", "High"),
                                ordered = TRUE)
  )

cat("Initial risk level distribution (MH SA):\n")
print(table(mh_sa$risk_level_initial, useNA = "ifany"))

# ============================================================================
# 2. SUICIDE RISK RE-SCREENS (MH SRS)
# ============================================================================

srs <- read_excel("~/Sue Goldie Dropbox/Jacob Jameson/DBH data/MH SRS.xlsx") %>%
  clean_names() %>%
  rename(
    raw_risk_level     = final_b_risk_level_b_including_any_change_based_on_clinical_judgment_if_applicable_u,
    high_risk_flag1    = documentation_follow_up_create_resolve_urgent_issue_red_flag_in_emr_high_risk,
    high_risk_flag2    = documentation_follow_up_create_suicide_risk_treatment_plan_high_risk,
    moderate_risk_flag = documentation_follow_up_include_suicide_risk_reduction_interventions_in_appropriate_treatment_plan_moderate_risk,
    low_risk_flag      = documentation_follow_up_n_a_only_applicable_if_low_risk
  )

# Normalize flag columns to logical
srs <- srs %>%
  mutate(
    high_risk_flag1    = tolower(str_trim(as.character(high_risk_flag1)))    %in% c("true", "t", "yes", "y", "1"),
    high_risk_flag2    = tolower(str_trim(as.character(high_risk_flag2)))    %in% c("true", "t", "yes", "y", "1"),
    moderate_risk_flag = tolower(str_trim(as.character(moderate_risk_flag))) %in% c("true", "t", "yes", "y", "1"),
    low_risk_flag      = tolower(str_trim(as.character(low_risk_flag)))      %in% c("true", "t", "yes", "y", "1")
  )

# Derive risk level from structured text field
srs <- srs %>%
  mutate(
    initial_risk_level = case_when(
      str_to_lower(str_trim(raw_risk_level)) == "low suicide risk"      ~ "Low",
      str_to_lower(str_trim(raw_risk_level)) == "moderate suicide risk" ~ "Moderate",
      str_to_lower(str_trim(raw_risk_level)) == "high suicide risk"     ~ "High",
      TRUE ~ NA_character_
    )
  )

# Derive risk level from documentation flags (priority: high > moderate > low)
srs <- srs %>%
  mutate(
    flag_based_risk = case_when(
      high_risk_flag1 | high_risk_flag2 ~ "High",
      moderate_risk_flag                ~ "Moderate",
      low_risk_flag                     ~ "Low",
      TRUE                              ~ NA_character_
    )
  )

# Coalesce: prefer structured text, fall back to flag-based
srs <- srs %>%
  mutate(
    risk_level = coalesce(initial_risk_level, flag_based_risk),
    risk_level = factor(risk_level, levels = c("Low", "Moderate", "High"), ordered = TRUE)
  ) %>%
  select(-raw_risk_level, -high_risk_flag1, -high_risk_flag2,
         -moderate_risk_flag, -low_risk_flag,
         -initial_risk_level, -flag_based_risk) %>%
  arrange(master_id, admission_date, evaluation_date)

cat("\nSRS risk level distribution:\n")
print(table(srs$risk_level, useNA = "ifany"))

# Keep first and last re-screen per admission episode
srs <- srs %>%
  group_by(master_id, admission_date) %>%
  filter(evaluation_date == min(evaluation_date) |
           evaluation_date == max(evaluation_date)) %>%
  mutate(
    srs_first = as.integer(evaluation_date == min(evaluation_date)),
    srs_last  = as.integer(evaluation_date == max(evaluation_date))
  ) %>%
  ungroup()

# Collapse to one row per admission with first/last risk levels and timing
srs <- srs %>%
  select(master_id, admission_date, evaluation_date, risk_level, srs_first, srs_last) %>%
  group_by(master_id, admission_date) %>%
  reframe(
    risk_level_srs_first = risk_level[srs_first == 1][1],
    risk_level_srs_last  = risk_level[srs_last == 1][1],
    days_first_srs = as.numeric(difftime(evaluation_date[srs_first == 1][1],
                                         admission_date, units = "days")),
    days_last_srs  = as.numeric(difftime(evaluation_date[srs_last == 1][1],
                                         admission_date, units = "days")),
    .groups = "drop"
  ) %>%
  unique()

# ============================================================================
# 3. PROGRESS NOTES — PSYCHOTHERAPY MODALITIES (MH PN)
# ============================================================================

pn <- read_excel("~/Sue Goldie Dropbox/Jacob Jameson/DBH data/MH PN.xlsx")
pn <- clean_names(pn)

# Rename therapy modality columns to short names
pn <- pn %>%
  rename(
    pn_evaluation_date       = evaluation_date,
    therapist                = staff_signature_1,
    act                      = evidence_based_modalities_employed_act,
    cbt                      = evidence_based_modalities_employed_cbt,
    dbt                      = evidence_based_modalities_employed_dbt,
    motivational_interviewing = evidence_based_modalities_employed_motivational_interviewing,
    mindfulness              = evidence_based_modalities_employed_mindfulness_techniques,
    stages_of_change         = evidence_based_modalities_employed_stages_of_change,
    family_systems           = evidence_based_modalities_employed_family_systems,
    trauma_informed          = evidence_based_modalities_employed_trauma_informed_strategies_inc_emdr
  ) %>%
  select(therapist, pn_evaluation_date, master_id,
         act:trauma_informed, session_type, time_service_started_ended)

# Restrict to patients with a suicide assessment
pn <- pn %>%
  filter(master_id %in% unique(mh_sa$master_id))

# Merge admission/discharge dates from MH SA
pn <- merge(
  select(mh_sa, master_id, admission_date, discharge_date),
  pn,
  by = "master_id", all.x = TRUE, suffixes = c("_mhsa", "_pn")
)

# Keep only progress notes within the admission–discharge window
pn <- pn %>%
  filter(pn_evaluation_date >= admission_date & pn_evaluation_date <= discharge_date)

# Keep the first progress note per admission (initial therapy session)
pn <- pn %>%
  arrange(master_id, admission_date, pn_evaluation_date) %>%
  group_by(master_id, admission_date, discharge_date) %>%
  slice_head(n = 1) %>%
  ungroup()

# --- Parse session duration from time range string ---
get_duration_minutes <- function(x) {
  if (suppressWarnings(!is.na(as.numeric(x)))) return(NA_real_)
  parts <- strsplit(x, "/")[[1]]
  if (length(parts) != 2) return(NA_real_)
  start_time <- parse_date_time(parts[1], orders = "ymd HMS z", tz = "UTC")
  end_time   <- parse_date_time(parts[2], orders = "ymd HMS z", tz = "UTC")
  if (is.na(start_time) | is.na(end_time)) return(NA_real_)
  return(as.numeric(difftime(end_time, start_time, units = "mins")))
}

pn$therapy_duration_minutes <- sapply(pn$time_service_started_ended, get_duration_minutes)

# Categorize session duration
pn <- pn %>%
  mutate(
    therapy_duration_category = case_when(
      therapy_duration_minutes < 30 ~ "00-30m",
      therapy_duration_minutes >= 30 & therapy_duration_minutes < 45  ~ "30-45m",
      therapy_duration_minutes >= 45 & therapy_duration_minutes <= 60 ~ "45-60m",
      therapy_duration_minutes > 60 ~ "60m+",
      TRUE ~ NA_character_
    ),
    therapy_duration_category = factor(therapy_duration_category,
                                       levels = c("00-30m", "30-45m", "45-60m", "60m+"),
                                       ordered = TRUE)
  )

# Derive delivery method and session mode from session_type
pn <- pn %>%
  mutate(
    delivery_method = case_when(
      str_detect(session_type, "^Telehealth") ~ "Telehealth",
      TRUE                                    ~ "In-person"
    ),
    session_mode = case_when(
      str_detect(session_type, "Individual")         ~ "Individual",
      str_detect(session_type, "Family")             ~ "Family",
      str_detect(session_type, "Collateral Contact") ~ "Collateral Contact",
      TRUE                                           ~ NA_character_
    ),
    delivery_method = factor(delivery_method, levels = c("In-person", "Telehealth")),
    session_mode    = factor(session_mode,    levels = c("Individual", "Family", "Collateral Contact"))
  ) %>%
  select(-session_type, -time_service_started_ended, -therapy_duration_minutes)

# ============================================================================
# 4. MERGE ALL DATA SOURCES
# ============================================================================

mh_sa <- merge(mh_sa, srs, 
               by = c("master_id", "admission_date"), all.x = TRUE)

mh_sa <- merge(mh_sa, pn,  
               by = c("master_id", "admission_date", "discharge_date"),
               all.x = TRUE)

# Drop episodes without a follow-up risk re-screen
mh_sa <- mh_sa %>%
  filter(!is.na(risk_level_srs_first))

# Compute days from intake assessment to first therapy session
mh_sa <- mh_sa %>%
  mutate(intake_to_pn = pmax(
    as.numeric(difftime(pn_evaluation_date, 
                        evaluation_date, 
                        units = "days")), 0))

# ============================================================================
# 5. FEATURE CLEANING
# ============================================================================

# Keep relevant columns (intake assessment + SRS + PN features)
mh_sa <- mh_sa[, c(1:88, 155:159, 161:174)]

# Convert logical columns to integer (0/1)
mh_sa <- mh_sa %>%
  mutate(across(where(is.logical), as.integer))

# Parse ordinal C-SSRS component scores from text labels
cssrs_ordinal_cols <- c(
  "how_many_times_have_you_had_these_thoughts",
  "when_you_have_the_thoughts_how_long_do_they_last",
  "frequency_month",
  "could_can_you_stop_thinking_about_killing_yourself_or_wanting_to_die_if_you_want_to",
  "are_there_things",
  "duration_month",
  "what_sort_of_reasons",
  "controllability_month",
  "deterrents_month"
)

mh_sa <- mh_sa %>%
  mutate(across(all_of(cssrs_ordinal_cols), ~ as.numeric(str_extract(.x, "\\d+")))) %>%
  mutate(across(all_of(cssrs_ordinal_cols), ~ replace_na(.x, 0)))

# Fill missing therapy indicators with 0 (no modality delivered)
therapy_cols <- c("act", "cbt", "dbt", "mindfulness",
                  "motivational_interviewing", "stages_of_change",
                  "family_systems", "trauma_informed")

mh_sa <- mh_sa %>%
  mutate(across(all_of(therapy_cols), ~ replace_na(.x, 0)))

# Drop columns with >5% missing values
cols_high_missing <- names(which(colMeans(is.na(mh_sa)) > 0.05))
if (length(cols_high_missing) > 0) {
  cat(sprintf("\nDropping %d columns with >5%% missing: %s\n",
              length(cols_high_missing), paste(cols_high_missing, collapse = ", ")))
  mh_sa <- mh_sa[, !(names(mh_sa) %in% cols_high_missing)]
}

# Require non-missing initial and follow-up risk levels
mh_sa <- mh_sa %>%
  filter(!is.na(risk_level_initial), !is.na(risk_level_srs_first))

# Extract therapist last name from signature field
mh_sa$therapist_name <- sub(",.*", "", mh_sa$therapist)
mh_sa <- mh_sa %>%
  filter(!is.na(therapist_name) & therapist_name != "")

# Recode access to lethal means as binary
mh_sa <- mh_sa %>%
  mutate(
    access_to_lethal_methods_does_patient_have_access_to_means_including_firearms_in_the_home = case_when(
      access_to_lethal_methods_does_patient_have_access_to_means_including_firearms_in_the_home == "Yes" ~ 1,
      access_to_lethal_methods_does_patient_have_access_to_means_including_firearms_in_the_home == "No"  ~ 0,
      TRUE ~ NA_real_
    )
  )

# Fill remaining missing values in binary clinical fields with 0
mh_sa <- mh_sa %>%
  mutate(across(6:66, ~ replace_na(.x, 0)))

# Require non-negative SRS timing
mh_sa <- mh_sa %>%
  filter((is.na(days_last_srs)  | days_last_srs  >= 0) &
           (is.na(days_first_srs) | days_first_srs >= 0))

# Extract temporal features from progress note date
mh_sa <- mh_sa %>%
  mutate(
    pn_year       = format(pn_evaluation_date, "%Y"),
    pn_month      = format(pn_evaluation_date, "%m"),
    pn_hour       = hour(pn_evaluation_date),
    pn_time_block = case_when(
      pn_hour <  4 ~ "00:00-03:59",
      pn_hour <  8 ~ "04:00-07:59",
      pn_hour < 12 ~ "08:00-11:59",
      pn_hour < 16 ~ "12:00-15:59",
      pn_hour < 20 ~ "16:00-19:59",
      TRUE         ~ "20:00-23:59"
    )
  ) %>%
  select(-pn_evaluation_date)

# Drop raw date columns no longer needed
mh_sa <- mh_sa %>%
  select(-date, -evaluation_date) %>%
  mutate(session_mode = ifelse(is.na(session_mode), "Individual", session_mode))

# ============================================================================
# 6. DEMOGRAPHICS, PROGRAM, AND DIAGNOSIS
# ============================================================================

demo <- read_excel("~/Sue Goldie Dropbox/Jacob Jameson/DBH data/demo.xlsx") %>%
  clean_names() %>%
  filter(division == "MH") %>%
  mutate(program = ifelse(php_days > 0, "PHP", program)) %>%
  select(master_id, admission_date, location, program, age_group, prim_mh_dx, sex_fs)

mh_sa <- merge(mh_sa, demo, by = c("master_id", "admission_date"), all.x = TRUE)

# Map ICD diagnoses to analytic categories
mh_sa <- mh_sa %>%
  mutate(
    dx_group = case_when(
      str_detect(prim_mh_dx, regex("depress|mood disorder", ignore_case = TRUE))
      ~ "Depressive Disorder",
      str_detect(prim_mh_dx, regex("bipolar|cyclothymic", ignore_case = TRUE))
      ~ "Bipolar Disorder",
      str_detect(prim_mh_dx, regex("anxiety|panic|phobia", ignore_case = TRUE))
      ~ "Anxiety Disorder",
      str_detect(prim_mh_dx, regex("trauma|stress|ptsd|reactive attachment", ignore_case = TRUE))
      ~ "Trauma-Related Disorder",
      str_detect(prim_mh_dx, regex("use disorder|substance|alcohol|cannabis|opioid|cocaine|amphetamine|hallucinogen|tobacco", ignore_case = TRUE))
      ~ "Substance Use Disorder",
      str_detect(prim_mh_dx, regex("adhd|autism|neurodevelopmental|intellectual", ignore_case = TRUE))
      ~ "Neurodevelopmental Disorder",
      str_detect(prim_mh_dx, regex("personality", ignore_case = TRUE))
      ~ "Personality Disorder",
      str_detect(prim_mh_dx, regex("schizo|psychotic|delusional", ignore_case = TRUE))
      ~ "Psychotic Disorder",
      str_detect(prim_mh_dx, regex("anorexia|bulimia|eating|feeding", ignore_case = TRUE))
      ~ "Eating Disorder",
      str_detect(prim_mh_dx, regex("disruptive|conduct|oppositional|impulse|explosive", ignore_case = TRUE))
      ~ "Disruptive/Impulse Disorder",
      str_detect(prim_mh_dx, regex("obsessive|ocd|trichotillomania", ignore_case = TRUE))
      ~ "OCD and Related",
      str_detect(prim_mh_dx, regex("insomnia|sleep", ignore_case = TRUE))
      ~ "Sleep Disorder",
      str_detect(prim_mh_dx, regex("somatic", ignore_case = TRUE))
      ~ "Somatic Symptom Disorder",
      str_detect(prim_mh_dx, regex("relational|parent child|relationship", ignore_case = TRUE))
      ~ "Relational/Other V-code",
      str_detect(prim_mh_dx, regex("no dx", ignore_case = TRUE))
      ~ "No Diagnosis",
      TRUE ~ "Other/Unspecified"
    )
  ) %>%
  select(-prim_mh_dx)

# Require non-missing program and location
mh_sa <- mh_sa %>%
  filter(!is.na(program) & program != "",
         !is.na(location) & location != "")

# ============================================================================
# 7. DERIVE ANALYTIC VARIABLES AND APPLY FINAL INCLUSION CRITERIA
# ============================================================================

# Derive convenience variables
mh_sa <- mh_sa %>%
  mutate(
    adolescent = as.integer(age_group == "Adolescent"),
    male       = as.integer(sex_fs == "Male")
  )

# Restrict to moderate-to-high risk (primary analytic cohort)
data <- mh_sa %>%
  filter(risk_level_initial != "Low")

# Define primary outcome: improvement in suicide risk stratum
risk_order <- c("Low" = 1, "Moderate" = 2, "High" = 3)

data <- data %>%
  mutate(
    risk_initial_num = risk_order[as.character(risk_level_initial)],
    risk_first_num   = risk_order[as.character(risk_level_srs_first)],
    improve          = as.integer(risk_first_num < risk_initial_num),
    risk_change      = risk_initial_num - risk_first_num
  )

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("ANALYTIC COHORT SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("  Total treatment episodes: %d\n", nrow(data)))
cat(sprintf("  Unique patients: %d\n", length(unique(data$master_id))))
cat(sprintf("  Unique therapists: %d\n", length(unique(data$therapist_name))))
cat(sprintf("  Unique facilities: %d\n", length(unique(data$location))))
cat(sprintf("  Moderate risk: %d (%.1f%%)\n",
            sum(data$risk_level_initial == "Moderate"),
            100 * mean(data$risk_level_initial == "Moderate")))
cat(sprintf("  High risk: %d (%.1f%%)\n",
            sum(data$risk_level_initial == "High"),
            100 * mean(data$risk_level_initial == "High")))
cat(sprintf("  Improvement rate: %.1f%%\n", 100 * mean(data$improve, na.rm = TRUE)))
cat(rep("=", 70), "\n", sep = "")