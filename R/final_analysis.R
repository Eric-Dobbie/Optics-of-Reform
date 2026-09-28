# ==============================================================================
# final_analysis.R
# Confirmatory analysis pipeline for the full experiment
# "The Optics of Reform: Diversification and Trust in Policing"
#
# This is the pre-registered analysis script, written to run on the Qualtrics
# export as delivered (see Data/Pilot/ for the column structure). Point DATA_PATH
# at the final dataset and set SMOKE_TEST <- FALSE to produce the real tables and
# figures.
#
# Design: 3-arm RCT (Control / Race / Gender) block-randomized by race x party.
#   H1: Black respondents report higher trust under Race Treatment vs. Control
#   H2: Women report higher trust under Gender Treatment vs. Control
#
# Usage:
#   Rscript R/final_analysis.R                      # uses config below
#   Rscript R/final_analysis.R path/to/final.csv    # override data path
#
# SMOKE_TEST = TRUE runs every stage end-to-end and prints only a PASS/FAIL
# checklist (no estimates), writing outputs to a scratch directory. This is the
# mode used to validate the script against the pilot before the full data lands.
# ==============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

DATA_PATH  <- "Data/Pilot/Trust_Policing_Pilot_June 17, 2026_17.06.csv"
SMOKE_TEST <- TRUE                    # TRUE = validation run, suppress estimates
OUTPUT_DIR <- if (SMOKE_TEST) file.path(tempdir(), "oor_smoke") else "R/output"

# Allow a data path override from the command line
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && nzchar(.args[1])) DATA_PATH <- .args[1]

# Exclusion / design parameters (pre-registered)
MIN_DURATION_SEC <- 120               # exclude respondents faster than 2 minutes
ATTENTION_ANSWER <- "Somewhat agree"  # correct response to the Reading check

# ── 1. Setup ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(estimatr)
  library(broom)
})

# Optional packages: the pipeline degrades gracefully if these are absent so the
# structural validation still completes. They are required for the real run.
HAS_MODELSUMMARY <- requireNamespace("modelsummary", quietly = TRUE)
HAS_GRF          <- requireNamespace("grf",          quietly = TRUE)
HAS_MEDIATION    <- requireNamespace("mediation",    quietly = TRUE)

dir.create(file.path(OUTPUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# Lightweight validation harness -----------------------------------------------
.stage_log <- list()
run_stage <- function(name, fn) {
  res <- tryCatch({ val <- fn(); list(ok = TRUE, value = val, msg = "") },
                  error = function(e) list(ok = FALSE, value = NULL,
                                           msg = conditionMessage(e)))
  .stage_log[[name]] <<- res
  status <- if (res$ok) "PASS" else "FAIL"
  cat(sprintf("  [%s] %-42s %s\n", status, name,
              if (!res$ok) paste("->", res$msg) else ""))
  invisible(res$value)   # never auto-print stage internals (avoids leaking estimates)
}

# In SMOKE_TEST mode we never print estimates, only structural facts.
report <- function(...) if (!SMOKE_TEST) cat(...)

# ── 2. Data loading and cleaning ──────────────────────────────────────────────

likert_map <- c(
  "Strongly disagree" = 1, "Disagree" = 2, "Somewhat disagree" = 3,
  "Neither agree nor disagree" = 4, "Somewhat agree" = 5,
  "Agree" = 6, "Strongly agree" = 7
)
recode_likert <- function(x) unname(likert_map[trimws(x)])

primary_outcomes   <- c("Trust_General", "Trust_Self", "Trust_Hiring")
secondary_outcomes <- c("Procedural_Justice", "Police_resources")
likert_vars        <- c(primary_outcomes, secondary_outcomes,
                        "Women in police", "Minority_police")

clean_data <- run_stage("load & clean", function() {
  raw <- read_csv(DATA_PATH, show_col_types = FALSE)
  # Qualtrics exports carry two extra header rows (question text + import IDs)
  if (nrow(raw) > 0 && raw[[1]][1] == "Start Date") raw <- raw[-c(1, 2), ]

  df <- raw |>
    mutate(
      across(all_of(likert_vars), recode_likert),
      arm      = factor(na_if(Treatment, ""),
                        levels = c("Control", "Race", "Gender")),
      duration = as.numeric(`Duration (in seconds)`),
      # Treatment indicators
      treat_race   = as.integer(arm == "Race"),
      treat_gender = as.integer(arm == "Gender"),
      # Subgroup flags (H1 / H2)
      black = as.integer(Race == "Black or African-American"),
      woman = as.integer(Gender == "Female"),
      # Covariates
      age        = suppressWarnings(as.numeric(Age)),
      democrat   = as.integer(Party_ID == "Democrat"),
      republican = as.integer(Party_ID == "Republican"),
      ideology   = recode(`Ideology Rating`,
                          "Very liberal" = 1, "Liberal" = 2, "Somewhat liberal" = 3,
                          "Moderate" = 4, "Somewhat conservative" = 5,
                          "Conservative" = 6, "Very conservative" = 7,
                          .default = NA_real_),
      educ_years = recode(Education,
                          "Less than high school" = 10, "High school graduate" = 12,
                          "Some college" = 14, "2 year degree" = 14,
                          "4 year degree" = 16, "Professional degree" = 18,
                          "Doctorate" = 20, .default = NA_real_),
      uof_contact = as.integer(Police_experience == "Yes"),
      att_women   = `Women in police`,
      att_poc     = Minority_police,
      cynicism    = 8 - att_poc,
      # Manipulation / attention checks
      attn_pass   = as.integer(trimws(Reading) == ATTENTION_ANSWER),
      block       = paste(Race, Party_ID, sep = "_")
    )

  # Pre-registered exclusions applied to the ANALYSIS sample
  df_analysis <- df |>
    filter(!is.na(arm), !is.na(Trust_General),
           duration >= MIN_DURATION_SEC, attn_pass == 1)

  list(full = df, analysis = df_analysis)
})

df_full     <- clean_data$full
df          <- clean_data$analysis
report(sprintf("Analysis N = %d (from %d raw)\n", nrow(df), nrow(df_full)))

# Covariate set used throughout
covariates <- c("age", "ideology", "educ_years", "uof_contact",
                "att_women", "att_poc", "democrat", "republican")
cov_str    <- paste(covariates, collapse = " + ")

# ── 3. Balance check (Table 1) ────────────────────────────────────────────────

run_stage("balance table", function() {
  bal_covs <- c("age", "ideology", "educ_years", "uof_contact",
                "att_women", "att_poc", "democrat", "republican", "black", "woman")
  bal_models <- map(bal_covs, ~ lm(reformulate("arm", .x), data = df))
  names(bal_models) <- bal_covs

  # Omnibus balance F-tests: predict each treatment indicator from all covariates
  omni <- map_dfr(c("treat_race", "treat_gender"), function(tv) {
    m  <- lm(reformulate(bal_covs, tv), data = df)
    fs <- summary(m)$fstatistic
    tibble(arm = tv, F = fs[1],
           p = pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
  })

  if (!SMOKE_TEST && HAS_MODELSUMMARY) {
    modelsummary::modelsummary(
      bal_models, output = file.path(OUTPUT_DIR, "tables/table1_balance.tex"),
      stars = TRUE, title = "Covariate Balance Across Treatment Arms",
      gof_map = c("nobs", "r.squared"), fmt = 2
    )
  }
  list(n_models = length(bal_models), omnibus = omni)
})

# ── 4. Primary OLS: H1 and H2 (Table 2) ──────────────────────────────────────
# Trust ~ treat_race + treat_gender + black + woman
#         + treat_race:black (H1) + treat_gender:woman (H2) + covariates
# HC2 robust SEs.

fit_primary <- function(outcome, data) {
  f <- reformulate(
    c("treat_race", "treat_gender", "black", "woman",
      "treat_race:black", "treat_gender:woman", covariates),
    response = outcome)
  lm_robust(f, data = data, se_type = "HC2")
}

models_all <- run_stage("primary + secondary OLS", function() {
  m <- map(c(primary_outcomes, secondary_outcomes), fit_primary, data = df)
  names(m) <- c(primary_outcomes, secondary_outcomes)

  if (!SMOKE_TEST && HAS_MODELSUMMARY) {
    coef_map <- c(treat_race = "Race Treatment", treat_gender = "Gender Treatment",
                  black = "Black respondent", woman = "Woman respondent",
                  "treat_race:black" = "Race x Black (H1)",
                  "treat_gender:woman" = "Gender x Woman (H2)")
    modelsummary::modelsummary(
      m[primary_outcomes],
      output = file.path(OUTPUT_DIR, "tables/table2_primary.tex"),
      coef_map = coef_map, stars = TRUE,
      title = "Effect of Diversity Framing on Trust (Primary Outcomes)",
      gof_map = c("nobs", "r.squared"), fmt = 3)
    modelsummary::modelsummary(
      m[secondary_outcomes],
      output = file.path(OUTPUT_DIR, "tables/table3_secondary.tex"),
      coef_map = coef_map, stars = TRUE,
      title = "Effect on Procedural Justice and Resource Support",
      gof_map = c("nobs", "r.squared"), fmt = 3)
  }
  m
})

# ── 5. Subgroup-only estimates (H1 within Black, H2 within women) ─────────────

run_stage("subgroup OLS (H1/H2)", function() {
  m_h1 <- lm_robust(reformulate(c("treat_race", covariates), "Trust_General"),
                    data = filter(df, black == 1), se_type = "HC2")
  m_h2 <- lm_robust(reformulate(c("treat_gender", covariates), "Trust_General"),
                    data = filter(df, woman == 1), se_type = "HC2")
  list(h1_n = m_h1$nobs, h2_n = m_h2$nobs)
})

# ── 6. Multiple testing correction (Holm, primary family) ─────────────────────

run_stage("Holm correction", function() {
  p_h1 <- map_dbl(models_all[primary_outcomes],
                  ~ tidy(.x) |> filter(term == "treat_race:black") |> pull(p.value))
  p_h2 <- map_dbl(models_all[primary_outcomes],
                  ~ tidy(.x) |> filter(term == "treat_gender:woman") |> pull(p.value))
  holm <- tibble(
    outcome    = rep(primary_outcomes, 2),
    hypothesis = c(rep("H1", 3), rep("H2", 3)),
    p_raw      = c(p_h1, p_h2),
    p_holm     = p.adjust(c(p_h1, p_h2), method = "holm")
  )
  if (!SMOKE_TEST)
    write_csv(holm, file.path(OUTPUT_DIR, "tables/holm_correction.csv"))
  list(n_tests = nrow(holm))
})

# ── 7. Heterogeneous treatment effects (causal forest) ───────────────────────

run_stage("HTE causal forest", function() {
  if (!HAS_GRF) stop("grf not installed (optional; required for real run)")
  hte_covs <- c("age", "ideology", "educ_years", "uof_contact",
                "att_women", "att_poc", "democrat", "republican",
                "black", "woman", "cynicism")
  cf_df <- df |> select(all_of(c("Trust_General", "treat_race", hte_covs))) |> drop_na()
  X  <- as.matrix(cf_df[hte_covs])
  cf <- grf::causal_forest(X = X, Y = cf_df$Trust_General,
                           W = cf_df$treat_race, num.trees = 2000, seed = 42)
  blp <- grf::best_linear_projection(
    cf, A = X[, c("black", "woman", "democrat", "uof_contact", "cynicism")])
  if (!SMOKE_TEST) {
    ate <- grf::average_treatment_effect(cf)
    saveRDS(list(blp = blp, ate = ate),
            file.path(OUTPUT_DIR, "tables/hte_race_trustgeneral.rds"))
  }
  list(n = nrow(cf_df))
})

# ── 8. Mechanism test: mediation via procedural justice ──────────────────────

run_stage("mediation (procedural justice)", function() {
  if (!HAS_MEDIATION) stop("mediation not installed (optional; required for real run)")
  df_black <- filter(df, black == 1)
  if (nrow(df_black) < 20) stop("subgroup too small for mediation in this dataset")
  med_m <- lm(Procedural_Justice ~ treat_race, data = df_black)
  out_m <- lm(Trust_General ~ treat_race + Procedural_Justice, data = df_black)
  med <- mediation::mediate(med_m, out_m, treat = "treat_race",
                            mediator = "Procedural_Justice",
                            robustSE = TRUE, sims = if (SMOKE_TEST) 50 else 1000)
  if (!SMOKE_TEST)
    saveRDS(summary(med), file.path(OUTPUT_DIR, "tables/mediation_h1.rds"))
  list(n = nrow(df_black))
})

# ── 9. Robustness: block FE and full (unexcluded) sample ─────────────────────

run_stage("robustness specs", function() {
  # Block fixed effects
  fe <- map(primary_outcomes, function(y) {
    f <- reformulate(c("treat_race", "treat_gender", "black", "woman",
                       "treat_race:black", "treat_gender:woman",
                       covariates, "block"), response = y)
    lm_robust(f, data = df, se_type = "HC2")
  })
  # Full sample (no exclusions): re-derive on df_full rows with an arm + outcome
  df_all <- df_full |> filter(!is.na(arm), !is.na(Trust_General))
  full <- map(primary_outcomes, fit_primary, data = df_all)
  list(fe_n = length(fe), full_n = nrow(df_all))
})

# ── 10. Coefficient plot (Figure 2) ──────────────────────────────────────────

run_stage("coefficient plot", function() {
  cd <- map_dfr(primary_outcomes, function(y) {
    tidy(models_all[[y]]) |>
      filter(term %in% c("treat_race:black", "treat_gender:woman")) |>
      mutate(outcome = y,
             label = recode(term, "treat_race:black" = "Race x Black (H1)",
                            "treat_gender:woman" = "Gender x Woman (H2)"))
  })
  p <- ggplot(cd, aes(estimate, outcome, color = label, shape = label)) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                    position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "Treatment Effect Interactions on Trust Outcomes",
         subtitle = "OLS with HC2 robust SEs | 95% CI",
         x = "Coefficient (scale points)", y = NULL,
         color = NULL, shape = NULL) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
  ggsave(file.path(OUTPUT_DIR, "figures/coefficient_plot.png"),
         p, width = 8, height = 5, dpi = 150)
  list(rows = nrow(cd))
})

# ── 11. Validation summary ────────────────────────────────────────────────────

cat("\n================= PIPELINE VALIDATION =================\n")
n_pass <- sum(map_lgl(.stage_log, "ok"))
n_tot  <- length(.stage_log)
cat(sprintf("Stages passed: %d / %d\n", n_pass, n_tot))
if (SMOKE_TEST) {
  cat("Mode: SMOKE TEST (estimates suppressed; outputs in a scratch dir)\n")
  cat(sprintf("Optional pkgs -> modelsummary:%s  grf:%s  mediation:%s\n",
              HAS_MODELSUMMARY, HAS_GRF, HAS_MEDIATION))
}
failed <- names(.stage_log)[!map_lgl(.stage_log, "ok")]
if (length(failed)) cat("Failed stages:", paste(failed, collapse = ", "), "\n")
cat("======================================================\n")

# Non-zero exit if any *non-optional* stage failed
optional_stages <- c("HTE causal forest", "mediation (procedural justice)")
hard_fail <- failed[!failed %in% optional_stages]
if (length(hard_fail) > 0) quit(status = 1)
