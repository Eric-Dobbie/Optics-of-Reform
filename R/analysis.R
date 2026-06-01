# ==============================================================================
# analysis.R
# Pre-Analysis Plan — Final Analysis Script
# "The Optics of Reform: Diversification and Trust in Policing"
#
# Design: 3-arm RCT (Control / Race / Gender) | N ≈ 1,200 | Prolific
# Primary hypotheses:
#   H1: Black respondents report higher trust under Race Treatment vs. Control
#   H2: Women report higher trust under Gender Treatment vs. Control
# ==============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────

required_pkgs <- c(
  "tidyverse", "estimatr", "modelsummary", "grf",
  "mediation", "BradleyTerry2", "sandwich", "lmtest",
  "ggplot2", "dotwhisker", "scales", "knitr"
)
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(tidyverse)
library(estimatr)
library(modelsummary)
library(grf)
library(mediation)
library(BradleyTerry2)
library(ggplot2)
library(dotwhisker)

dir.create("R/output", showWarnings = FALSE)
dir.create("R/output/tables", showWarnings = FALSE)
dir.create("R/output/figures", showWarnings = FALSE)

# ── 2. Data Loading and Cleaning ──────────────────────────────────────────────

# Replace path with actual data file path when available.
# df <- read_csv("data/optics_of_reform_raw.csv")

# Synthetic data for script validation (100 rows, correct schema)
set.seed(42)
n_synth <- 100
df_raw <- tibble(
  # Identifiers
  ResponseId       = paste0("R", 1:n_synth),
  # Treatment (0 = Control, 1 = Race, 2 = Gender)
  treatment_raw    = sample(0:2, n_synth, replace = TRUE),
  # Demographics
  race             = sample(c("White","Black","Hispanic","Asian","Other"),
                            n_synth, replace = TRUE, prob = c(.5,.3,.1,.05,.05)),
  gender           = sample(c("Man","Woman","Nonbinary","Prefer not"),
                            n_synth, replace = TRUE, prob = c(.45,.48,.04,.03)),
  age              = sample(18:80, n_synth, replace = TRUE),
  party_id         = sample(c("Democrat","Republican","Independent","Lean Dem",
                              "Lean Rep","Neither"), n_synth, replace = TRUE),
  ideology         = sample(1:7, n_synth, replace = TRUE),
  education        = sample(1:6, n_synth, replace = TRUE),
  # Pre-treatment attitudes (7-pt Likert)
  att_diversity    = sample(1:7, n_synth, replace = TRUE),
  att_poc          = sample(1:7, n_synth, replace = TRUE),
  # Use-of-force contact
  uof_contact      = sample(c(0L, 1L), n_synth, replace = TRUE),
  # Primary outcomes (7-pt Likert)
  Trust_General    = sample(1:7, n_synth, replace = TRUE),
  Trust_Self       = sample(1:7, n_synth, replace = TRUE),
  Trust_Hiring     = sample(1:7, n_synth, replace = TRUE),
  # Secondary outcomes
  Procedural_Justice = sample(1:7, n_synth, replace = TRUE),
  Police_Resources   = sample(1:7, n_synth, replace = TRUE),
  # Open-ended
  Explain          = sample(c("More fair","Skeptical","No change",
                              "Promising","Does not matter"), n_synth, replace = TRUE),
  # Timing
  survey_seconds   = runif(n_synth, 60, 900)
)

# ── Variable construction ──────────────────────────────────────────────────────

df <- df_raw |>
  # Exclusion: < 120 seconds (inattentive)
  filter(survey_seconds >= 120) |>
  mutate(
    # Treatment indicators
    treat_race   = as.integer(treatment_raw == 1),
    treat_gender = as.integer(treatment_raw == 2),
    treat_any    = as.integer(treatment_raw > 0),
    treatment    = factor(treatment_raw, levels = 0:2,
                          labels = c("Control", "Race", "Gender")),

    # Subgroup flags
    black        = as.integer(race == "Black"),
    woman        = as.integer(gender == "Woman"),

    # Party collapse
    democrat     = as.integer(party_id %in% c("Democrat", "Lean Dem")),
    republican   = as.integer(party_id %in% c("Republican", "Lean Rep")),

    # Covariate: cynicism proxy (reverse att_poc)
    cynicism     = 8L - att_poc,

    # Block ID for FE robustness
    block        = paste(race, party_id, sep = "_")
  )

cat(sprintf("Analysis N (after exclusions): %d\n", nrow(df)))
cat(sprintf("Excluded (< 120s): %d\n", nrow(df_raw) - nrow(df)))

# ── 3. Balance Check (Table 1) ────────────────────────────────────────────────

balance_covariates <- c("age", "ideology", "education", "uof_contact",
                        "att_diversity", "att_poc", "democrat", "republican",
                        "black", "woman")

balance_models <- map(balance_covariates, ~ {
  fmla <- as.formula(paste(.x, "~ treatment"))
  lm(fmla, data = df)
})
names(balance_models) <- balance_covariates

modelsummary(
  balance_models,
  output = "R/output/tables/table1_balance.tex",
  stars = TRUE,
  title = "Covariate Balance Across Treatment Arms",
  gof_map = c("nobs", "r.squared"),
  fmt = 2
)

# Omnibus F-test: regress treatment dummy on all covariates
omnibus_fmla <- as.formula(
  paste("treat_race ~", paste(balance_covariates, collapse = " + "))
)
omnibus_race   <- lm(omnibus_fmla, data = df)
omnibus_gender <- update(omnibus_race, treat_gender ~ .)
cat("\n--- Omnibus Balance F-test ---\n")
cat("Race arm:   F =", round(summary(omnibus_race)$fstatistic[1], 3),
    " p =", round(pf(summary(omnibus_race)$fstatistic[1],
                     summary(omnibus_race)$fstatistic[2],
                     summary(omnibus_race)$fstatistic[3],
                     lower.tail = FALSE), 4), "\n")
cat("Gender arm: F =", round(summary(omnibus_gender)$fstatistic[1], 3),
    " p =", round(pf(summary(omnibus_gender)$fstatistic[1],
                     summary(omnibus_gender)$fstatistic[2],
                     summary(omnibus_gender)$fstatistic[3],
                     lower.tail = FALSE), 4), "\n")

# ── 4. Primary OLS (Table 2) ─────────────────────────────────────────────────
# Equation from PAP:
#   Trust = β0 + β1*D_race + β2*D_gender + β3*Black + β4*Woman
#         + β5*(Black × D_race) + β6*(Woman × D_gender) + γX + ε
# HC2 robust SEs via estimatr::lm_robust

covariates <- c("age", "ideology", "education", "uof_contact",
                "att_diversity", "att_poc", "democrat", "republican")
cov_str    <- paste(covariates, collapse = " + ")

outcomes_primary   <- c("Trust_General", "Trust_Self", "Trust_Hiring")
outcomes_secondary <- c("Procedural_Justice", "Police_Resources")
outcomes_all       <- c(outcomes_primary, outcomes_secondary)

fit_primary <- function(outcome, data) {
  fmla <- as.formula(paste0(
    outcome, " ~ treat_race + treat_gender + black + woman +
     treat_race:black + treat_gender:woman + ", cov_str
  ))
  lm_robust(fmla, data = data, se_type = "HC2")
}

models_primary <- map(outcomes_all, fit_primary, data = df)
names(models_primary) <- outcomes_all

modelsummary(
  models_primary[outcomes_primary],
  output = "R/output/tables/table2_primary.tex",
  coef_map = c(
    "treat_race"               = "Race Treatment",
    "treat_gender"             = "Gender Treatment",
    "black"                    = "Black respondent",
    "woman"                    = "Woman respondent",
    "treat_race:black"         = "Race × Black (H1)",
    "treat_gender:woman"       = "Gender × Woman (H2)"
  ),
  stars = TRUE,
  title  = "Effect of Diversity Framing on Trust in Policing (Primary Outcomes)",
  gof_map = c("nobs", "r.squared"),
  fmt = 3
)

cat("\nTable 2 (primary outcomes) saved.\n")

# ── 5. Secondary Outcomes ─────────────────────────────────────────────────────

modelsummary(
  models_primary[outcomes_secondary],
  output = "R/output/tables/table3_secondary.tex",
  coef_map = c(
    "treat_race"       = "Race Treatment",
    "treat_gender"     = "Gender Treatment",
    "treat_race:black" = "Race × Black",
    "treat_gender:woman" = "Gender × Woman"
  ),
  stars = TRUE,
  title  = "Effect of Diversity Framing on Procedural Justice and Resource Support",
  gof_map = c("nobs", "r.squared"),
  fmt = 3
)

# ── 6. Multiple Testing Correction (Holm) ─────────────────────────────────────

extract_pval <- function(model, coef_name) {
  tidy_m <- broom::tidy(model)
  tidy_m$p.value[tidy_m$term == coef_name]
}

# Primary family: H1 interaction across 3 primary outcomes; H2 interaction across 3 outcomes
p_h1 <- map_dbl(models_primary[outcomes_primary], extract_pval, "treat_race:black")
p_h2 <- map_dbl(models_primary[outcomes_primary], extract_pval, "treat_gender:woman")

p_primary_family <- c(p_h1, p_h2)
p_holm           <- p.adjust(p_primary_family, method = "holm")

holm_tbl <- tibble(
  Outcome     = rep(outcomes_primary, 2),
  Hypothesis  = c(rep("H1 (Race × Black)", 3), rep("H2 (Gender × Woman)", 3)),
  p_raw       = p_primary_family,
  p_holm      = p_holm,
  sig_holm    = p_holm < 0.05
)

cat("\n=== Holm-corrected p-values (primary family) ===\n")
print(holm_tbl)
write_csv(holm_tbl, "R/output/tables/holm_correction.csv")

# ── 7. Heterogeneous Treatment Effects (GRF) ─────────────────────────────────

hte_covariates <- c("age", "ideology", "education", "uof_contact",
                    "att_diversity", "att_poc", "democrat", "republican",
                    "black", "woman", "cynicism")

run_causal_forest <- function(outcome, treat_var, data, covars) {
  df_cf <- data |>
    select(all_of(c(outcome, treat_var, covars))) |>
    drop_na()

  Y <- df_cf[[outcome]]
  W <- df_cf[[treat_var]]
  X <- as.matrix(df_cf[covars])

  cf <- causal_forest(X = X, Y = Y, W = W, num.trees = 2000, seed = 42)

  # Best linear projection of CATEs on covariates
  blp <- best_linear_projection(cf, A = X[, c("black", "woman", "democrat",
                                               "uof_contact", "cynicism"),
                                          drop = FALSE])
  list(forest = cf, blp = blp)
}

# Race treatment HTE on Trust_General
hte_race_trust <- run_causal_forest(
  outcome   = "Trust_General",
  treat_var = "treat_race",
  data      = df,
  covars    = hte_covariates
)
cat("\n=== GRF Best Linear Projection: Race Treatment, Trust_General ===\n")
print(hte_race_trust$blp)

# Gender treatment HTE on Trust_General
hte_gender_trust <- run_causal_forest(
  outcome   = "Trust_General",
  treat_var = "treat_gender",
  data      = df,
  covars    = hte_covariates
)
cat("\n=== GRF Best Linear Projection: Gender Treatment, Trust_General ===\n")
print(hte_gender_trust$blp)

# CATE distribution plot (Race treatment)
cates_race <- predict(hte_race_trust$forest)$predictions
p_cate <- ggplot(tibble(CATE = cates_race), aes(x = CATE)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribution of Estimated CATEs",
       subtitle = "Race Treatment → Trust_General (Causal Forest)",
       x = "Estimated CATE (scale points)", y = "Count") +
  theme_bw(base_size = 13)

ggsave("R/output/figures/cate_distribution_race.png",
       p_cate, width = 7, height = 4, dpi = 150)

# ── 8. Mechanism Tests — Mediation Analysis ───────────────────────────────────

run_mediation <- function(treatment_var, mediator, outcome, data,
                          n_sims = 1000, seed = 42) {
  set.seed(seed)
  med_fmla    <- as.formula(paste(mediator, "~", treatment_var))
  outcome_fmla <- as.formula(paste(outcome, "~", treatment_var, "+", mediator))

  med_model <- lm(med_fmla, data = data)
  out_model <- lm(outcome_fmla, data = data)

  mediate(med_model, out_model,
          treat = treatment_var, mediator = mediator,
          robustSE = TRUE, sims = n_sims)
}

# Mechanism 1: Identity Activation via Procedural Justice
# H1: Race treatment → Procedural_Justice → Trust_General (Black respondents)
df_black <- filter(df, black == 1)

cat("\n--- Mediation: Race Treatment → Procedural Justice → Trust_General (Black) ---\n")
if (nrow(df_black) >= 20) {
  med_h1 <- run_mediation("treat_race", "Procedural_Justice",
                           "Trust_General", df_black)
  print(summary(med_h1))
  capture.output(summary(med_h1),
                 file = "R/output/tables/mediation_h1.txt")
} else {
  cat("  [Insufficient n for mediation — run with real data]\n")
}

# H2: Gender treatment → Procedural_Justice → Trust_General (Women)
df_women <- filter(df, woman == 1)
cat("\n--- Mediation: Gender Treatment → Procedural Justice → Trust_General (Women) ---\n")
if (nrow(df_women) >= 20) {
  med_h2 <- run_mediation("treat_gender", "Procedural_Justice",
                           "Trust_General", df_women)
  print(summary(med_h2))
  capture.output(summary(med_h2),
                 file = "R/output/tables/mediation_h2.txt")
} else {
  cat("  [Insufficient n for mediation — run with real data]\n")
}

# ── 9. Open-Ended Response Analysis (Bradley-Terry) ──────────────────────────

# Placeholder: expects a pairwise comparison dataset produced by LLM scoring.
# Each row represents one pairwise comparison between two respondents'
# open-ended Explain responses, with a winner column indicating which was
# judged more trusting.
#
# bt_data should have columns: player1, player2, win1
# (win1 = 1 if player1's response showed more trust, 0 otherwise)
#
# To generate bt_data from LLM output, see companion script R/bt_prep.R (TBD)

run_bradley_terry <- function(bt_data) {
  BTm(
    outcome  = cbind(win1, 1 - win1),
    player1  = player1,
    player2  = player2,
    data     = bt_data
  )
}

# Stub: create synthetic pairwise comparison data
players     <- paste0("P", 1:nrow(df))
n_pairs     <- 200
bt_synth    <- tibble(
  player1 = sample(players, n_pairs, replace = TRUE),
  player2 = sample(players, n_pairs, replace = TRUE),
  win1    = rbinom(n_pairs, 1, 0.5)
) |> filter(player1 != player2)

cat("\n--- Bradley-Terry Model (synthetic data) ---\n")
if (nrow(bt_synth) > 0) {
  tryCatch({
    bt_model <- BTm(
      outcome = cbind(win1, 1 - win1),
      player1 = factor(player1, levels = players),
      player2 = factor(player2, levels = players),
      data    = bt_synth
    )
    cat("  Bradley-Terry model fitted on synthetic pairwise data.\n")
    # Merge latent trust scores back to df for treatment effect estimation
    # bt_scores <- BTabilities(bt_model)
    # df <- df |> left_join(as_tibble(bt_scores, rownames = "ResponseId"), by = "ResponseId")
    # lm_robust(ability ~ treat_race + treat_gender, data = df)
  }, error = function(e) {
    cat("  BT model error (expected with random synthetic data):", e$message, "\n")
  })
}

# ── 10. Robustness Checks ─────────────────────────────────────────────────────

# 10a. Block fixed effects
fit_block_fe <- function(outcome) {
  fmla <- as.formula(paste0(
    outcome, " ~ treat_race + treat_gender + black + woman +
     treat_race:black + treat_gender:woman + ",
    cov_str, " + block"
  ))
  lm_robust(fmla, data = df, se_type = "HC2")
}

models_fe <- map(outcomes_primary, fit_block_fe)
names(models_fe) <- outcomes_primary

modelsummary(
  models_fe,
  output = "R/output/tables/table_robustness_fe.tex",
  coef_map = c(
    "treat_race"       = "Race Treatment",
    "treat_gender"     = "Gender Treatment",
    "treat_race:black" = "Race × Black (H1)",
    "treat_gender:woman" = "Gender × Woman (H2)"
  ),
  stars = TRUE,
  title  = "Primary Results with Block Fixed Effects (Robustness)",
  gof_map = c("nobs", "r.squared"),
  fmt = 3
)

# 10b. Full sample (no exclusions)
fit_full <- function(outcome) {
  fmla <- as.formula(paste0(
    outcome, " ~ treat_race + treat_gender + black + woman +
     treat_race:black + treat_gender:woman + ", cov_str
  ))
  lm_robust(fmla, data = df_raw |>
              mutate(
                treat_race   = as.integer(treatment_raw == 1),
                treat_gender = as.integer(treatment_raw == 2),
                black        = as.integer(race == "Black"),
                woman        = as.integer(gender == "Woman"),
                democrat     = as.integer(party_id %in% c("Democrat","Lean Dem")),
                republican   = as.integer(party_id %in% c("Republican","Lean Rep")),
                cynicism     = 8L - att_poc
              ),
            se_type = "HC2")
}

models_full <- map(outcomes_primary, fit_full)
names(models_full) <- outcomes_primary

modelsummary(
  models_full,
  output = "R/output/tables/table_robustness_fullsample.tex",
  coef_map = c(
    "treat_race"         = "Race Treatment",
    "treat_gender"       = "Gender Treatment",
    "treat_race:black"   = "Race × Black (H1)",
    "treat_gender:woman" = "Gender × Woman (H2)"
  ),
  stars = TRUE,
  title  = "Primary Results Including Excluded Respondents (Robustness)",
  gof_map = c("nobs", "r.squared"),
  fmt = 3
)

cat("\nRobustness tables saved.\n")

# ── 11. Coefficient Plot (Figure 2) ───────────────────────────────────────────

# Extract coefficients for H1 and H2 interactions across primary outcomes
coef_plot_data <- map_dfr(outcomes_primary, function(y) {
  m   <- models_primary[[y]]
  est <- broom::tidy(m) |>
    filter(term %in% c("treat_race:black", "treat_gender:woman")) |>
    mutate(
      outcome  = y,
      label    = recode(term,
                        "treat_race:black"   = "Race × Black (H1)",
                        "treat_gender:woman" = "Gender × Woman (H2)")
    )
  est
})

p_coef <- ggplot(coef_plot_data,
                 aes(x = estimate, y = outcome,
                     color = label, shape = label)) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = 0.5),
                  size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "Treatment Effect Interactions on Primary Trust Outcomes",
    subtitle = "OLS with HC2 robust SEs | 95% confidence intervals",
    x        = "Coefficient estimate (scale points)",
    y        = "Outcome",
    color    = "Interaction term",
    shape    = "Interaction term"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("R/output/figures/coefficient_plot.png",
       p_coef, width = 8, height = 5, dpi = 150)
cat("Coefficient plot saved.\n")

cat("\n=== Analysis complete. All outputs in R/output/ ===\n")
