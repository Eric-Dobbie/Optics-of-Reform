# ==============================================================================
# pilot_power_analysis.R
# Post-pilot MDE / power update using observed pilot data
# "The Optics of Reform: Diversification and Trust in Policing"
#
# Reads the Qualtrics pilot export, maps 7-point Likert responses to a numeric
# scale, estimates the empirical SD of each primary outcome, and recomputes the
# minimum detectable effect (MDE) and required sample size for the full study.
#
# Companion to R/power_analysis.R (pre-pilot assumptions).
# ==============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────

required_pkgs <- c("tidyverse")
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages(library(tidyverse))

dir.create("R/output", showWarnings = FALSE)
dir.create("R/output/figures", showWarnings = FALSE)
dir.create("R/output/tables", showWarnings = FALSE)

# ── 2. MDE helpers (shared with R/power_analysis.R) ───────────────────────────

mde_two_sample <- function(n_per_arm, sd, alpha = 0.05, power = 0.80) {
  (qnorm(1 - alpha / 2) + qnorm(power)) * sd * sqrt(2 / n_per_arm)
}

n_required <- function(mde, sd, alpha = 0.05, power = 0.80) {
  ceiling(2 * ((qnorm(1 - alpha / 2) + qnorm(power)) * sd / mde)^2)
}

# ── 3. Load pilot data ────────────────────────────────────────────────────────
# Qualtrics export: row 1 = variable names, rows 2-3 = question text / import IDs

pilot_path <- "Data/Pilot/Trust_Policing_Pilot_June 17, 2026_17.06.csv"

pilot_raw <- read_csv(pilot_path, skip = 0, show_col_types = FALSE)
# Drop the two Qualtrics metadata header rows (kept as data rows 1-2)
pilot_raw <- pilot_raw[-c(1, 2), ]

cat(sprintf("Raw pilot rows (excl. metadata headers): %d\n", nrow(pilot_raw)))

# ── 4. Recode 7-point Likert text → numeric ───────────────────────────────────

likert_map <- c(
  "Strongly disagree"          = 1,
  "Disagree"                   = 2,
  "Somewhat disagree"          = 3,
  "Neither agree nor disagree" = 4,
  "Somewhat agree"             = 5,
  "Agree"                      = 6,
  "Strongly agree"             = 7
)

recode_likert <- function(x) unname(likert_map[trimws(x)])

primary_outcomes   <- c("Trust_General", "Trust_Self", "Trust_Hiring")
secondary_outcomes <- c("Procedural_Justice", "Police_resources")
likert_vars        <- c(primary_outcomes, secondary_outcomes,
                        "Women in police", "Minority_police")

pilot <- pilot_raw |>
  mutate(across(all_of(likert_vars), recode_likert)) |>
  mutate(
    arm        = na_if(Treatment, ""),
    duration   = as.numeric(`Duration (in seconds)`),
    black      = as.integer(Race == "Black or African-American"),
    woman      = as.integer(Gender == "Female")
  ) |>
  # Analysis sample: assigned to an arm and provided the primary outcome
  filter(!is.na(arm), !is.na(Trust_General))

cat(sprintf("Analysis sample (assigned arm + primary outcome): %d\n", nrow(pilot)))
cat("\nArm sizes:\n"); print(table(pilot$arm))
cat(sprintf("\nObserved %% Black: %.1f%% | %% Women: %.1f%%\n",
            100 * mean(pilot$black, na.rm = TRUE),
            100 * mean(pilot$woman, na.rm = TRUE)))

# ── 5. Empirical SD of outcomes (pooled and by arm) ──────────────────────────

sd_pooled <- pilot |>
  summarise(across(all_of(c(primary_outcomes, secondary_outcomes)),
                   ~ sd(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "outcome", values_to = "sd_pooled")

mean_pooled <- pilot |>
  summarise(across(all_of(c(primary_outcomes, secondary_outcomes)),
                   ~ mean(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "outcome", values_to = "mean")

# Within-arm (residual) SD — the relevant quantity for treatment-effect power
resid_sd <- map_dfr(c(primary_outcomes, secondary_outcomes), function(y) {
  fit <- lm(reformulate("arm", response = y), data = pilot)
  tibble(outcome = y, sd_within = sd(residuals(fit)))
})

sd_summary <- mean_pooled |>
  left_join(sd_pooled, by = "outcome") |>
  left_join(resid_sd, by = "outcome")

cat("\n=== Observed outcome means and SDs (pilot) ===\n")
print(as.data.frame(sd_summary), digits = 3)

# Use within-arm SD of primary outcomes as the power input (conservative: max)
sd_primary_within <- resid_sd |>
  filter(outcome %in% primary_outcomes) |>
  pull(sd_within)
sigma_hat     <- max(sd_primary_within)   # conservative
sigma_hat_avg <- mean(sd_primary_within)  # average across primary outcomes

cat(sprintf(
  "\nPrimary-outcome within-arm SD: avg = %.3f, max (conservative) = %.3f\n",
  sigma_hat_avg, sigma_hat
))

# ── 6. Updated MDE for the full study (N = 1,200, 3 arms) ─────────────────────

n_per_arm_full <- 400

updated_mde <- tibble(
  basis = c("Pre-pilot assumption (SD=1.5)",
            "Pilot avg within-arm SD",
            "Pilot max within-arm SD (conservative)"),
  sd    = c(1.5, sigma_hat_avg, sigma_hat)
) |>
  mutate(
    MDE_main = mde_two_sample(n_per_arm_full, sd),
    d_main   = MDE_main / sd
  )

cat("\n=== Updated MAIN-EFFECT MDE (400 per arm) ===\n")
print(as.data.frame(updated_mde), digits = 3)

# ── 7. Updated subgroup MDE using observed subgroup proportions ──────────────

p_black <- mean(pilot$black, na.rm = TRUE)
p_woman <- mean(pilot$woman, na.rm = TRUE)

# Expected subgroup n per arm in full study (N=1,200, 400/arm),
# using observed proportions. NOTE: full study oversamples Black to >=30%.
n_black_arm_observed <- round(400 * p_black)
n_black_arm_target   <- round(400 * 0.30)   # design oversample floor
n_woman_arm_observed <- round(400 * p_woman)

subgroup_mde <- tibble(
  subgroup  = c("Black (observed prop.)", "Black (30% oversample target)",
                "Women (observed prop.)"),
  n_per_arm = c(n_black_arm_observed, n_black_arm_target, n_woman_arm_observed),
  sd        = sigma_hat
) |>
  mutate(
    MDE = mde_two_sample(n_per_arm, sd),
    d   = MDE / sd
  )

cat("\n=== Updated SUBGROUP MDE (conservative SD) ===\n")
print(as.data.frame(subgroup_mde), digits = 3)

# ── 8. Required N to detect a target effect ──────────────────────────────────

target_effects <- c(0.25, 0.30, 0.40, 0.50)  # in scale points
required_n_tbl <- tibble(
  target_MDE_points = target_effects,
  sd                = sigma_hat
) |>
  mutate(
    n_per_arm   = n_required(target_MDE_points, sd),
    N_total_3arm = n_per_arm * 3
  )

cat("\n=== Required N per arm to detect target effects (conservative SD) ===\n")
print(as.data.frame(required_n_tbl), digits = 3)

write_csv(updated_mde,    "R/output/tables/pilot_updated_mde.csv")
write_csv(subgroup_mde,   "R/output/tables/pilot_subgroup_mde.csv")
write_csv(required_n_tbl, "R/output/tables/pilot_required_n.csv")
write_csv(sd_summary,     "R/output/tables/pilot_outcome_sds.csv")

# ── 9. Figure: MDE by N at the pilot-estimated SD ────────────────────────────

mde_curve <- tibble(N_total = seq(300, 2000, by = 25)) |>
  mutate(
    n_per_arm = N_total / 3,
    MDE_avg   = mde_two_sample(n_per_arm, sigma_hat_avg),
    MDE_max   = mde_two_sample(n_per_arm, sigma_hat)
  ) |>
  pivot_longer(c(MDE_avg, MDE_max),
               names_to = "sd_basis", values_to = "MDE") |>
  mutate(sd_basis = recode(sd_basis,
                           MDE_avg = sprintf("Avg SD = %.2f", sigma_hat_avg),
                           MDE_max = sprintf("Max SD = %.2f", sigma_hat)))

p_pilot_mde <- ggplot(mde_curve, aes(N_total, MDE,
                                      color = sd_basis, linetype = sd_basis)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 1200, linetype = "dotted", color = "grey50") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title    = "Updated MDE by Sample Size (Pilot-Estimated SD)",
    subtitle = "Equal 3-arm allocation | alpha = 0.05, power = 0.80 | Dotted line: N = 1,200",
    x = "Total N", y = "MDE (scale points)",
    color = "Within-arm SD", linetype = "Within-arm SD"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("R/output/figures/pilot_mde_by_n.png", p_pilot_mde,
       width = 8, height = 5, dpi = 150)

# ── 10. Summary ───────────────────────────────────────────────────────────────

cat("\n================= PILOT POWER SUMMARY =================\n")
cat(sprintf("Pilot analysis N: %d (Control/Race/Gender)\n", nrow(pilot)))
cat(sprintf("Primary outcome within-arm SD: %.2f-%.2f (avg %.2f)\n",
            min(sd_primary_within), max(sd_primary_within), sigma_hat_avg))
cat(sprintf("At planned N=1,200 (400/arm), main-effect MDE = %.2f-%.2f points (d=%.2f-%.2f)\n",
            mde_two_sample(400, sigma_hat_avg), mde_two_sample(400, sigma_hat),
            mde_two_sample(400, sigma_hat_avg) / sigma_hat_avg,
            mde_two_sample(400, sigma_hat) / sigma_hat))
cat(sprintf("Observed share Black = %.0f%%, Women = %.0f%%\n",
            100 * p_black, 100 * p_woman))
cat("Outputs written to R/output/tables/ and R/output/figures/\n")
cat("=======================================================\n")
