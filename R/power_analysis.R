# ==============================================================================
# power_analysis.R
# Minimum Detectable Effect (MDE) and power analysis
# "The Optics of Reform: Diversification and Trust in Policing"
#
# Design: 3-arm RCT (Control / Race Treatment / Gender Treatment)
#   N = 1,200 | Block-randomized by race x party
#   Primary outcomes: 7-point Likert trust items
#   alpha = 0.05, power = 0.80
# ==============================================================================

# ── 1. Setup ──────────────────────────────────────────────────────────────────

required_pkgs <- c("pwr", "tidyverse", "knitr", "kableExtra")
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(pwr)
library(tidyverse)

# ── 2. Helpers ────────────────────────────────────────────────────────────────

#' Compute two-sample MDE (scale points) given n per arm, SD, alpha, power
mde_two_sample <- function(n_per_arm, sd = 1.5, alpha = 0.05, power = 0.80) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(power)
  (z_alpha + z_beta) * sd * sqrt(2 / n_per_arm)
}

#' Compute required n per arm given target MDE, SD, alpha, power
n_required <- function(mde, sd = 1.5, alpha = 0.05, power = 0.80) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(power)
  ceiling(2 * ((z_alpha + z_beta) * sd / mde)^2)
}

# ── 3. Main Effect MDE — Pre-Pilot Estimates ──────────────────────────────────

# Scenario: one treatment arm vs. control, 400 per arm, vary SD
n_per_arm_main <- 400

sd_assumptions <- c(1.0, 1.5, 2.0)
mde_main <- tibble(
  SD    = sd_assumptions,
  MDE   = map_dbl(SD, ~ mde_two_sample(n_per_arm_main, sd = .x)),
  d     = MDE / SD,   # standardized effect size
  N_arm = n_per_arm_main,
  N_total = n_per_arm_main * 3
)

cat("=== Main Effect MDE (400 per arm, alpha=0.05, power=0.80) ===\n")
print(mde_main)

# ── 4. Subgroup MDEs ──────────────────────────────────────────────────────────

# H1: Black respondents (≥30% oversample → ~120 Black respondents per relevant arm)
n_black_per_arm <- 120
# H2: Women respondents (~50% of sample → ~200 women per relevant arm)
n_women_per_arm <- 200

subgroup_mde <- tibble(
  Subgroup     = c("Black respondents (H1)", "Women respondents (H2)"),
  n_per_arm    = c(n_black_per_arm, n_women_per_arm),
  SD           = 1.5,
  MDE          = map_dbl(n_per_arm, ~ mde_two_sample(.x, sd = 1.5)),
  d            = MDE / SD
)

cat("\n=== Subgroup MDEs (SD=1.5, alpha=0.05, power=0.80) ===\n")
print(subgroup_mde)

# ── 5. Sensitivity Table — MDE Across N and SD ───────────────────────────────

sample_sizes <- c(600, 800, 1000, 1200, 1500)
sd_vals      <- c(1.0, 1.5, 2.0)

sensitivity_tbl <- expand_grid(
  N_total = sample_sizes,
  SD      = sd_vals
) |>
  mutate(
    n_per_arm = N_total / 3,
    MDE       = map2_dbl(n_per_arm, SD, ~ mde_two_sample(.x, sd = .y)),
    d         = MDE / SD
  )

cat("\n=== MDE Sensitivity Table ===\n")
sensitivity_wide <- sensitivity_tbl |>
  mutate(label = paste0("SD=", SD)) |>
  select(N_total, n_per_arm, label, MDE) |>
  pivot_wider(names_from = label, values_from = MDE)
print(sensitivity_wide)

# ── 6. Power Curves ───────────────────────────────────────────────────────────

n_seq  <- seq(50, 600, by = 10)   # n per arm
effect_sizes <- c(0.20, 0.30, 0.40)  # standardized Cohen's d

power_curve_data <- expand_grid(
  n_per_arm = n_seq,
  d         = effect_sizes,
  SD        = 1.5
) |>
  mutate(
    power = map2_dbl(d, n_per_arm, ~ {
      pwr.t.test(n = .y, d = .x, sig.level = 0.05, type = "two.sample")$power
    })
  )

p_power_curve <- ggplot(power_curve_data, aes(x = n_per_arm, y = power,
                                               color = factor(d), linetype = factor(d))) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 400, linetype = "dotted", color = "grey60") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2",
                     labels = c("d = 0.20 (small)", "d = 0.30", "d = 0.40")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"),
                        labels = c("d = 0.20 (small)", "d = 0.30", "d = 0.40")) +
  labs(
    title    = "Power Curves by Effect Size and Sample Size (per arm)",
    subtitle = "Two-sided t-test, α = 0.05 | SD = 1.5 | Vertical line: planned N = 400/arm",
    x        = "N per arm",
    y        = "Statistical power",
    color    = "Effect size",
    linetype = "Effect size"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("R/output/power_curve.png", p_power_curve, width = 8, height = 5, dpi = 150)
cat("\nPower curve saved to R/output/power_curve.png\n")

# ── 7. MDE-by-N Plot ─────────────────────────────────────────────────────────

mde_by_n <- expand_grid(
  N_total = seq(300, 2000, by = 50),
  SD      = c(1.0, 1.5, 2.0)
) |>
  mutate(
    n_per_arm = N_total / 3,
    MDE       = map2_dbl(n_per_arm, SD, ~ mde_two_sample(.x, sd = .y)),
    SD_label  = paste0("SD = ", SD)
  )

p_mde_n <- ggplot(mde_by_n, aes(x = N_total, y = MDE,
                                  color = SD_label, linetype = SD_label)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 1200, linetype = "dotted", color = "grey50") +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title    = "Minimum Detectable Effect by Total Sample Size",
    subtitle = "Equal allocation across 3 arms | α = 0.05, power = 0.80 | Vertical: N = 1,200",
    x        = "Total N",
    y        = "MDE (scale points)",
    color    = "Outcome SD",
    linetype = "Outcome SD"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("R/output/mde_by_n.png", p_mde_n, width = 8, height = 5, dpi = 150)
cat("MDE-by-N plot saved to R/output/mde_by_n.png\n")

# ── 8. Post-Pilot MDE Update (Placeholder) ───────────────────────────────────

# After the pilot study is run, replace pilot_sd with the observed SD from pilot data.
# pilot_sd <- sd(pilot_data$Trust_General, na.rm = TRUE)
#
# post_pilot_mde <- mde_two_sample(n_per_arm = 400, sd = pilot_sd)
# cat("\nPost-pilot MDE (Trust_General):", round(post_pilot_mde, 3), "\n")
#
# If post_pilot_mde exceeds a meaningful threshold, recalculate required N:
# n_adj <- n_required(mde = post_pilot_mde * 0.9, sd = pilot_sd)
# cat("Adjusted N per arm for same power:", n_adj, "\n")

cat("\n[Post-pilot section is a placeholder — fill in pilot_sd after piloting.]\n")

# ── 9. Summary Output ─────────────────────────────────────────────────────────

cat("\n=== Summary ===\n")
cat(sprintf(
  "Planned N = 1,200 (400 per arm)\n  Main effect MDE (SD=1.5): %.3f scale points (d = %.2f)\n",
  mde_two_sample(400, 1.5), mde_two_sample(400, 1.5) / 1.5
))
cat(sprintf(
  "  H1 subgroup MDE (n=120, SD=1.5): %.3f scale points (d = %.2f)\n",
  mde_two_sample(120, 1.5), mde_two_sample(120, 1.5) / 1.5
))
cat(sprintf(
  "  H2 subgroup MDE (n=200, SD=1.5): %.3f scale points (d = %.2f)\n",
  mde_two_sample(200, 1.5), mde_two_sample(200, 1.5) / 1.5
))
