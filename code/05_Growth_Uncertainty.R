# 05_Growth_Uncertainty.R
# Explore uncertainty in the dW and dFL growth models from 02_Growth.R.
#
# Covers:
#   1. Fixed-effect confidence intervals (Wald)
#   2. Variance partitioning: stream RE vs residual (ICC)
#   3. Lognormal retransformation correction: exp(sigma^2 / 2)
#   4. Negative-growth exclusion: how many fish dropped, directional bias
#   5. dW vs dFL consistency: compare dW model predictions to
#      L-W-derived dW from the dFL model
#   6. Parametric bootstrap of mean population growth (propagating
#      fixed-effect and RE uncertainty into the production rate)

library(tidyverse)
library(lme4)
library(patchwork)

# ============================================================
# Load saved results
# ============================================================

grw <- readRDS("growth_results.RDS")

growth_mod_FL <- grw$growth_mod_FL
growth_mod_W  <- grw$growth_mod_W
fe_W          <- grw$fe_W
sigma_W       <- grw$sigma_W
fe_FL         <- grw$fe_FL
sigma_FL      <- grw$sigma_FL
stream_re_W   <- grw$stream_re_W
stream_re_FL  <- grw$stream_re_FL
lw_a          <- grw$lw_a
lw_b          <- grw$lw_b
growth_dat    <- grw$growth_dat    # all matched pairs (includes neg growth)
growth_dat_W  <- grw$growth_dat_W  # subset used to fit dW model
growth_dat_FL <- grw$growth_dat_FL # subset used to fit dFL model

# ============================================================
# 1. Fixed-effect confidence intervals (Wald)
# ============================================================

cat("=== 1. Fixed-effect 95% CIs (Wald) ===\n\n")

ci_FL <- confint(growth_mod_FL, method = "Wald", parm = "beta_")
ci_W  <- confint(growth_mod_W,  method = "Wald", parm = "beta_")

cat("--- dFL model: log(dFL_day) ~ log(FL_aug) + (1|Location) ---\n")
print(round(cbind(estimate = fe_FL, ci_FL), 4))

cat("\n--- dW model: log(dW_day) ~ log(W_aug) + (1|Location) ---\n")
print(round(cbind(estimate = fe_W, ci_W), 4))

# ============================================================
# 2. Variance partitioning: ICC
# ============================================================

cat("\n=== 2. Variance partitioning ===\n\n")

vc_FL <- as.data.frame(VarCorr(growth_mod_FL))
vc_W  <- as.data.frame(VarCorr(growth_mod_W))

stream_var_FL  <- vc_FL$vcov[vc_FL$grp == "Location"]
resid_var_FL   <- vc_FL$vcov[vc_FL$grp == "Residual"]
icc_FL         <- stream_var_FL / (stream_var_FL + resid_var_FL)

stream_var_W   <- vc_W$vcov[vc_W$grp == "Location"]
resid_var_W    <- vc_W$vcov[vc_W$grp == "Residual"]
icc_W          <- stream_var_W / (stream_var_W + resid_var_W)

cat("dFL model:\n")
cat(sprintf("  Stream variance (log scale):   %.4f  (SD = %.4f)\n",
            stream_var_FL, sqrt(stream_var_FL)))
cat(sprintf("  Residual variance (log scale): %.4f  (SD = %.4f = sigma_FL)\n",
            resid_var_FL, sqrt(resid_var_FL)))
cat(sprintf("  ICC (stream / total):          %.3f\n\n", icc_FL))

cat("dW model:\n")
cat(sprintf("  Stream variance (log scale):   %.4f  (SD = %.4f)\n",
            stream_var_W, sqrt(stream_var_W)))
cat(sprintf("  Residual variance (log scale): %.4f  (SD = %.4f = sigma_W)\n",
            resid_var_W, sqrt(resid_var_W)))
cat(sprintf("  ICC (stream / total):          %.3f\n\n", icc_W))

# ============================================================
# 3. Lognormal retransformation correction
# ============================================================

cat("=== 3. Lognormal retransformation correction ===\n\n")
cat("The expected value of exp(X) when X ~ N(mu, sigma^2) is exp(mu + sigma^2/2).\n")
cat("The correction factor exp(sigma^2/2) inflates the naive back-transform.\n\n")

corr_FL <- exp(sigma_FL^2 / 2)
corr_W  <- exp(sigma_W^2  / 2)

cat(sprintf("dFL model: sigma_FL = %.4f  →  exp(sigma^2/2) = %.4f  (%.1f%% inflation)\n",
            sigma_FL, corr_FL, (corr_FL - 1) * 100))
cat(sprintf("dW  model: sigma_W  = %.4f  →  exp(sigma^2/2) = %.4f  (%.1f%% inflation)\n\n",
            sigma_W,  corr_W,  (corr_W  - 1) * 100))

cat("Interpretation: using exp(mu) alone would underestimate mean growth by these factors.\n\n")

# ============================================================
# 4. Negative-growth exclusion
# ============================================================

cat("=== 4. Negative-growth exclusion ===\n\n")

# All matched pairs (includes negative growers)
growth_all_W <- growth_dat %>%
  filter(!is.na(dW), W_aug > 0, W_oct > 0)

n_total_W    <- nrow(growth_all_W)
n_neg_W      <- sum(growth_all_W$dW_day <= 0)
n_pos_W      <- sum(growth_all_W$dW_day > 0)
pct_neg_W    <- 100 * n_neg_W / n_total_W

mean_all_W   <- mean(growth_all_W$dW_day)
mean_pos_W   <- mean(growth_all_W$dW_day[growth_all_W$dW_day > 0])
bias_pct     <- 100 * (mean_pos_W - mean_all_W) / abs(mean_all_W)

cat(sprintf("Fish with both W_aug and W_oct measured: %d\n", n_total_W))
cat(sprintf("  Positive dW_day (included in model):   %d\n", n_pos_W))
cat(sprintf("  Zero/negative dW_day (excluded):       %d  (%.1f%%)\n\n",
            n_neg_W, pct_neg_W))

cat(sprintf("Mean dW_day across ALL fish (raw):       %.5f g day-1\n", mean_all_W))
cat(sprintf("Mean dW_day for POSITIVE growers only:   %.5f g day-1\n", mean_pos_W))
cat(sprintf("Upward bias from exclusion:              %.1f%%\n\n", bias_pct))

cat("Note: the log-normal correction (section 3) partially compensates,\n")
cat("but the exclusion bias and the retransformation correction are separate issues.\n\n")

# By stream
cat("--- Negative growth by stream ---\n")
growth_all_W %>%
  group_by(Location) %>%
  summarise(
    n_total    = n(),
    n_neg      = sum(dW_day <= 0),
    pct_neg    = round(100 * n_neg / n_total, 1),
    mean_all   = round(mean(dW_day), 5),
    mean_pos   = round(mean(dW_day[dW_day > 0]), 5),
    .groups    = "drop"
  ) %>%
  arrange(desc(pct_neg)) %>%
  print()

# ============================================================
# 5. dW vs dFL consistency via L-W
# ============================================================

cat("\n=== 5. dW vs dFL consistency ===\n\n")
cat("If W = exp(lw_a) * FL^lw_b, then by chain rule:\n")
cat("  dW/dt ≈ exp(lw_a) * lw_b * FL^(lw_b - 1) * dFL/dt\n")
cat("This converts dFL model predictions to g day-1 for comparison.\n\n")

# Use fish in growth_dat_W (have both FL and W measurements)
check_df <- growth_dat_W %>%
  mutate(
    # Predicted log(dW_day) from dW model (using stream RE)
    re_W_val    = stream_re_W[Location],
    re_W_val    = ifelse(is.na(re_W_val), 0, re_W_val),
    mu_dW       = fe_W["(Intercept)"] + re_W_val + fe_W["logW_aug"] * log(W_aug),
    pred_dW_mod = exp(mu_dW + sigma_W^2 / 2),

    # Predicted log(dFL_day) from dFL model, then convert via L-W chain rule
    re_FL_val   = stream_re_FL[Location],
    re_FL_val   = ifelse(is.na(re_FL_val), 0, re_FL_val),
    mu_dFL      = fe_FL["(Intercept)"] + re_FL_val + fe_FL["logFL_aug"] * log(FL_aug),
    pred_dFL    = exp(mu_dFL + sigma_FL^2 / 2),   # mm day-1 (lognormal mean)
    pred_dW_LW  = exp(lw_a) * lw_b * FL_aug^(lw_b - 1) * pred_dFL,  # g day-1

    # Observed dW_day for reference
    obs_dW      = dW_day
  )

cat("--- Per-stream mean predicted dW (g day-1): dW model vs dFL-via-L-W ---\n")
check_df %>%
  group_by(Location) %>%
  summarise(
    n             = n(),
    obs_mean      = round(mean(obs_dW),        5),
    dW_model_mean = round(mean(pred_dW_mod),   5),
    dFL_LW_mean   = round(mean(pred_dW_LW),    5),
    ratio_dW_LW   = round(mean(pred_dW_mod) / mean(pred_dW_LW), 3),
    .groups = "drop"
  ) %>%
  print()

# Scatter: dW model vs dFL-via-L-W at fish level
p_consist <- ggplot(check_df, aes(x = pred_dW_LW, y = pred_dW_mod,
                                   colour = Location)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_log10() + scale_y_log10() +
  labs(x = expression("Predicted dW via dFL model + L-W  (g day"^{-1}*")"),
       y = expression("Predicted dW from dW model  (g day"^{-1}*")"),
       title    = "Consistency: dW model vs dFL model + L-W conversion",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; dashed = 1:1; both use lognormal mean correction",
       colour   = NULL) +
  theme_bw(base_size = 12)

print(p_consist)
ggsave("plot_dW_vs_dFL_consistency.png", p_consist, width = 6, height = 5, dpi = 300)

# ============================================================
# 6. Parametric bootstrap of mean population growth
# ============================================================
# Propagates uncertainty in fixed effects and residual sigma.
# (Stream REs treated as fixed — limited recaptures make full
#  Bayesian posterior impractical here.)

cat("\n=== 6. Parametric bootstrap: uncertainty in mean growth rate ===\n\n")
cat("Samples B realisations of (beta_0, beta_W) from their approximate\n")
cat("joint normal distribution (Wald), then computes the population-level\n")
cat("mean daily growth across the dep1 fish for each sample.\n\n")

# Load dep1 fish weights (via bioenergetics results for convenience)
bio       <- readRDS("bioenergetics_results.RDS")
fish_dep1 <- bio$fish_bioen

B     <- 5000
set.seed(42)

# Approximate covariance matrix of fixed effects (Wald)
vcov_W <- vcov(growth_mod_W)

# Draw B realisations of (Intercept, logW_aug) from joint normal
beta_draws <- MASS::mvrnorm(B, mu = fe_W, Sigma = vcov_W)

streams_all <- sort(unique(fish_dep1$Stream))

boot_means <- lapply(streams_all, function(s) {
  w_s  <- fish_dep1$W_aug[fish_dep1$Stream == s &
                           !is.na(fish_dep1$W_aug) & fish_dep1$W_aug > 0]
  if (length(w_s) == 0) return(NULL)
  re_s <- if (!is.na(stream_re_W[s])) stream_re_W[s] else 0

  # For each bootstrap draw, compute mean dW_day across dep1 fish
  means_b <- apply(beta_draws, 1, function(b) {
    mu  <- b[1] + re_s + b[2] * log(w_s)
    mean(exp(mu + sigma_W^2 / 2))   # lognormal mean, averaged over fish
  })

  data.frame(
    Stream  = s,
    boot_id = seq_len(B),
    dW_mean = means_b
  )
}) %>% bind_rows()

boot_summ <- boot_means %>%
  group_by(Stream) %>%
  summarise(
    point_est = mean(exp(fe_W["(Intercept)"] +
                         (if (!is.na(stream_re_W[Stream[1]])) stream_re_W[Stream[1]] else 0) +
                         fe_W["logW_aug"] * log(
                           fish_dep1$W_aug[fish_dep1$Stream == Stream[1] &
                                           !is.na(fish_dep1$W_aug) & fish_dep1$W_aug > 0])
                         + sigma_W^2 / 2)),
    boot_med  = round(median(dW_mean), 5),
    boot_lo95 = round(quantile(dW_mean, 0.025), 5),
    boot_hi95 = round(quantile(dW_mean, 0.975), 5),
    cv_pct    = round(100 * sd(dW_mean) / mean(dW_mean), 1),
    .groups   = "drop"
  )

cat("--- Bootstrap CI on mean dW_day (g day-1) per stream ---\n")
cat("(B = 5000; fixed effects sampled from Wald joint normal; RE treated as fixed)\n\n")
print(boot_summ %>% select(Stream, boot_med, boot_lo95, boot_hi95, cv_pct) %>%
        rename(`median` = boot_med, `lo 95%` = boot_lo95,
               `hi 95%` = boot_hi95, `CV (%)` = cv_pct))

# Plot bootstrap distributions
p_boot <- ggplot(boot_means, aes(x = dW_mean, y = reorder(Stream, dW_mean, median))) +
  geom_violin(fill = "steelblue", alpha = 0.45, colour = NA) +
  geom_pointrange(data = boot_summ,
                  aes(x = boot_med, xmin = boot_lo95, xmax = boot_hi95,
                      y = reorder(Stream, boot_med)),
                  linewidth = 0.6) +
  labs(x = expression("Mean dW day"^{-1}*"  (g fish"^{-1}*" day"^{-1}*")"),
       y = NULL,
       title    = "Parametric bootstrap: mean daily weight growth per stream",
       subtitle = "B = 5000 draws from Wald distribution of fixed effects; RE fixed") +
  theme_bw(base_size = 12)

print(p_boot)
ggsave("plot_growth_bootstrap.png", p_boot, width = 7, height = 5, dpi = 300)

cat("\nDone.\n")
