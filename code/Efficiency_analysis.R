# Efficiency_analysis.R
# Trophic transfer efficiency: does temperature predict how much GPP reaches fish?
#
# Motivation: the proposed mechanism is that high food availability (GPP) compensates
# for warm-temperature metabolic costs. The mechanistically direct operationalization
# is trophic transfer efficiency = fish energy demand / GPP supply.
#
# Two response variables:
#   (P+R) / GPP  -- total fish metabolic demand as a fraction of GPP
#   P     / GPP  -- growth production only as a fraction of GPP
#
# Regression: ln(efficiency) ~ Temp_z + Temp_z^2   [n=7, 3 parameters]
# No concavity constraint -- let data determine shape.
#
# Diagnostic: residuals from temperature model vs GPP tests food-compensation
# directly: do high-GPP streams sit above the temperature-only curve?

library(tidyverse)
library(lubridate)
library(R2jags)
library(ggrepel)
library(patchwork)

# ============================================================
# 1. Load saved production results
# ============================================================

res           <- readRDS("production_results.RDS")
prod_kJ_m2    <- res$prod_kJ_m2    # kJ m-2 day-1, matrix [n_iter x Nstream]
resp_kJ_m2    <- res$resp_kJ_m2    # kJ m-2 day-1, matrix [n_iter x Nstream]
streams       <- res$streams
interval_days <- res$interval_days  # Stream, date_d1, date_d2, days
n_iter        <- nrow(prod_kJ_m2)
Nstream       <- length(streams)

# ============================================================
# 2. Daily GPP per stream (kJ m-2 day-1)
# ============================================================

oxy_to_kJ <- 14.1

gpp_name_map <- c(Buf = "Buffalo", CC = "CCT", CL = "CLT", cow = "Cow",
                  dry = "Dry", Jerry = "Jerry", Pint = "Pintler", Plimp = "Plimpton")

temp_dat <- read.csv("data/temp.csv") %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y"))

deploy_days <- temp_dat %>%
  group_by(Stream) %>%
  summarise(deploy_days = as.numeric(max(Date) - min(Date)) + 1, .groups = "drop") %>%
  mutate(Stream = recode(Stream, CC = "CCT"))

gpp_post <- readRDS("data/GPP_EST/IntegratedPP_posteriors.RDS") %>%
  mutate(Stream = gpp_name_map[site]) %>%
  filter(!is.na(Stream), Stream %in% streams) %>%
  left_join(deploy_days, by = "Stream") %>%
  mutate(GPP_kJ = IntegratedPP / deploy_days * oxy_to_kJ)   # kJ m-2 day-1

# ============================================================
# 3. Mean temperature over production interval per stream
# ============================================================

temp_full <- read.csv("data/temp.csv") %>%
  mutate(Date   = as.Date(Date, "%m/%d/%Y"),
         Stream = recode(Stream, CC = "CCT"))

stream_temp <- temp_full %>%
  left_join(interval_days %>% select(Stream, date_d1, date_d2), by = "Stream") %>%
  filter(Date >= date_d1 & Date <= date_d2) %>%
  group_by(Stream) %>%
  summarise(mean_temp = mean(Temp, na.rm = TRUE), .groups = "drop")

temp_vec <- stream_temp %>%
  arrange(match(Stream, streams)) %>%
  pull(mean_temp)

cat("\n--- Mean temperature over production interval ---\n")
print(stream_temp)

# ============================================================
# 4. Efficiency posteriors per stream
# ============================================================

set.seed(42)

eff_post <- lapply(streams, function(s) {
  gpp_s  <- gpp_post %>% filter(Stream == s) %>% pull(GPP_kJ)
  prod_s <- prod_kJ_m2[, s]
  resp_s <- resp_kJ_m2[, s]

  gpp_matched <- sample(gpp_s, n_iter, replace = TRUE)

  data.frame(
    Stream    = s,
    P_kJ      = prod_s,
    R_kJ      = resp_s,
    GPP_kJ    = gpp_matched,
    P_GPP     = prod_s / gpp_matched,
    PR_GPP    = (prod_s + resp_s) / gpp_matched,
    lnPR_GPP  = log((prod_s + resp_s) / gpp_matched),
    lnP_GPP   = log(prod_s / gpp_matched)
  )
}) %>%
  bind_rows()

eff_summary <- eff_post %>%
  group_by(Stream) %>%
  summarise(
    GPP_med      = median(GPP_kJ),
    P_med        = median(P_kJ),
    PR_GPP_med   = median(PR_GPP),
    PR_GPP_lo95  = quantile(PR_GPP, 0.025),
    PR_GPP_hi95  = quantile(PR_GPP, 0.975),
    P_GPP_med    = median(P_GPP),
    P_GPP_lo95   = quantile(P_GPP, 0.025),
    P_GPP_hi95   = quantile(P_GPP, 0.975),
    lnPR_GPP_med = median(lnPR_GPP),
    lnPR_GPP_sd  = sd(lnPR_GPP),
    lnP_GPP_med  = median(lnP_GPP),
    lnP_GPP_sd   = sd(lnP_GPP),
    .groups = "drop"
  ) %>%
  left_join(stream_temp, by = "Stream")

cat("\n--- Trophic transfer efficiency summary ---\n")
eff_summary %>%
  select(Stream, mean_temp, GPP_med,
         P_GPP_med, P_GPP_lo95, P_GPP_hi95,
         PR_GPP_med, PR_GPP_lo95, PR_GPP_hi95) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  arrange(desc(mean_temp)) %>%
  print()

# ============================================================
# 5. Bayesian regression: ln(efficiency) ~ Temp_z + Temp_z^2
# ============================================================
# GPP is absorbed into the response; temperature is the sole predictor.
# No concavity constraint -- let data determine shape.
# Uncertainty in efficiency enters as observation error (posterior SD per stream).

temp_z <- (temp_vec - mean(temp_vec)) / sd(temp_vec)
names(temp_z) <- streams

# Regression model (shared for both response variables)
sink("Efficiency_regression.txt")
cat("
model {
  for(s in 1:Nstream) {
    y_obs[s] ~ dnorm(mu[s], pow(sigma_y[s], -2))
    mu[s] <- beta0 + beta_T * Temp_z[s] + beta_T2 * pow(Temp_z[s], 2)
  }
  beta0   ~ dnorm(0, 0.01)
  beta_T  ~ dnorm(0, 0.01)
  beta_T2 ~ dnorm(0, 0.01)
}
", fill = TRUE)
sink()

# --- (P+R)/GPP regression ---
lnPR_med <- eff_summary %>% arrange(match(Stream, streams)) %>% pull(lnPR_GPP_med)
lnPR_sd  <- eff_summary %>% arrange(match(Stream, streams)) %>% pull(lnPR_GPP_sd)

cat("\nRunning (P+R)/GPP efficiency regression...\n")
fit_PR <- jags(
  data               = list(y_obs   = lnPR_med,
                             sigma_y = lnPR_sd,
                             Temp_z  = temp_z,
                             Nstream = Nstream),
  inits              = NULL,
  parameters.to.save = c("beta0", "beta_T", "beta_T2"),
  model.file         = "Efficiency_regression.txt",
  n.chains = 3, n.thin = 5, n.iter = 50000, n.burnin = 25000,
  DIC = FALSE
)

coef_PR  <- fit_PR$BUGSoutput$sims.matrix[, c("beta0", "beta_T", "beta_T2")]
n_reg    <- nrow(coef_PR)

cat("\n--- Regression: ln((P+R)/GPP) ~ Temp_z + Temp_z^2 ---\n")
cat("n =", Nstream, "streams, 3 parameters (", Nstream - 3, "residual df )\n\n")
print(t(apply(coef_PR, 2, function(x)
  c(median = round(median(x), 3),
    lo95   = round(quantile(x, 0.025), 3),
    hi95   = round(quantile(x, 0.975), 3)))))

# --- P/GPP regression ---
lnP_med <- eff_summary %>% arrange(match(Stream, streams)) %>% pull(lnP_GPP_med)
lnP_sd  <- eff_summary %>% arrange(match(Stream, streams)) %>% pull(lnP_GPP_sd)

cat("\nRunning P/GPP efficiency regression...\n")
fit_P <- jags(
  data               = list(y_obs   = lnP_med,
                             sigma_y = lnP_sd,
                             Temp_z  = temp_z,
                             Nstream = Nstream),
  inits              = NULL,
  parameters.to.save = c("beta0", "beta_T", "beta_T2"),
  model.file         = "Efficiency_regression.txt",
  n.chains = 3, n.thin = 5, n.iter = 50000, n.burnin = 25000,
  DIC = FALSE
)

coef_P <- fit_P$BUGSoutput$sims.matrix[, c("beta0", "beta_T", "beta_T2")]

cat("\n--- Regression: ln(P/GPP) ~ Temp_z + Temp_z^2 ---\n")
cat("n =", Nstream, "streams, 3 parameters (", Nstream - 3, "residual df )\n\n")
print(t(apply(coef_P, 2, function(x)
  c(median = round(median(x), 3),
    lo95   = round(quantile(x, 0.025), 3),
    hi95   = round(quantile(x, 0.975), 3)))))

# ============================================================
# 6. Plots
# ============================================================

Temp_z_to_C <- function(z) z * sd(temp_vec) + mean(temp_vec)

n_grid     <- 150
Tz_seq     <- seq(min(temp_z) - 0.4, max(temp_z) + 0.4, length.out = n_grid)
Tc_seq     <- Temp_z_to_C(Tz_seq)

pred_summary <- function(pred_mat)
  data.frame(
    med  = apply(pred_mat, 1, median),
    lo50 = apply(pred_mat, 1, quantile, 0.25),
    hi50 = apply(pred_mat, 1, quantile, 0.75),
    lo95 = apply(pred_mat, 1, quantile, 0.025),
    hi95 = apply(pred_mat, 1, quantile, 0.975)
  )

# Predicted curves on log scale, then back-transform
make_pred_df <- function(coefs) {
  pred_mat <- outer(Tz_seq, seq_len(nrow(coefs)), function(t, i)
    coefs[i, "beta0"] + coefs[i, "beta_T"] * t + coefs[i, "beta_T2"] * t^2)
  cbind(Temp_C = Tc_seq, pred_summary(pred_mat)) %>%
    mutate(across(c(med, lo50, hi50, lo95, hi95), exp))
}

pred_PR_df <- make_pred_df(coef_PR)
pred_P_df  <- make_pred_df(coef_P)

# --- Plot 1: (P+R)/GPP vs temperature, colored by GPP ---
p1 <- ggplot(eff_summary,
             aes(x = mean_temp, y = PR_GPP_med,
                 colour = GPP_med, label = Stream)) +
  geom_ribbon(data = pred_PR_df,
              aes(x = Temp_C, ymin = lo95, ymax = hi95),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
  geom_ribbon(data = pred_PR_df,
              aes(x = Temp_C, ymin = lo50, ymax = hi50),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.4) +
  geom_line(data = pred_PR_df, aes(x = Temp_C, y = med),
            inherit.aes = FALSE, colour = "grey20", linewidth = 1) +
  geom_errorbar(aes(ymin = PR_GPP_lo95, ymax = PR_GPP_hi95),
                width = 0, linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.2, colour = "grey20") +
  scale_colour_gradient(low = "steelblue", high = "tomato",
                        name = expression("GPP (kJ m"^{-2}*" day"^{-1}*")")) +
  labs(
    x        = "Mean temperature (C)",
    y        = "(P + R) / GPP",
    title    = "Total fish energy demand as a fraction of GPP",
    subtitle = "Colour = GPP; streams above curve exceed temperature-only prediction"
  ) +
  theme_bw(base_size = 12)

# --- Plot 2: P/GPP vs temperature, colored by GPP ---
p2 <- ggplot(eff_summary,
             aes(x = mean_temp, y = P_GPP_med,
                 colour = GPP_med, label = Stream)) +
  geom_ribbon(data = pred_P_df,
              aes(x = Temp_C, ymin = lo95, ymax = hi95),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
  geom_ribbon(data = pred_P_df,
              aes(x = Temp_C, ymin = lo50, ymax = hi50),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.4) +
  geom_line(data = pred_P_df, aes(x = Temp_C, y = med),
            inherit.aes = FALSE, colour = "grey20", linewidth = 1) +
  geom_errorbar(aes(ymin = P_GPP_lo95, ymax = P_GPP_hi95),
                width = 0, linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.2, colour = "grey20") +
  scale_colour_gradient(low = "steelblue", high = "tomato",
                        name = expression("GPP (kJ m"^{-2}*" day"^{-1}*")")) +
  labs(
    x        = "Mean temperature (C)",
    y        = "P / GPP",
    title    = "Growth production as a fraction of GPP",
    subtitle = "Colour = GPP"
  ) +
  theme_bw(base_size = 12)

# --- Plot 3: residuals from temperature model vs GPP (food-compensation diagnostic) ---
# Residual = observed ln((P+R)/GPP) - predicted by temperature-only model.
# Positive residual: stream transfers GPP to fish more efficiently than temperature
# alone predicts. If the food-compensation mechanism is real, high-GPP streams
# should trend positive. NOTE: this is a visual diagnostic, not an independent test --
# the same data informed the regression.

b_PR <- apply(coef_PR, 2, median)

eff_summary <- eff_summary %>%
  mutate(
    Temp_z_s    = (mean_temp - mean(temp_vec)) / sd(temp_vec),
    lnPR_fitted = b_PR["beta0"] +
                  b_PR["beta_T"]  * Temp_z_s +
                  b_PR["beta_T2"] * Temp_z_s^2,
    lnPR_resid  = lnPR_GPP_med - lnPR_fitted
  )

p3 <- ggplot(eff_summary,
             aes(x = GPP_med, y = lnPR_resid,
                 colour = mean_temp, label = Stream)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(size = 3.2, colour = "grey20") +
  scale_colour_gradient(low = "steelblue", high = "tomato",
                        name = "Mean temp (C)") +
  labs(
    x        = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
    y        = "Residual ln((P+R)/GPP)",
    title    = "Food-compensation diagnostic",
    subtitle = "Above zero: more efficient than temperature alone predicts (visual diagnostic only)"
  ) +
  theme_bw(base_size = 12)

print(p1)
print(p2)
print(p3)

(p1 + p2) / p3 +
  plot_annotation(
    title = "Trophic transfer efficiency",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )
ggsave("plot_efficiency_combined.png",   width = 12, height = 10, dpi = 300)
ggsave("plot_efficiency_PR_GPP.png", p1, width = 7,  height = 5,  dpi = 300)
ggsave("plot_efficiency_P_GPP.png",  p2, width = 7,  height = 5,  dpi = 300)
ggsave("plot_efficiency_residuals.png", p3, width = 6, height = 5, dpi = 300)

# ============================================================
# 7. Save
# ============================================================

saveRDS(list(
  eff_post    = eff_post,
  eff_summary = eff_summary,
  coef_PR     = coef_PR,
  coef_P      = coef_P,
  temp_z      = temp_z,
  temp_vec    = temp_vec
), "efficiency_results.RDS")

cat("\nDone. Results saved to efficiency_results.RDS\n")
