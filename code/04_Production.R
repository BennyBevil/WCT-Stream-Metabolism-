# 04_Production.R
# Westslope cutthroat trout — secondary production synthesis
# Aug–Oct 2025, 8 isolated streams
#
# Loads results from:
#   abundance_results.RDS     (01_Abundance.R)
#   bioenergetics_results.RDS (03_Bioenergetics.R — fb4package per-fish simulations)
#
# Workflow:
#   1. Load upstream results
#   2. Monte Carlo production, respiration, and consumption (N_posterior × per-fish fb4 means)
#      Units: kJ m-2 day-1
#   3. Summaries and diagnostic plots
#   4. Fish P, R, C vs GPP (BaMM posteriors)
#   5. Regression: log(P) ~ log(GPP) + Temp + log(GPP):Temp  (lm, weighted)
#   6. Regression: log(C) ~ log(GPP)  (lm, weighted)
#   7. Marginal and interaction effect plots (P model)
#   8. Save combined outputs

setwd("C:/Users/benny/OneDrive/Desktop/Masters.working")

library(tidyverse)
library(MuMIn)
library(patchwork)
library(ggrepel)

# ============================================================
# Toggles
# ============================================================

# Abundance model: simple (single q per stream) vs size-bin (q ~ FL)
USE_BINS <- TRUE

# Depletion occasion used for abundance in the production MC:
#   "dep1" — August standing crop (fish present at start of interval)
#   "dep2" — October standing crop (fish present at end of interval)
#   "mean" — arithmetic mean of Dep1 and Dep2 posteriors (average standing crop)
N_OCCASION <- "dep1"


# ============================================================
# 1. Load upstream results
# ============================================================

abund <- readRDS("abundance_results.RDS")
list2env(abund, envir = environment())
# Provides: N_post_d1, N_post_d2, N_post_bin_d1, N_post_bin_d2,
#           N_total_bin_d1, N_total_bin_d2,
#           bins, Nbin, bin_mids,
#           streams, Nstream, stream_area, interval_days, dep1, dep2

bio <- readRDS("bioenergetics_results.RDS")
fish_bioen  <- bio$fish_bioen     # per-fish: Stream, FL_aug, W_aug, W_final,
                                  #           cum_growth_kJ, cum_resp_kJ, cum_cons_kJ, converged
stream_temp <- bio$stream_temp    # data frame: Stream, mean_temp
temp_vec    <- bio$temp_vec       # named vector: mean temp per stream (aligned to streams)
oxy_to_kJ   <- bio$OXYCAL_J / 1000   # J g-1 O2 -> kJ g-1 O2

grw <- readRDS("growth_results.RDS")
fe_W        <- grw$fe_W          # fixed effects: (Intercept), logW_aug
sigma_W     <- grw$sigma_W       # residual SD on log scale (log g day-1)
stream_re_W <- grw$stream_re_W   # named vector: stream random effects (log g day-1)

# Caloric density for converting growth in wet mass to energy
CALORIC_DENSITY_KJ_PER_G <- 5.7  # kJ g-1 wet mass; salmonid tissue; KEY ASSUMPTION

# Total-N posteriors [n_iter × Nstream] for the abundance comparison plot
N_post_active_d1 <- if (USE_BINS) N_total_bin_d1 else N_post_d1
N_post_active_d2 <- if (USE_BINS) N_total_bin_d2 else N_post_d2

n_iter <- nrow(N_post_active_d1)

# N posteriors used in the production MC — simple [n_iter × Nstream]
N_mc <- switch(N_OCCASION,
  dep1 = N_post_active_d1,
  dep2 = N_post_active_d2,
  mean = (N_post_active_d1 + N_post_active_d2) / 2,
  stop("N_OCCASION must be 'dep1', 'dep2', or 'mean'")
)

# Binned N posteriors used in the production MC — [n_iter × Nstream × Nbin]
N_mc_bin <- switch(N_OCCASION,
  dep1 = N_post_bin_d1,
  dep2 = N_post_bin_d2,
  mean = (N_post_bin_d1 + N_post_bin_d2) / 2,
  stop("N_OCCASION must be 'dep1', 'dep2', or 'mean'")
)

cat(sprintf("\nAbundance model: %s  |  Occasion: %s\n",
            if (USE_BINS) "size-bin" else "simple", N_OCCASION))

# ============================================================
# 2. Monte Carlo production, respiration, and consumption
# ============================================================
# Production (P): dep2 N × resampled mean individual growth from dep1 catch.
#   Per-fish expected daily growth (kJ day-1) from dW model (02_Growth.R):
#     dW_kJ_day = exp(mu_log_dW_day + sigma_W²/2) × CALORIC_DENSITY
#   For each MC iteration the dep1 catch is treated as a finite sample of the
#   stream population. Mean individual growth for that iteration is drawn from
#     Normal(mean(dW_kJ_day_s), sd(dW_kJ_day_s) / sqrt(round(N_mc[iter, s])))
#   which is the CLT approximation to bootstrapping N fish from the dep1 pool.
#   This propagates both N uncertainty (dep2 posterior) and sampling uncertainty
#   in mean individual growth (dep1 catch size).
#   Growth model gives a per-day rate, so no /days_s for production.
# Respiration (R) and consumption (C): fb4 interval totals ÷ days_s.

# Helper: per-fish expected daily growth (kJ day-1) from dW model
dW_kJ_day <- function(w_vec, re) {
  mu <- fe_W["(Intercept)"] + re + fe_W["logW_aug"] * log(w_vec)
  exp(mu + sigma_W^2 / 2) * CALORIC_DENSITY_KJ_PER_G
}

# Helper: iteration-wise mean growth via CLT resampling (vectorized)
# Returns vector of length n_iter; where N_vec == 0 returns 0.
resample_mean <- function(dW_vec, N_vec) {
  m  <- mean(dW_vec)
  s  <- if (length(dW_vec) > 1) sd(dW_vec) else 0
  N_pos <- pmax(round(N_vec), 1L)
  out <- rnorm(length(N_vec), m, s / sqrt(N_pos))
  out[round(N_vec) == 0] <- 0
  out
}

prod_kJ_m2 <- matrix(0, nrow = n_iter, ncol = Nstream,
                     dimnames = list(NULL, streams))
resp_kJ_m2 <- matrix(0, nrow = n_iter, ncol = Nstream,
                     dimnames = list(NULL, streams))
cons_kJ_m2 <- matrix(0, nrow = n_iter, ncol = Nstream,
                     dimnames = list(NULL, streams))

if (!USE_BINS) {

  # --- Simple model ---
  for (s in seq_along(streams)) {
    sname  <- streams[s]
    area_s <- stream_area[s]
    days_s <- interval_days %>% filter(Stream == sname) %>% pull(days)

    fish_s <- fish_bioen %>%
      filter(Stream == sname, converged == TRUE,
             !is.na(cum_resp_kJ), !is.na(cum_cons_kJ))
    if (nrow(fish_s) == 0) next

    re_s <- if (!is.na(stream_re_W[sname])) stream_re_W[sname] else 0
    w_s  <- fish_s$W_aug[!is.na(fish_s$W_aug) & fish_s$W_aug > 0]
    if (length(w_s) == 0) next

    # P: iteration-wise resampled mean growth × dep2 N
    dW_s <- dW_kJ_day(w_s, re_s)
    prod_kJ_m2[, s] <- N_mc[, s] * resample_mean(dW_s, N_mc[, s]) / area_s

    # R and C: fb4 interval totals → daily rate
    resp_kJ_m2[, s] <- N_mc[, s] * mean(fish_s$cum_resp_kJ) / area_s / days_s
    cons_kJ_m2[, s] <- N_mc[, s] * mean(fish_s$cum_cons_kJ) / area_s / days_s
  }

} else {

  # --- Size-bin model ---
  fish_bioen <- fish_bioen %>%
    mutate(fl_bin = as.integer(cut(FL_aug, breaks = bins,
                                   include.lowest = TRUE, right = FALSE)))

  for (s in seq_along(streams)) {
    sname  <- streams[s]
    area_s <- stream_area[s]
    days_s <- interval_days %>% filter(Stream == sname) %>% pull(days)

    fish_s <- fish_bioen %>%
      filter(Stream == sname, converged == TRUE,
             !is.na(cum_resp_kJ), !is.na(cum_cons_kJ))
    if (nrow(fish_s) == 0) next

    re_s <- if (!is.na(stream_re_W[sname])) stream_re_W[sname] else 0
    w_s  <- fish_s$W_aug[!is.na(fish_s$W_aug) & fish_s$W_aug > 0]

    # Stream-level fallback (used when a bin has no dep1 fish)
    dW_s_all       <- if (length(w_s) > 0) dW_kJ_day(w_s, re_s) else numeric(0)
    fallback_resp  <- mean(fish_s$cum_resp_kJ)
    fallback_cons  <- mean(fish_s$cum_cons_kJ)

    for (b in seq_len(Nbin)) {
      fish_sb <- fish_s %>% filter(fl_bin == b)

      N_bin_sb <- N_mc_bin[, s, b]
      N_bin_sb[is.na(N_bin_sb)] <- 0

      if (nrow(fish_sb) == 0 || all(is.na(fish_sb$W_aug) | fish_sb$W_aug <= 0)) {
        if (length(dW_s_all) == 0) next   # no fish anywhere; skip
        warning(sprintf("Stream %s bin %d (%g mm): no dep1 fish — resampling from stream pool.",
                        sname, b, bin_mids[b]))
        dW_b <- dW_s_all
        mr   <- fallback_resp
        mc   <- fallback_cons
      } else {
        w_sb <- fish_sb$W_aug[!is.na(fish_sb$W_aug) & fish_sb$W_aug > 0]
        dW_b <- dW_kJ_day(w_sb, re_s)
        mr   <- mean(fish_sb$cum_resp_kJ)
        mc   <- mean(fish_sb$cum_cons_kJ)
      }

      prod_kJ_m2[, s] <- prod_kJ_m2[, s] + N_bin_sb * resample_mean(dW_b, N_bin_sb) / area_s
      resp_kJ_m2[, s] <- resp_kJ_m2[, s] + N_bin_sb * mr / area_s / days_s
      cons_kJ_m2[, s] <- cons_kJ_m2[, s] + N_bin_sb * mc / area_s / days_s
    }
  }
}

# ============================================================
# 3. Summaries
# ============================================================

post_summary <- function(df, val_col) {
  df %>%
    group_by(Stream) %>%
    summarise(med  = median(.data[[val_col]], na.rm = TRUE),
              lo50 = quantile(.data[[val_col]], 0.25,  na.rm = TRUE),
              hi50 = quantile(.data[[val_col]], 0.75,  na.rm = TRUE),
              lo95 = quantile(.data[[val_col]], 0.025, na.rm = TRUE),
              hi95 = quantile(.data[[val_col]], 0.975, na.rm = TRUE),
              .groups = "drop")
}

prod_long <- as.data.frame(prod_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "P_kJ_m2")
resp_long <- as.data.frame(resp_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "R_kJ_m2")
cons_long <- as.data.frame(cons_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "C_kJ_m2")

prod_ps <- post_summary(prod_long, "P_kJ_m2") %>%
  rename(P_median = med, P_lo50 = lo50, P_hi50 = hi50, P_lo95 = lo95, P_hi95 = hi95)
resp_ps <- post_summary(resp_long, "R_kJ_m2") %>%
  rename(R_median = med, R_lo50 = lo50, R_hi50 = hi50, R_lo95 = lo95, R_hi95 = hi95)
cons_ps <- post_summary(cons_long, "C_kJ_m2") %>%
  rename(C_median = med, C_lo50 = lo50, C_hi50 = hi50, C_lo95 = lo95, C_hi95 = hi95)

pr_ratio_long <- data.frame(
  Stream   = rep(streams, each = n_iter),
  P_kJ_m2  = as.vector(prod_kJ_m2),
  R_kJ_m2  = as.vector(resp_kJ_m2),
  C_kJ_m2  = as.vector(cons_kJ_m2)
) %>% mutate(PR = P_kJ_m2 / R_kJ_m2)

prod_summary <- prod_ps %>%
  left_join(resp_ps, by = "Stream") %>%
  left_join(cons_ps, by = "Stream") %>%
  left_join(interval_days, by = "Stream") %>%
  left_join(
    pr_ratio_long %>%
      group_by(Stream) %>%
      summarise(PR_median = median(PR,           na.rm = TRUE),
                PR_lo95   = quantile(PR, 0.025,  na.rm = TRUE),
                PR_hi95   = quantile(PR, 0.975,  na.rm = TRUE),
                .groups   = "drop"),
    by = "Stream"
  ) %>%
  arrange(desc(P_median))

cat("\n--- Abundance comparison: Dep1 vs Dep2 (median, 95% CRI) ---\n")
fmt_post <- function(mat) apply(mat, 2, function(x)
  sprintf("%d (%d\u2013%d)", round(median(x)),
          round(quantile(x, 0.025)), round(quantile(x, 0.975))))
N_comparison <- data.frame(
  Stream     = streams,
  `Dep1 Aug` = fmt_post(N_post_active_d1),
  `Dep2 Oct` = fmt_post(N_post_active_d2),
  check.names = FALSE
)
print(N_comparison)

cat("\n--- Per-fish fb4 summary by stream (converged fish only; kJ fish-1 interval-1) ---\n")
fish_bioen %>%
  filter(converged == TRUE) %>%
  group_by(Stream) %>%
  summarise(
    n          = n(),
    p_median   = round(median(p_value,       na.rm = TRUE), 3),
    cons_med   = round(median(cum_cons_kJ,   na.rm = TRUE), 2),
    resp_med   = round(median(cum_resp_kJ,   na.rm = TRUE), 2),
    growth_med = round(median(cum_growth_kJ, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

cat("\n--- Production, respiration, and consumption summary (kJ m-2 day-1) ---\n")
prod_summary %>%
  select(Stream, days, P_median, P_lo95, P_hi95,
         R_median, C_median, PR_median) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

# ============================================================
# 4. Diagnostic plots
# ============================================================

# Production posterior
p1 <- ggplot(prod_long,
             aes(x = P_kJ_m2, y = reorder(Stream, P_kJ_m2, median))) +
  geom_violin(fill = "steelblue", alpha = 0.5, colour = NA) +
  geom_pointrange(data = prod_ps,
                  aes(x = P_median, xmin = P_lo95, xmax = P_hi95,
                      y = reorder(Stream, P_median)),
                  linewidth = 0.6) +
  geom_linerange(data = prod_ps,
                 aes(x = P_median, xmin = P_lo50, xmax = P_hi50,
                     y = reorder(Stream, P_median)),
                 linewidth = 1.5) +
  labs(x = expression("Production (kJ m"^{-2}*" day"^{-1}*")"),
       y = NULL,
       title = "WCT secondary production, Aug\u2013Oct 2025",
       subtitle = "Median with 50% (thick) and 95% (thin) posterior intervals") +
  theme_bw(base_size = 12)

# Abundance posteriors — Dep1 vs Dep2 (total N, whichever model is active)
N_long <- bind_rows(
  as.data.frame(N_post_active_d1) %>%
    pivot_longer(everything(), names_to = "Stream", values_to = "N") %>%
    mutate(Occasion = "Depletion 1 (Aug)"),
  as.data.frame(N_post_active_d2) %>%
    pivot_longer(everything(), names_to = "Stream", values_to = "N") %>%
    mutate(Occasion = "Depletion 2 (Oct)")
)

N_ps <- N_long %>%
  group_by(Stream, Occasion) %>%
  summarise(med  = median(N),
            lo50 = quantile(N, 0.25),
            hi50 = quantile(N, 0.75),
            lo95 = quantile(N, 0.025),
            hi95 = quantile(N, 0.975),
            .groups = "drop")

stream_order_N <- N_ps %>%
  filter(Occasion == "Depletion 1 (Aug)") %>%
  arrange(med) %>%
  pull(Stream)

N_long <- N_long %>% mutate(Stream = factor(Stream, levels = stream_order_N))
N_ps   <- N_ps   %>% mutate(Stream = factor(Stream, levels = stream_order_N))

p2 <- ggplot(N_long, aes(x = N, y = Stream, fill = Occasion)) +
  geom_violin(alpha = 0.45, colour = NA, position = position_dodge(width = 0.8)) +
  geom_pointrange(data = N_ps,
                  aes(x = med, xmin = lo95, xmax = hi95, colour = Occasion),
                  position = position_dodge(width = 0.8), linewidth = 0.6,
                  show.legend = FALSE) +
  geom_linerange(data = N_ps,
                 aes(x = med, xmin = lo50, xmax = hi50, colour = Occasion),
                 position = position_dodge(width = 0.8), linewidth = 1.5,
                 show.legend = FALSE) +
  scale_fill_manual(values   = c("Depletion 1 (Aug)" = "tomato",
                                  "Depletion 2 (Oct)" = "steelblue")) +
  scale_colour_manual(values = c("Depletion 1 (Aug)" = "tomato",
                                  "Depletion 2 (Oct)" = "steelblue")) +
  labs(x = "Estimated abundance", y = NULL, fill = NULL,
       title = "Abundance posteriors \u2014 Depletion 1 vs 2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

print(p1)
print(p2)

# ============================================================
# 5. Fish P, R, C vs GPP
# ============================================================

gpp_name_map <- c(
  Buf   = "Buffalo",
  CC    = "CCT",
  CL    = "CLT",
  cow   = "Cow",
  dry   = "Dry",
  Jerry = "Jerry",
  Pint  = "Pintler",
  Plimp = "Plimpton",
  Hen   = "Henry's Fork"
)

# IntegratedPP is already a daily rate (g O2 m-2 day-1) from the BaMM model.
# One posterior file per site in Data/GPP_EST/ (e.g. "Pint_IND_BaMM.RDS"); the
# site code is the filename prefix and matches the keys in gpp_name_map.
gpp_files <- list.files("Data/GPP_EST", pattern = "_IND_BaMM\\.RDS$",
                        full.names = TRUE)

gpp_post <- purrr::map_dfr(gpp_files, function(f) {
  site_code <- sub("_IND_BaMM$", "", tools::file_path_sans_ext(basename(f)))
  readRDS(f) %>%
    transmute(site = site_code, IntegratedPP = Integrated.PP)
}) %>%
  mutate(Stream = gpp_name_map[site]) %>%
  filter(!is.na(Stream), Stream %in% streams) %>%
  mutate(GPP_kJ = IntegratedPP * oxy_to_kJ)

cat("\n--- GPP summary (daily rate, kJ m-2 day-1) ---\n")
gpp_post %>%
  group_by(Stream) %>%
  summarise(GPP_kJ_med  = round(median(GPP_kJ), 2),
            GPP_kJ_lo95 = round(quantile(GPP_kJ, 0.025), 2),
            GPP_kJ_hi95 = round(quantile(GPP_kJ, 0.975), 2),
            .groups = "drop") %>%
  print()

set.seed(42)

ratio_post <- lapply(streams, function(s) {
  gpp_s  <- gpp_post %>% filter(Stream == s) %>% pull(GPP_kJ)
  prod_s <- prod_kJ_m2[, s]
  resp_s <- resp_kJ_m2[, s]
  cons_s <- cons_kJ_m2[, s]
  if (length(gpp_s) == 0) return(NULL)
  gpp_s_matched <- sample(gpp_s, n_iter, replace = TRUE)
  data.frame(
    Stream    = s,
    P_kJ_m2   = prod_s,
    R_kJ_m2   = resp_s,
    C_kJ_m2   = cons_s,
    PR_kJ_m2  = prod_s + resp_s,
    GPP_kJ_m2 = gpp_s_matched,
    P_GPP     = prod_s / gpp_s_matched,
    R_GPP     = resp_s / gpp_s_matched,
    C_GPP     = cons_s / gpp_s_matched,
    PR_GPP    = (prod_s + resp_s) / gpp_s_matched
  )
}) %>% bind_rows()

cat("\n--- P/GPP, C/GPP, and (P+R)/GPP ratio posterior summaries ---\n")
ratio_post %>%
  group_by(Stream) %>%
  summarise(P_GPP_med  = round(median(P_GPP),  4),
            P_GPP_lo95 = round(quantile(P_GPP,  0.025), 4),
            P_GPP_hi95 = round(quantile(P_GPP,  0.975), 4),
            C_GPP_med  = round(median(C_GPP),  4),
            C_GPP_lo95 = round(quantile(C_GPP,  0.025), 4),
            C_GPP_hi95 = round(quantile(C_GPP,  0.975), 4),
            PR_GPP_med = round(median(PR_GPP), 4),
            .groups = "drop") %>%
  arrange(desc(C_GPP_med)) %>%
  print()

# Plots A1–A3: C, P, and P+R vs GPP — separate scatter plots
prod_summ_gpp <- ratio_post %>%
  group_by(Stream) %>%
  summarise(P_med   = median(P_kJ_m2),   P_lo   = quantile(P_kJ_m2,   0.025),
            P_hi    = quantile(P_kJ_m2,   0.975),
            C_med   = median(C_kJ_m2),   C_lo   = quantile(C_kJ_m2,   0.025),
            C_hi    = quantile(C_kJ_m2,   0.975),
            PR_med  = median(PR_kJ_m2),  PR_lo  = quantile(PR_kJ_m2,  0.025),
            PR_hi   = quantile(PR_kJ_m2,  0.975),
            GPP_med = median(GPP_kJ_m2), GPP_lo = quantile(GPP_kJ_m2, 0.025),
            GPP_hi  = quantile(GPP_kJ_m2, 0.975),
            .groups = "drop")

flux_scatter <- function(val_med, val_lo, val_hi, colour, y_lab, title_str, file_str) {
  df <- prod_summ_gpp %>%
    transmute(Stream, GPP_med, GPP_lo, GPP_hi,
              val_med = {{ val_med }}, val_lo = {{ val_lo }}, val_hi = {{ val_hi }})
  p <- ggplot(df, aes(x = GPP_med, y = val_med, label = Stream)) +
    geom_errorbar(aes(ymin = val_lo, ymax = val_hi),
                  width = 0, linewidth = 0.4, colour = colour) +
    geom_errorbar(aes(xmin = GPP_lo, xmax = GPP_hi),
                  width = 0, linewidth = 0.4, colour = colour, orientation = "y") +
    geom_point(colour = colour, size = 3) +
    geom_text_repel(size = 3.2, colour = "grey30") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
         y = y_lab,
         title = title_str,
         subtitle = "Log\u2081\u2080 \u2013 log\u2081\u2080 scale") +
    theme_bw(base_size = 12)
  ggsave(file_str, p, width = 6, height = 5, dpi = 300)
  p
}

pA1 <- flux_scatter(C_med,  C_lo,  C_hi,  "forestgreen",
                    expression("C (kJ m"^{-2}*" day"^{-1}*")"),
                    "Fish consumption vs GPP", "plot_C_vs_GPP.png")
pA2 <- flux_scatter(P_med,  P_lo,  P_hi,  "steelblue",
                    expression("P (kJ m"^{-2}*" day"^{-1}*")"),
                    "Fish production vs GPP", "plot_P_vs_GPP.png")
pA3 <- flux_scatter(PR_med, PR_lo, PR_hi, "tomato",
                    expression("P + R (kJ m"^{-2}*" day"^{-1}*")"),
                    "Fish total demand (P+R) vs GPP", "plot_PR_vs_GPP.png")

print(pA1)
print(pA2)
print(pA3)

# ============================================================
# 6. Regression: log(P) ~ log(GPP) + Temp + log(GPP):Temp  (lm, weighted)
# ============================================================
# WARNING: n = 8 streams, 4 parameters — 4 residual df. Exploratory only.
# No quadratic Temp term: streams span a narrow temperature range.
# Weights = 1 / posterior_variance of log(P) per stream.

gpp_medians      <- sapply(streams, function(s)
  median(gpp_post %>% filter(Stream == s) %>% pull(GPP_kJ)))

log_gpp_medians  <- log(gpp_medians)
log_gpp_mean_ref <- mean(log_gpp_medians)
log_gpp_sd_ref   <- sd(log_gpp_medians)
temp_z           <- (temp_vec - mean(temp_vec)) / sd(temp_vec)
GPP_z_vec        <- (log_gpp_medians - log_gpp_mean_ref) / log_gpp_sd_ref

P_log_med_vec <- apply(prod_kJ_m2, 2, function(x) median(log(x), na.rm = TRUE))
P_log_sd_vec  <- apply(prod_kJ_m2, 2, function(x) sd(log(x),     na.rm = TRUE))
C_log_med_vec <- apply(cons_kJ_m2, 2, function(x) median(log(x), na.rm = TRUE))
C_log_sd_vec  <- apply(cons_kJ_m2, 2, function(x) sd(log(x),     na.rm = TRUE))

reg_df <- data.frame(
  P_log  = P_log_med_vec,
  C_log  = C_log_med_vec,
  GPP_z  = GPP_z_vec,
  Temp_z = temp_z,
  w_P    = 1 / P_log_sd_vec^2,
  w_C    = 1 / C_log_sd_vec^2
)

P_lm <- lm(P_log ~ GPP_z + Temp_z + GPP_z:Temp_z,
           weights = w_P, data = reg_df)

cat("\n--- log(P) ~ log(GPP) + Temp + log(GPP):Temp  [n = 8] ---\n")
cat("NOTE: n=8 streams, 4 parameters — 4 residual df. Exploratory only.\n\n")
print(summary(P_lm))

coef_ci_P <- as.data.frame(cbind(est = coef(P_lm), confint(P_lm)))
colnames(coef_ci_P) <- c("est", "lo95", "hi95")
cat("\n--- P model coefficients (estimate, 95% CI) ---\n")
print(round(coef_ci_P, 3))

# Model selection: all subsets of the global P model, ranked by AICc
# Interaction only appears when both main effects are present (subset argument).
options(na.action = "na.fail")
dredge_P <- dredge(P_lm,
                   subset = !(`GPP_z:Temp_z`) | (GPP_z & Temp_z),
                   rank   = "AICc")
cat("\n--- P model selection (AICc, n = 8) ---\n")
print(dredge_P)
options(na.action = "na.omit")   # restore default

# ============================================================
# 7. Regression: log(C) ~ log(GPP)  (lm, weighted)
# ============================================================

C_lm <- lm(C_log ~ GPP_z, weights = w_C, data = reg_df)

cat("\n--- log(C) ~ log(GPP)  [n = 8] ---\n")
cat("NOTE: n = 8 streams. Exploratory.\n\n")
print(summary(C_lm))

coef_ci_C <- as.data.frame(cbind(est = coef(C_lm), confint(C_lm)))
colnames(coef_ci_C) <- c("est", "lo95", "hi95")
cat("\n--- C model coefficients (estimate, 95% CI) ---\n")
print(round(coef_ci_C, 3))

# ============================================================
# 8. Effect plots
# ============================================================

n_grid <- 150

GPP_z_to_kJ <- function(z) exp(z * log_gpp_sd_ref + log_gpp_mean_ref)
Temp_z_to_C <- function(z) z * sd(temp_vec) + mean(temp_vec)

GPP_z_seq    <- seq(min(GPP_z_vec) - 0.3, max(GPP_z_vec) + 0.3, length.out = n_grid)
GPP_kJ_seq   <- GPP_z_to_kJ(GPP_z_seq)
Temp_z_seq   <- seq(min(temp_z)    - 0.3, max(temp_z)    + 0.3, length.out = n_grid)
Temp_C_seq   <- Temp_z_to_C(Temp_z_seq)

obs_df <- data.frame(
  Stream   = streams,
  P_obs    = apply(prod_kJ_m2, 2, median, na.rm = TRUE),
  GPP_z    = GPP_z_vec,
  Temp_z   = temp_z,
  Temp_raw = temp_vec
) %>% mutate(GPP_kJ = GPP_z_to_kJ(GPP_z))

# Helper: predict on log scale then back-transform, returning data frame
lm_pred <- function(mod, newdata) {
  pr <- predict(mod, newdata = newdata, interval = "confidence")
  data.frame(med = exp(pr[, "fit"]), lo95 = exp(pr[, "lwr"]), hi95 = exp(pr[, "upr"]))
}

# Coefficient forest plot — P model
coef_df_P <- coef_ci_P %>%
  tibble::rownames_to_column("term") %>%
  mutate(term = factor(term, levels = rev(rownames(coef_ci_P))))

pC_reg <- ggplot(coef_df_P, aes(x = est, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo95, xmax = hi95), linewidth = 0.6) +
  labs(x = "Coefficient (log P per SD)", y = NULL,
       title    = "log(P) ~ log(GPP) + Temp + log(GPP):Temp  [n = 8]",
       subtitle = "Estimate \u00B1 95% CI (exploratory)") +
  theme_bw(base_size = 12)

# Fitted vs observed
fitted_df <- obs_df %>%
  mutate(P_fit = exp(fitted(P_lm)))

pD <- ggplot(fitted_df, aes(x = P_fit, y = P_obs, label = Stream)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = Temp_raw), size = 3) +
  geom_text_repel(size = 3.2, colour = "grey30") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_gradient(low = "steelblue", high = "tomato", name = "Mean temp (\u00B0C)") +
  labs(x = expression("Fitted P (kJ m"^{-2}*" day"^{-1}*")"),
       y = expression("Observed P (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "Fitted vs observed production",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080") +
  theme_bw(base_size = 12)

# Plot E: marginal effect of GPP (Temp at mean = 0)
pred_E_df <- cbind(GPP_kJ = GPP_kJ_seq,
  lm_pred(P_lm, data.frame(GPP_z = GPP_z_seq, Temp_z = 0)))

pE <- ggplot(pred_E_df, aes(x = GPP_kJ)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = med), colour = "steelblue", linewidth = 1) +
  geom_point(data = obs_df, aes(x = GPP_kJ, y = P_obs, colour = Temp_raw), size = 3) +
  geom_text_repel(data = obs_df, aes(x = GPP_kJ, y = P_obs, label = Stream),
                  size = 3.2, colour = "grey30") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_gradient(low = "grey70", high = "tomato", name = "Mean temp (\u00B0C)") +
  labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
       y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "Marginal effect of GPP on production",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; temperature held at its mean") +
  theme_bw(base_size = 12)

# Plot F: marginal effect of Temp (GPP at mean = 0)
pred_F_df <- cbind(Temp_C = Temp_C_seq,
  lm_pred(P_lm, data.frame(GPP_z = 0, Temp_z = Temp_z_seq)))

pF <- ggplot(pred_F_df, aes(x = Temp_C)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "tomato", alpha = 0.25) +
  geom_line(aes(y = med), colour = "tomato", linewidth = 1) +
  geom_point(data = obs_df, aes(x = Temp_raw, y = P_obs, colour = GPP_kJ), size = 3) +
  geom_text_repel(data = obs_df, aes(x = Temp_raw, y = P_obs, label = Stream),
                  size = 3.2, colour = "grey30") +
  scale_y_log10() +
  scale_colour_gradient(low = "grey70", high = "steelblue",
                        name = expression("GPP (kJ m"^{-2}*")")) +
  labs(x = "Mean temperature (\u00B0C)",
       y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "Marginal effect of temperature on production",
       subtitle = "Log\u2081\u2080 y-axis; GPP held at its mean") +
  theme_bw(base_size = 12)

# Plot G: GPP effect at low / mean / high temperature
temp_levels <- c(-1, 0, 1)
temp_labels <- paste0(round(Temp_z_to_C(temp_levels), 1), " \u00B0C")

pred_G_df <- bind_rows(lapply(seq_along(temp_levels), function(j) {
  cbind(GPP_kJ = GPP_kJ_seq, Temp_level = temp_labels[j],
        lm_pred(P_lm, data.frame(GPP_z = GPP_z_seq, Temp_z = temp_levels[j])))
})) %>% mutate(Temp_level = factor(Temp_level, levels = temp_labels),
               across(c(med, lo95, hi95), as.numeric))

pG <- ggplot(pred_G_df, aes(x = as.numeric(GPP_kJ), colour = Temp_level,
                             fill = Temp_level)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med), linewidth = 1) +
  geom_point(data = obs_df, aes(x = GPP_kJ, y = P_obs),
             colour = "black", size = 2.5, inherit.aes = FALSE) +
  geom_text_repel(data = obs_df, aes(x = GPP_kJ, y = P_obs, label = Stream),
                  size = 3, colour = "grey30", inherit.aes = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("steelblue", "forestgreen", "tomato"),
                      name = "Temperature") +
  scale_fill_manual(values   = c("steelblue", "forestgreen", "tomato"),
                    name = "Temperature") +
  labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
       y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "GPP effect at low / mean / high temperature",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; shading = 95% CI") +
  theme_bw(base_size = 12)

# Plot H: Temp effect at low / mean / high GPP
gpp_levels <- c(-1, 0, 1)
gpp_labels <- paste0(round(GPP_z_to_kJ(gpp_levels)), " kJ m\u207B\u00B2")

pred_H_df <- bind_rows(lapply(seq_along(gpp_levels), function(j) {
  cbind(Temp_C = Temp_C_seq, GPP_level = gpp_labels[j],
        lm_pred(P_lm, data.frame(GPP_z = gpp_levels[j], Temp_z = Temp_z_seq)))
})) %>% mutate(GPP_level = factor(GPP_level, levels = gpp_labels),
               across(c(med, lo95, hi95), as.numeric))

pH <- ggplot(pred_H_df, aes(x = as.numeric(Temp_C), colour = GPP_level,
                             fill = GPP_level)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med), linewidth = 1) +
  geom_point(data = obs_df, aes(x = Temp_raw, y = P_obs),
             colour = "black", size = 2.5, inherit.aes = FALSE) +
  geom_text_repel(data = obs_df, aes(x = Temp_raw, y = P_obs, label = Stream),
                  size = 3, colour = "grey30", inherit.aes = FALSE) +
  scale_y_log10() +
  scale_colour_manual(values = c("steelblue", "forestgreen", "tomato"),
                      name = "GPP level") +
  scale_fill_manual(values   = c("steelblue", "forestgreen", "tomato"),
                    name = "GPP level") +
  labs(x = "Mean temperature (\u00B0C)",
       y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "Temperature effect at low / mean / high GPP",
       subtitle = "Log\u2081\u2080 y-axis; shading = 95% CI") +
  theme_bw(base_size = 12)

# C ~ GPP scatter with regression line
obs_C_df <- data.frame(
  Stream = streams,
  C_obs  = apply(cons_kJ_m2, 2, median, na.rm = TRUE),
  GPP_kJ = gpp_medians,
  Temp_C = temp_vec
)

pred_C_df <- cbind(GPP_kJ = GPP_kJ_seq,
  lm_pred(C_lm, data.frame(GPP_z = GPP_z_seq)))

pC_scatter <- ggplot(pred_C_df, aes(x = GPP_kJ)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "forestgreen", alpha = 0.25) +
  geom_line(aes(y = med), colour = "forestgreen", linewidth = 1) +
  geom_point(data = obs_C_df, aes(x = GPP_kJ, y = C_obs, colour = Temp_C), size = 3) +
  geom_text_repel(data = obs_C_df, aes(x = GPP_kJ, y = C_obs, label = Stream),
                  size = 3.2, colour = "grey30") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_gradient(low = "grey70", high = "tomato", name = "Mean temp (\u00B0C)") +
  labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
       y = expression("C (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "Fish consumption vs GPP",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; shading = 95% CI") +
  theme_bw(base_size = 12)

coef_df_C <- coef_ci_C %>%
  tibble::rownames_to_column("term") %>%
  mutate(term = factor(term, levels = rev(rownames(coef_ci_C))))

pC_coef <- ggplot(coef_df_C, aes(x = est, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo95, xmax = hi95), linewidth = 0.6) +
  labs(x = "Coefficient (log C per SD log GPP)", y = NULL,
       title    = "log(C) ~ log(GPP)  [n = 8]",
       subtitle = "Estimate \u00B1 95% CI") +
  theme_bw(base_size = 12)

print(pC_reg)
print(pD)
print(pE)
print(pF)
print(pG)
print(pH)
print(pC_scatter)
print(pC_coef)

ggsave("plot_regression_coefs.png", pC_reg,     width = 6, height = 4, dpi = 300)
ggsave("plot_fitted_vs_obs.png",    pD,          width = 6, height = 5, dpi = 300)
ggsave("plot_effect_GPP.png",        pE,         width = 6, height = 5, dpi = 300)
ggsave("plot_effect_Temp.png",       pF,         width = 6, height = 5, dpi = 300)
ggsave("plot_effect_GPP_x_Temp.png", pG,         width = 7, height = 5, dpi = 300)
ggsave("plot_effect_Temp_x_GPP.png", pH,         width = 7, height = 5, dpi = 300)
ggsave("plot_C_vs_GPP.png",          pC_scatter, width = 6, height = 5, dpi = 300)
ggsave("plot_C_GPP_coefs.png",       pC_coef,    width = 5, height = 3, dpi = 300)

(pE + pF) / (pG + pH) + plot_annotation(
  title = "Regression effects: Fish production ~ GPP + Temp + GPP:Temp",
  theme = theme(plot.title = element_text(size = 13, face = "bold"))
)
ggsave("plot_effects_combined.png", width = 13, height = 10, dpi = 300)

# ============================================================
# 9. Save combined outputs
# ============================================================

saveRDS(list(
  N_post_d1     = N_post_d1,
  N_post_d2     = N_post_d2,
  prod_kJ_m2    = prod_kJ_m2,
  resp_kJ_m2    = resp_kJ_m2,
  cons_kJ_m2    = cons_kJ_m2,
  prod_summary  = prod_summary,
  ratio_post    = ratio_post,
  P_lm          = P_lm,         # log(P) ~ log(GPP) + Temp + Temp^2 + log(GPP):Temp
  C_lm          = C_lm,         # log(C) ~ log(GPP)
  streams       = streams,
  stream_area   = stream_area,
  interval_days = interval_days,
  CALORIC_DENSITY_KJ_PER_G = CALORIC_DENSITY_KJ_PER_G
), "production_results.RDS")

write_csv(prod_summary, "production_summary.csv")

cat("\nDone. Results saved to production_results.RDS and production_summary.csv\n")
