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
#   5. Regression: P ~ log(GPP) + Temp + log(GPP):Temp  (NATURAL-scale P;
#      identity link so non-positive production, e.g. Jerry, is admissible)
#   6. Regression: log(C) ~ log(GPP)  (C always positive, kept on log scale)
#   7. Marginal and interaction effect plots (P model)
#   8. Save combined outputs

library(tidyverse)
library(MuMIn)
library(patchwork)
library(ggrepel)
library(R2jags)

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
fe_FL        <- grw$fe_FL         # fixed effects: (Intercept), logFL_aug  (log mm day-1)
sigma_FL     <- grw$sigma_FL      # residual SD on log scale (log mm day-1)
stream_re_FL <- grw$stream_re_FL  # named vector: stream random effects (log mm day-1)
pred_W       <- grw$pred_W        # hierarchical L-W: pred_W(FL, stream, occasion) -> g
lw_b         <- grw$lw_b          # L-W exponent (retained for reference)

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
# Production (P): N × resampled mean individual growth from dep1 catch.
#   Per-fish expected daily growth derived from the dFL model (02_Growth.R):
#     E[dFL_day] = exp(mu_log_dFL_day + sigma_FL²/2)       [mm day-1, lognormal mean]
#   Length is projected over the interval and mapped to weight via the
#   hierarchical L-W model (02_Growth.R):
#     FL_oct = FL_aug + dFL_day × days
#   Two non-negative production variants are computed (both reported):
#     "gross"  — occasion-specific L-W (W_oct via Oct "dep2" curve, which carries
#                the global ~-10% Aug->Oct condition penalty), but per-fish
#                increments are clamped at zero: dW = max(W_oct - W_aug, 0).
#                This is the Allen/Benke gross-production convention.
#     "struct" — structural (length-based) growth only: FL_oct is converted with
#                the Aug "dep1" L-W curve (no condition penalty), so
#                dW = W_oct - W_aug >= 0 always (since FL_oct >= FL_aug).
#   In both cases W_aug = pred_W(FL_aug, stream, "dep1").
#   Converted to energy: dW/days × CALORIC_DENSITY_KJ_PER_G  [kJ day-1]
#   dFL is used (not direct dW) because measured dW is noisy: 33% of fish
#   showed apparent negative dW vs only 9% negative dFL.
# Respiration (R) and consumption (C): fb4 interval totals ÷ days_s.

# Helper: per-fish expected daily growth (kJ day-1) via dFL model + hierarchical L-W.
# fl_vec: FL_aug (mm), stream: stream name, re: stream RE (log mm day-1), days: interval length.
# variant: "gross" (occasion-specific L-W, increments clamped >= 0) or
#          "struct" (Aug L-W curve only, always >= 0).
dW_kJ_day <- function(fl_vec, stream, re, days, variant = c("gross", "struct")) {
  variant <- match.arg(variant)
  dFL_day <- exp(fe_FL["(Intercept)"] + re +
                 fe_FL["logFL_aug"] * log(fl_vec) + sigma_FL^2 / 2)  # mm day-1
  FL_oct  <- fl_vec + dFL_day * days
  W_aug   <- pred_W(fl_vec, stream, "dep1")   # g, Aug L-W relationship
  if (variant == "gross") {
    W_oct <- pred_W(FL_oct, stream, "dep2")   # g, Oct L-W (with condition penalty)
    dW    <- pmax(W_oct - W_aug, 0)           # clamp negative per-fish increments
  } else {
    W_oct <- pred_W(FL_oct, stream, "dep1")   # g, Aug L-W (no condition penalty)
    dW    <- W_oct - W_aug                    # always >= 0 (FL_oct >= FL_aug)
  }
  (dW / days) * CALORIC_DENSITY_KJ_PER_G                              # kJ day-1
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

prod_gross_kJ_m2  <- matrix(0, nrow = n_iter, ncol = Nstream,
                            dimnames = list(NULL, streams))
prod_struct_kJ_m2 <- matrix(0, nrow = n_iter, ncol = Nstream,
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

    re_s   <- if (!is.na(stream_re_FL[sname])) stream_re_FL[sname] else 0
    valid_s <- !is.na(fish_s$FL_aug) & fish_s$FL_aug > 0   # W predicted from FL via L-W
    fl_s   <- fish_s$FL_aug[valid_s]
    if (length(fl_s) == 0) next

    # P: iteration-wise resampled mean growth × N (gross and structural variants)
    dW_s_gross  <- dW_kJ_day(fl_s, sname, re_s, days_s, "gross")
    dW_s_struct <- dW_kJ_day(fl_s, sname, re_s, days_s, "struct")
    prod_gross_kJ_m2[, s]  <- N_mc[, s] * resample_mean(dW_s_gross,  N_mc[, s]) / area_s
    prod_struct_kJ_m2[, s] <- N_mc[, s] * resample_mean(dW_s_struct, N_mc[, s]) / area_s

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

    re_s    <- if (!is.na(stream_re_FL[sname])) stream_re_FL[sname] else 0
    valid_s  <- !is.na(fish_s$FL_aug) & fish_s$FL_aug > 0   # W predicted from FL via L-W
    fl_s    <- fish_s$FL_aug[valid_s]

    # Stream-level fallback (used when a bin has no dep1 fish)
    dW_s_all_gross  <- if (length(fl_s) > 0) dW_kJ_day(fl_s, sname, re_s, days_s, "gross")  else numeric(0)
    dW_s_all_struct <- if (length(fl_s) > 0) dW_kJ_day(fl_s, sname, re_s, days_s, "struct") else numeric(0)
    fallback_resp <- mean(fish_s$cum_resp_kJ)
    fallback_cons <- mean(fish_s$cum_cons_kJ)

    for (b in seq_len(Nbin)) {
      fish_sb <- fish_s %>% filter(fl_bin == b)

      N_bin_sb <- N_mc_bin[, s, b]
      N_bin_sb[is.na(N_bin_sb)] <- 0

      valid_sb <- !is.na(fish_sb$FL_aug) & fish_sb$FL_aug > 0   # W predicted from FL via L-W
      if (nrow(fish_sb) == 0 || !any(valid_sb)) {
        if (length(dW_s_all_gross) == 0) next   # no fish anywhere; skip
        warning(sprintf("Stream %s bin %d (%g mm): no dep1 fish — resampling from stream pool.",
                        sname, b, bin_mids[b]))
        dW_b_gross  <- dW_s_all_gross
        dW_b_struct <- dW_s_all_struct
        mr   <- fallback_resp
        mc   <- fallback_cons
      } else {
        fl_sb <- fish_sb$FL_aug[valid_sb]
        dW_b_gross  <- dW_kJ_day(fl_sb, sname, re_s, days_s, "gross")
        dW_b_struct <- dW_kJ_day(fl_sb, sname, re_s, days_s, "struct")
        mr    <- mean(fish_sb$cum_resp_kJ)
        mc    <- mean(fish_sb$cum_cons_kJ)
      }

      prod_gross_kJ_m2[, s]  <- prod_gross_kJ_m2[, s]  + N_bin_sb * resample_mean(dW_b_gross,  N_bin_sb) / area_s
      prod_struct_kJ_m2[, s] <- prod_struct_kJ_m2[, s] + N_bin_sb * resample_mean(dW_b_struct, N_bin_sb) / area_s
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

# Two non-negative production variants: gross (occasion-specific L-W, clamped) and
# structural (length-based, Aug L-W). Both summarised and reported.
prod_gross_long <- as.data.frame(prod_gross_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "P_kJ_m2")
prod_struct_long <- as.data.frame(prod_struct_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "P_kJ_m2")
resp_long <- as.data.frame(resp_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "R_kJ_m2")
cons_long <- as.data.frame(cons_kJ_m2) %>%
  pivot_longer(everything(), names_to = "Stream", values_to = "C_kJ_m2")

prod_gross_ps <- post_summary(prod_gross_long, "P_kJ_m2") %>%
  rename(P_gross_median = med, P_gross_lo50 = lo50, P_gross_hi50 = hi50,
         P_gross_lo95 = lo95, P_gross_hi95 = hi95)
prod_struct_ps <- post_summary(prod_struct_long, "P_kJ_m2") %>%
  rename(P_struct_median = med, P_struct_lo50 = lo50, P_struct_hi50 = hi50,
         P_struct_lo95 = lo95, P_struct_hi95 = hi95)
resp_ps <- post_summary(resp_long, "R_kJ_m2") %>%
  rename(R_median = med, R_lo50 = lo50, R_hi50 = hi50, R_lo95 = lo95, R_hi95 = hi95)
cons_ps <- post_summary(cons_long, "C_kJ_m2") %>%
  rename(C_median = med, C_lo50 = lo50, C_hi50 = hi50, C_lo95 = lo95, C_hi95 = hi95)

# PR uses gross P (the primary-comparable variant): P+R demand vs respiration.
pr_ratio_long <- data.frame(
  Stream   = rep(streams, each = n_iter),
  P_kJ_m2  = as.vector(prod_gross_kJ_m2),
  R_kJ_m2  = as.vector(resp_kJ_m2),
  C_kJ_m2  = as.vector(cons_kJ_m2)
) %>% mutate(PR = P_kJ_m2 / R_kJ_m2)

prod_summary <- prod_gross_ps %>%
  left_join(prod_struct_ps, by = "Stream") %>%
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
  arrange(desc(P_gross_median))

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

cat("\n--- Production (gross & structural), respiration, consumption summary (kJ m-2 day-1) ---\n")
prod_summary %>%
  select(Stream, days, P_gross_median, P_struct_median,
         R_median, C_median, PR_median) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

# ============================================================
# 4. Diagnostic plots
# ============================================================

# Production posterior — one violin plot per variant (gross, structural)
prod_violin <- function(prod_long_v, prod_ps_v, med_col, lo50_col, hi50_col,
                        lo95_col, hi95_col, title_str) {
  ggplot(prod_long_v,
         aes(x = P_kJ_m2, y = reorder(Stream, P_kJ_m2, median))) +
    geom_violin(fill = "steelblue", alpha = 0.5, colour = NA) +
    geom_pointrange(data = prod_ps_v,
                    aes(x = .data[[med_col]], xmin = .data[[lo95_col]], xmax = .data[[hi95_col]],
                        y = reorder(Stream, .data[[med_col]])),
                    linewidth = 0.6) +
    geom_linerange(data = prod_ps_v,
                   aes(x = .data[[med_col]], xmin = .data[[lo50_col]], xmax = .data[[hi50_col]],
                       y = reorder(Stream, .data[[med_col]])),
                   linewidth = 1.5) +
    labs(x = expression("Production (kJ m"^{-2}*" day"^{-1}*")"),
         y = NULL,
         title = title_str,
         subtitle = "Median with 50% (thick) and 95% (thin) posterior intervals") +
    theme_bw(base_size = 12)
}

p1_gross <- prod_violin(prod_gross_long, prod_gross_ps,
                        "P_gross_median", "P_gross_lo50", "P_gross_hi50",
                        "P_gross_lo95", "P_gross_hi95",
                        "WCT gross secondary production, Aug\u2013Oct 2025")
p1_struct <- prod_violin(prod_struct_long, prod_struct_ps,
                         "P_struct_median", "P_struct_lo50", "P_struct_hi50",
                         "P_struct_lo95", "P_struct_hi95",
                         "WCT structural secondary production, Aug\u2013Oct 2025")

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

print(p1_gross)
print(p1_struct)
print(p2)

# ============================================================
# 5. Fish P, R, C vs GPP
# ============================================================

# GPP from the BaMM respiration model: per-stream posterior draws (column
# Integrated.PP) produced by BAMM/02_run_BaMM.R. Prefer the 2-stage fit
# (Data/GPP_EST/{site}_2stage_BaMM.RDS); fall back to the 1-stage fit
# (Data/GPP_EST/{site}_IND_BaMM.RDS) for streams where the 2-stage model failed
# to converge (e.g. cow — singular Hessian). Site codes are keyed
# case-insensitively because the RDS filenames vary in case ("Dry", "plimp").
gpp_name_map <- c(
  BUF   = "Buffalo",
  CC    = "CCT",
  CL    = "CLT",
  COW   = "Cow",
  DRY   = "Dry",
  HENRY = "Henry",
  JERRY = "Jerry",
  PINT  = "Pintler",
  PLIMP = "Plimpton"
)

# Integrated.PP is a daily rate in mg O2 m-2 d-1 (BaMM.tpl report line; confirmed
# from source). Divide by 1000 to convert mg -> g before applying oxy_to_kJ
# (kJ g-1 O2). Streams not in the production set (Henry, Pintler) drop out at the
# Stream %in% streams filter.
gpp_dir   <- "Data/GPP_EST"
gpp_sites <- sub("_IND_BaMM\\.RDS$", "",
                 basename(list.files(gpp_dir, pattern = "_IND_BaMM\\.RDS$")))

gpp_post <- lapply(gpp_sites, function(site) {
  f2 <- file.path(gpp_dir, paste0(site, "_2stage_BaMM.RDS"))
  f1 <- file.path(gpp_dir, paste0(site, "_IND_BaMM.RDS"))
  use_2stage <- file.exists(f2)
  draws <- readRDS(if (use_2stage) f2 else f1)
  data.frame(site         = site,
             Stream       = unname(gpp_name_map[toupper(site)]),
             gpp_model    = if (use_2stage) "2stage" else "1stage",
             IntegratedPP = draws$Integrated.PP)
}) %>%
  bind_rows() %>%
  filter(!is.na(Stream), Stream %in% streams) %>%
  mutate(GPP_kJ = IntegratedPP / 1000 * oxy_to_kJ)

cat("\n--- GPP source model per stream ---\n")
gpp_post %>% distinct(Stream, gpp_model) %>% arrange(Stream) %>% print()

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
  gpp_s    <- gpp_post %>% filter(Stream == s) %>% pull(GPP_kJ)
  prod_g_s <- prod_gross_kJ_m2[, s]
  prod_s_s <- prod_struct_kJ_m2[, s]
  resp_s   <- resp_kJ_m2[, s]
  cons_s   <- cons_kJ_m2[, s]
  if (length(gpp_s) == 0) return(NULL)
  gpp_s_matched <- sample(gpp_s, n_iter, replace = TRUE)
  data.frame(
    Stream        = s,
    P_gross_kJ_m2  = prod_g_s,
    P_struct_kJ_m2 = prod_s_s,
    R_kJ_m2       = resp_s,
    C_kJ_m2       = cons_s,
    PR_kJ_m2      = prod_g_s + resp_s,   # P+R demand uses gross P
    GPP_kJ_m2     = gpp_s_matched,
    P_gross_GPP   = prod_g_s / gpp_s_matched,
    P_struct_GPP  = prod_s_s / gpp_s_matched,
    R_GPP         = resp_s / gpp_s_matched,
    C_GPP         = cons_s / gpp_s_matched,
    PR_GPP        = (prod_g_s + resp_s) / gpp_s_matched
  )
}) %>% bind_rows()

cat("\n--- P/GPP (gross & structural), C/GPP, and (P+R)/GPP ratio posterior summaries ---\n")
ratio_post %>%
  group_by(Stream) %>%
  summarise(Pg_GPP_med = round(median(P_gross_GPP),  4),
            Ps_GPP_med = round(median(P_struct_GPP), 4),
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
  summarise(Pg_med  = median(P_gross_kJ_m2),  Pg_lo  = quantile(P_gross_kJ_m2,  0.025),
            Pg_hi   = quantile(P_gross_kJ_m2,  0.975),
            Ps_med  = median(P_struct_kJ_m2), Ps_lo  = quantile(P_struct_kJ_m2, 0.025),
            Ps_hi   = quantile(P_struct_kJ_m2, 0.975),
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
pA2g <- flux_scatter(Pg_med, Pg_lo, Pg_hi, "steelblue",
                     expression("P (kJ m"^{-2}*" day"^{-1}*")"),
                     "Fish gross production vs GPP", "plot_P_gross_vs_GPP.png")
pA2s <- flux_scatter(Ps_med, Ps_lo, Ps_hi, "steelblue",
                     expression("P (kJ m"^{-2}*" day"^{-1}*")"),
                     "Fish structural production vs GPP", "plot_P_struct_vs_GPP.png")
pA3 <- flux_scatter(PR_med, PR_lo, PR_hi, "tomato",
                    expression("P + R (kJ m"^{-2}*" day"^{-1}*")"),
                    "Fish total demand (P+R) vs GPP", "plot_PR_vs_GPP.png")

print(pA1)
print(pA2g)
print(pA2s)
print(pA3)

# ============================================================
# 6. Regression setup (shared by both JAGS and lm)
# ============================================================

gpp_medians      <- sapply(streams, function(s)
  median(gpp_post %>% filter(Stream == s) %>% pull(GPP_kJ)))

log_gpp_medians  <- log(gpp_medians)
log_gpp_mean_ref <- mean(log_gpp_medians)
log_gpp_sd_ref   <- sd(log_gpp_medians)
temp_z           <- (temp_vec - mean(temp_vec)) / sd(temp_vec)
GPP_z_vec        <- (log_gpp_medians - log_gpp_mean_ref) / log_gpp_sd_ref

# Both P variants are regressed on the LOG scale inside fit_P_regression()
# (gross and structural production are strictly positive). C stays on log scale.
C_log_med_vec <- apply(cons_kJ_m2, 2, function(x) median(log(x), na.rm = TRUE))

# Per-stream median production for each variant (named by stream)
P_gross_med_vec  <- apply(prod_gross_kJ_m2,  2, median, na.rm = TRUE)
P_struct_med_vec <- apply(prod_struct_kJ_m2, 2, median, na.rm = TRUE)

reg_df <- data.frame(
  C_log  = C_log_med_vec,
  GPP_z  = GPP_z_vec,
  Temp_z = temp_z
)

# Prediction grid and back-transform helpers (shared)
n_grid <- 150

GPP_z_to_kJ <- function(z) exp(z * log_gpp_sd_ref + log_gpp_mean_ref)
Temp_z_to_C <- function(z) z * sd(temp_vec) + mean(temp_vec)

GPP_z_seq  <- seq(min(GPP_z_vec) - 0.3, max(GPP_z_vec) + 0.3, length.out = n_grid)
GPP_kJ_seq <- GPP_z_to_kJ(GPP_z_seq)
Temp_z_seq <- seq(min(temp_z)    - 0.3, max(temp_z)    + 0.3, length.out = n_grid)
Temp_C_seq <- Temp_z_to_C(Temp_z_seq)

# Observed C data frame (shared; single estimate). Per-variant P observed
# data frames are built inside fit_P_regression().
obs_C_df <- data.frame(
  Stream = streams,
  C_obs  = apply(cons_kJ_m2, 2, median, na.rm = TRUE),
  GPP_kJ = gpp_medians,
  Temp_C = temp_vec
)

# --- Shared helpers --------------------------------------------------------
# Posterior summary table (median, 2.5%, 97.5%)
psumm <- function(draws_obj, pars) {
  t(sapply(pars, function(p) {
    x <- draws_obj$BUGSoutput$sims.list[[p]]
    c(med  = round(median(x), 3),
      lo95 = round(quantile(x, 0.025, names = FALSE), 3),
      hi95 = round(quantile(x, 0.975, names = FALSE), 3))
  }))
}

# lm prediction helper. trans back-transforms the response (exp for log models).
lm_pred <- function(mod, newdata, trans = identity) {
  pr <- predict(mod, newdata = newdata, interval = "confidence")
  data.frame(med = trans(pr[, "fit"]), lo95 = trans(pr[, "lwr"]), hi95 = trans(pr[, "upr"]))
}

# ============================================================
# 7-9. Regression: P (gross & structural) vs GPP + Temp
# ============================================================
# P model (per variant): log(P) ~ GPP + Temp + Temp^2 (<= 0) + GPP:Temp
#   Fit on the LOG scale (both variants strictly positive), mirroring the C model
#   priors known to converge. Wrapped in fit_P_regression() and called once per
#   variant (gross, structural). The C model (below) is fit once.
# n = 8 streams. Priors weakly informative.

# ------------------------------------------------------------
# fit_P_regression(): JAGS + lm fit and all P plots for one variant.
#   P_vec     : named per-stream median production (kJ m-2 day-1), all > 0
#   title_sfx : title annotation, e.g. "(gross)" / "(structural)"
#   file_sfx  : filename suffix, e.g. "_gross" / "_struct"
# Returns a list of fitted objects and plot objects.
# ------------------------------------------------------------
fit_P_regression <- function(P_vec, title_sfx, file_sfx) {
  stopifnot(all(P_vec > 0))

  reg_df_P <- data.frame(
    P_log  = log(P_vec),
    GPP_z  = GPP_z_vec,
    Temp_z = temp_z
  )

  obs_df <- data.frame(
    Stream   = streams,
    P_obs    = P_vec,
    GPP_z    = GPP_z_vec,
    Temp_z   = temp_z,
    Temp_raw = temp_vec
  ) %>% mutate(GPP_kJ = GPP_z_to_kJ(GPP_z))

  ## --- JAGS P model (log scale) ---
  jags_P_str <- "
  model {
    for (i in 1:N) {
      P_log[i] ~ dnorm(mu[i], tau)
      mu[i] <- alpha + beta_GPP * GPP_z[i] + beta_T * Temp_z[i] +
                beta_T2 * pow(Temp_z[i], 2) + beta_int * GPP_z[i] * Temp_z[i]
    }
    # Log-scale priors (mirror the C model; known to converge):
    alpha    ~ dnorm(0, 0.001)
    beta_GPP ~ dnorm(0, 0.01)
    beta_T   ~ dnorm(0, 0.01)
    beta_T2  ~ dnorm(0, 0.01) T(,0)   # constrained: quadratic <= 0
    beta_int ~ dnorm(0, 0.01)
    tau      ~ dgamma(0.001, 0.001)
    sigma    <- 1 / sqrt(tau)
  }
  "
  jags_P_data  <- list(N = nrow(reg_df_P), P_log = reg_df_P$P_log,
                       GPP_z = reg_df_P$GPP_z, Temp_z = reg_df_P$Temp_z)
  jags_P_inits <- function() list(alpha = 0, beta_GPP = 0, beta_T = 0,
                                   beta_T2 = -0.05, beta_int = 0, tau = 1)

  set.seed(42)
  draws_P <- jags(
    data               = jags_P_data,
    inits              = jags_P_inits,
    parameters.to.save = c("alpha", "beta_GPP", "beta_T", "beta_T2", "beta_int", "sigma"),
    model.file         = textConnection(jags_P_str),
    n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 5, DIC = FALSE
  )
  coef_draws <- draws_P$BUGSoutput$sims.list

  cat(sprintf("\n--- JAGS P model %s: log(P) ~ GPP + Temp + Temp^2 + GPP:Temp  [n = 8] ---\n",
              title_sfx))
  cat("NOTE: log link (back-transformed via exp); beta_T2 constrained <= 0.\n\n")
  print(psumm(draws_P, c("alpha", "beta_GPP", "beta_T", "beta_T2", "beta_int", "sigma")))
  cat("Rhat max:", round(max(draws_P$BUGSoutput$summary[, "Rhat"], na.rm = TRUE), 3), "\n")

  ## --- Posterior predictive helpers (log link -> exp back-transform) ---
  jags_pred_P <- function(GPP_z_grid, Temp_z_fixed) {
    base  <- as.vector(coef_draws$alpha) + as.vector(coef_draws$beta_T) * Temp_z_fixed +
              as.vector(coef_draws$beta_T2) * Temp_z_fixed^2
    slope <- as.vector(coef_draws$beta_GPP) + as.vector(coef_draws$beta_int) * Temp_z_fixed
    mu_mat <- sweep(outer(slope, GPP_z_grid), 1, base, "+")
    data.frame(
      GPP_kJ = GPP_z_to_kJ(GPP_z_grid),
      med    = exp(apply(mu_mat, 2, median)),
      lo95   = exp(apply(mu_mat, 2, quantile, 0.025)),
      hi95   = exp(apply(mu_mat, 2, quantile, 0.975))
    )
  }

  jags_pred_P_temp <- function(Temp_z_grid, GPP_z_fixed) {
    base <- as.vector(coef_draws$alpha) + as.vector(coef_draws$beta_GPP) * GPP_z_fixed
    lin  <- as.vector(coef_draws$beta_T) + as.vector(coef_draws$beta_int) * GPP_z_fixed
    quad <- as.vector(coef_draws$beta_T2)
    mu_mat <- sweep(outer(lin, Temp_z_grid) + outer(quad, Temp_z_grid^2), 1, base, "+")
    data.frame(
      Temp_C = Temp_z_to_C(Temp_z_grid),
      med    = exp(apply(mu_mat, 2, median)),
      lo95   = exp(apply(mu_mat, 2, quantile, 0.025)),
      hi95   = exp(apply(mu_mat, 2, quantile, 0.975))
    )
  }

  ## JAGS coefficient forest
  coef_pars_P <- c("alpha", "beta_GPP", "beta_T", "beta_T2", "beta_int")
  coef_summ_P <- as.data.frame(psumm(draws_P, coef_pars_P)) %>%
    tibble::rownames_to_column("term") %>%
    mutate(term = factor(term, levels = rev(coef_pars_P)))

  pJ_coef_P <- ggplot(coef_summ_P, aes(x = med, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lo95, xmax = hi95), linewidth = 0.6) +
    labs(x = "Posterior median (log P per SD)", y = NULL,
         title    = paste("JAGS: log(P) ~ GPP + Temp + Temp\u00B2 + GPP:Temp", title_sfx),
         subtitle = "Median \u00B1 95% CRI; Temp\u00B2 constrained \u2264 0") +
    theme_bw(base_size = 12)

  pred_E_jags <- jags_pred_P(GPP_z_seq, Temp_z_fixed = 0)
  pE_jags <- ggplot(pred_E_jags, aes(x = GPP_kJ)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = med), colour = "steelblue", linewidth = 1) +
    geom_point(data = obs_df, aes(x = GPP_kJ, y = P_obs, colour = Temp_raw), size = 3) +
    geom_text_repel(data = obs_df, aes(x = GPP_kJ, y = P_obs, label = Stream),
                    size = 3.2, colour = "grey30") +
    scale_x_log10() + scale_y_log10() +
    scale_colour_gradient(low = "grey70", high = "tomato", name = "Mean temp (\u00B0C)") +
    labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
         y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
         title    = paste("JAGS: marginal GPP effect on production", title_sfx),
         subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; Temp at mean; shading = 95% CRI") +
    theme_bw(base_size = 12)

  pred_F_jags <- jags_pred_P_temp(Temp_z_seq, GPP_z_fixed = 0)
  pF_jags <- ggplot(pred_F_jags, aes(x = Temp_C)) +
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
         title    = paste("JAGS: quadratic temperature effect on production", title_sfx),
         subtitle = "Log\u2081\u2080 P; GPP at mean; shading = 95% CRI") +
    theme_bw(base_size = 12)

  print(pJ_coef_P)
  print(pE_jags)
  print(pF_jags)
  ggsave(paste0("plot_jags_coefs_P", file_sfx, ".png"),     pJ_coef_P, width = 6, height = 4, dpi = 300)
  ggsave(paste0("plot_jags_effect_GPP", file_sfx, ".png"),  pE_jags,   width = 6, height = 5, dpi = 300)
  ggsave(paste0("plot_jags_effect_Temp", file_sfx, ".png"), pF_jags,   width = 6, height = 5, dpi = 300)

  jags_combined <- (pE_jags + pF_jags) / pJ_coef_P + plot_annotation(
    title = paste("JAGS regression effects", title_sfx),
    theme = theme(plot.title = element_text(size = 13, face = "bold")))
  ggsave(paste0("plot_jags_effects_combined", file_sfx, ".png"),
         jags_combined, width = 13, height = 10, dpi = 300)

  ## --- lm P model (log scale) ---
  P_lm <- lm(P_log ~ GPP_z + Temp_z + GPP_z:Temp_z, data = reg_df_P)

  options(na.action = "na.fail")
  dredge_P <- dredge(P_lm,
                     subset = !(`GPP_z:Temp_z`) | (GPP_z & Temp_z),
                     rank   = "AICc")
  options(na.action = "na.omit")

  cat(sprintf("\n--- lm P model %s: log(P) ~ GPP + Temp + GPP:Temp  [n = 8] ---\n",
              title_sfx))
  cat("NOTE: n = 8 streams, exploratory only.\n\n")
  print(summary(P_lm))

  coef_ci_P <- as.data.frame(cbind(est = coef(P_lm), confint(P_lm)))
  colnames(coef_ci_P) <- c("est", "lo95", "hi95")
  cat("\n--- P model coefficients (estimate, 95% CI) ---\n")
  print(round(coef_ci_P, 3))
  cat("\n--- P model selection (AICc, n = 8) ---\n")
  print(dredge_P)

  coef_df_P <- coef_ci_P %>%
    tibble::rownames_to_column("term") %>%
    mutate(term = factor(term, levels = rev(rownames(coef_ci_P))))

  pC_reg <- ggplot(coef_df_P, aes(x = est, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lo95, xmax = hi95), linewidth = 0.6) +
    labs(x = "Coefficient (log P per SD)", y = NULL,
         title    = paste("lm: log(P) ~ GPP + Temp + GPP:Temp", title_sfx),
         subtitle = "Estimate \u00B1 95% CI (exploratory)") +
    theme_bw(base_size = 12)

  fitted_df <- obs_df %>% mutate(P_fit = exp(fitted(P_lm)))
  pD <- ggplot(fitted_df, aes(x = P_fit, y = P_obs, label = Stream)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = Temp_raw), size = 3) +
    geom_text_repel(size = 3.2, colour = "grey30") +
    scale_x_log10() + scale_y_log10() +
    scale_colour_gradient(low = "steelblue", high = "tomato", name = "Mean temp (\u00B0C)") +
    labs(x = expression("Fitted P (kJ m"^{-2}*" day"^{-1}*")"),
         y = expression("Observed P (kJ m"^{-2}*" day"^{-1}*")"),
         title    = paste("lm: fitted vs observed production", title_sfx),
         subtitle = "Log\u2081\u2080\u2013log\u2081\u2080") +
    theme_bw(base_size = 12)

  pred_E_df <- cbind(GPP_kJ = GPP_kJ_seq,
    lm_pred(P_lm, data.frame(GPP_z = GPP_z_seq, Temp_z = 0), trans = exp))
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
         title    = paste("lm: marginal effect of GPP on production", title_sfx),
         subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; temperature held at its mean") +
    theme_bw(base_size = 12)

  pred_F_df <- cbind(Temp_C = Temp_C_seq,
    lm_pred(P_lm, data.frame(GPP_z = 0, Temp_z = Temp_z_seq), trans = exp))
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
         title    = paste("lm: marginal effect of temperature on production", title_sfx),
         subtitle = "Log\u2081\u2080 P; GPP held at its mean") +
    theme_bw(base_size = 12)

  temp_levels <- c(-1, 0, 1)
  temp_labels <- paste0(round(Temp_z_to_C(temp_levels), 1), " \u00B0C")
  pred_G_df <- bind_rows(lapply(seq_along(temp_levels), function(j) {
    cbind(GPP_kJ = GPP_kJ_seq, Temp_level = temp_labels[j],
          lm_pred(P_lm, data.frame(GPP_z = GPP_z_seq, Temp_z = temp_levels[j]), trans = exp))
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
    scale_colour_manual(values = c("steelblue", "forestgreen", "tomato"), name = "Temperature") +
    scale_fill_manual(values   = c("steelblue", "forestgreen", "tomato"), name = "Temperature") +
    labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
         y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
         title    = paste("lm: GPP effect at low / mean / high temperature", title_sfx),
         subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; shading = 95% CI") +
    theme_bw(base_size = 12)

  gpp_levels <- c(-1, 0, 1)
  gpp_labels <- paste0(round(GPP_z_to_kJ(gpp_levels)), " kJ m\u207B\u00B2")
  pred_H_df <- bind_rows(lapply(seq_along(gpp_levels), function(j) {
    cbind(Temp_C = Temp_C_seq, GPP_level = gpp_labels[j],
          lm_pred(P_lm, data.frame(GPP_z = gpp_levels[j], Temp_z = Temp_z_seq), trans = exp))
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
    scale_colour_manual(values = c("steelblue", "forestgreen", "tomato"), name = "GPP level") +
    scale_fill_manual(values   = c("steelblue", "forestgreen", "tomato"), name = "GPP level") +
    labs(x = "Mean temperature (\u00B0C)",
         y = expression("P (kJ m"^{-2}*" day"^{-1}*")"),
         title    = paste("lm: temperature effect at low / mean / high GPP", title_sfx),
         subtitle = "Log\u2081\u2080 P; shading = 95% CI") +
    theme_bw(base_size = 12)

  print(pC_reg)
  print(pD)
  print(pE)
  print(pF)
  print(pG)
  print(pH)
  ggsave(paste0("plot_lm_coefs_P", file_sfx, ".png"),       pC_reg, width = 6, height = 4, dpi = 300)
  ggsave(paste0("plot_lm_fitted_vs_obs", file_sfx, ".png"), pD,     width = 6, height = 5, dpi = 300)
  ggsave(paste0("plot_lm_effect_GPP", file_sfx, ".png"),    pE,     width = 6, height = 5, dpi = 300)
  ggsave(paste0("plot_lm_effect_Temp", file_sfx, ".png"),   pF,     width = 6, height = 5, dpi = 300)
  ggsave(paste0("plot_lm_GPP_x_Temp", file_sfx, ".png"),    pG,     width = 7, height = 5, dpi = 300)
  ggsave(paste0("plot_lm_Temp_x_GPP", file_sfx, ".png"),    pH,     width = 7, height = 5, dpi = 300)

  lm_combined <- (pE + pF) / (pG + pH) + plot_annotation(
    title = paste("lm regression effects: Fish production ~ GPP + Temp + GPP:Temp", title_sfx),
    theme = theme(plot.title = element_text(size = 13, face = "bold")))
  ggsave(paste0("plot_lm_effects_combined", file_sfx, ".png"),
         lm_combined, width = 13, height = 10, dpi = 300)

  list(draws_P = draws_P, coef_draws = coef_draws, P_lm = P_lm,
       dredge_P = dredge_P, reg_df_P = reg_df_P,
       pJ_coef_P = pJ_coef_P, pE_jags = pE_jags, pF_jags = pF_jags,
       pC_reg = pC_reg, pD = pD, pE = pE, pF = pF, pG = pG, pH = pH)
}

# Fit both variants ----------------------------------------------------------
P_reg_gross  <- fit_P_regression(P_gross_med_vec,  "(gross)",      "_gross")
P_reg_struct <- fit_P_regression(P_struct_med_vec, "(structural)", "_struct")

# ============================================================
# C model: log(C) ~ GPP  (run once)
# ============================================================
# --- JAGS C model ---
jags_C_str <- "
model {
  for (i in 1:N) {
    C_log[i] ~ dnorm(mu[i], tau)
    mu[i]    <- alpha + beta_GPP * GPP_z[i]
  }
  alpha    ~ dnorm(0, 0.001)
  beta_GPP ~ dnorm(0, 0.01) T(0,)   # constrained: C increases with GPP
  tau      ~ dgamma(0.001, 0.001)
  sigma    <- 1 / sqrt(tau)
}
"

jags_C_data  <- list(N = nrow(reg_df), C_log = reg_df$C_log, GPP_z = reg_df$GPP_z)
jags_C_inits <- function() list(alpha = 0, beta_GPP = 0.5, tau = 1)

set.seed(42)
draws_C <- jags(
  data               = jags_C_data,
  inits              = jags_C_inits,
  parameters.to.save = c("alpha", "beta_GPP", "sigma"),
  model.file         = textConnection(jags_C_str),
  n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 5, DIC = FALSE
)
cC_draws <- draws_C$BUGSoutput$sims.list

cat("\n--- JAGS C model: log(C) ~ GPP  [n = 8] ---\n")
cat("NOTE: beta_GPP constrained >= 0 (C cannot decrease with GPP).\n\n")
print(psumm(draws_C, c("alpha", "beta_GPP", "sigma")))
cat("Rhat max:", round(max(draws_C$BUGSoutput$summary[, "Rhat"], na.rm = TRUE), 3), "\n")

jags_pred_C <- function(GPP_z_grid) {
  base  <- as.vector(cC_draws$alpha)
  slope <- as.vector(cC_draws$beta_GPP)
  mu_mat <- sweep(outer(slope, GPP_z_grid), 1, base, "+")
  data.frame(
    GPP_kJ = GPP_z_to_kJ(GPP_z_grid),
    med    = exp(apply(mu_mat, 2, median)),
    lo95   = exp(apply(mu_mat, 2, quantile, 0.025)),
    hi95   = exp(apply(mu_mat, 2, quantile, 0.975))
  )
}

pred_C_jags <- jags_pred_C(GPP_z_seq)
pC_jags <- ggplot(pred_C_jags, aes(x = GPP_kJ)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "forestgreen", alpha = 0.25) +
  geom_line(aes(y = med), colour = "forestgreen", linewidth = 1) +
  geom_point(data = obs_C_df, aes(x = GPP_kJ, y = C_obs, colour = Temp_C), size = 3) +
  geom_text_repel(data = obs_C_df, aes(x = GPP_kJ, y = C_obs, label = Stream),
                  size = 3.2, colour = "grey30") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_gradient(low = "grey70", high = "tomato", name = "Mean temp (\u00B0C)") +
  labs(x = expression("GPP (kJ m"^{-2}*" day"^{-1}*")"),
       y = expression("C (kJ m"^{-2}*" day"^{-1}*")"),
       title    = "JAGS: fish consumption vs GPP",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; beta_GPP \u2265 0; shading = 95% CRI") +
  theme_bw(base_size = 12)
print(pC_jags)
ggsave("plot_jags_C_vs_GPP.png", pC_jags, width = 6, height = 5, dpi = 300)

# --- lm C model (log scale) ---
C_lm <- lm(C_log ~ GPP_z, data = reg_df)

cat("\n--- lm: log(C) ~ log(GPP)  [n = 8] ---\n")
cat("NOTE: n = 8 streams. Exploratory.\n\n")
print(summary(C_lm))

coef_ci_C <- as.data.frame(cbind(est = coef(C_lm), confint(C_lm)))
colnames(coef_ci_C) <- c("est", "lo95", "hi95")
cat("\n--- C model coefficients (estimate, 95% CI) ---\n")
print(round(coef_ci_C, 3))

# C ~ GPP scatter with lm line
pred_C_df <- cbind(GPP_kJ = GPP_kJ_seq,
                   lm_pred(C_lm, data.frame(GPP_z = GPP_z_seq), trans = exp))

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
       title    = "lm: fish consumption vs GPP",
       subtitle = "Log\u2081\u2080\u2013log\u2081\u2080; shading = 95% CI") +
  theme_bw(base_size = 12)

# Coefficient forest — C model
coef_df_C <- coef_ci_C %>%
  tibble::rownames_to_column("term") %>%
  mutate(term = factor(term, levels = rev(rownames(coef_ci_C))))

pC_coef <- ggplot(coef_df_C, aes(x = est, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo95, xmax = hi95), linewidth = 0.6) +
  labs(x = "Coefficient (log C per SD log GPP)", y = NULL,
       title    = "lm: log(C) ~ log(GPP)  [n = 8]",
       subtitle = "Estimate \u00B1 95% CI") +
  theme_bw(base_size = 12)

print(pC_scatter)
print(pC_coef)
ggsave("plot_lm_C_vs_GPP.png", pC_scatter, width = 6, height = 5, dpi = 300)
ggsave("plot_lm_coefs_C.png",  pC_coef,    width = 5, height = 3, dpi = 300)

# ============================================================
# 12. Save combined outputs
# ============================================================

saveRDS(list(
  N_post_d1         = N_post_d1,
  N_post_d2         = N_post_d2,
  prod_gross_kJ_m2  = prod_gross_kJ_m2,
  prod_struct_kJ_m2 = prod_struct_kJ_m2,
  prod_kJ_m2        = prod_gross_kJ_m2,   # back-compat alias (= gross variant)
  resp_kJ_m2        = resp_kJ_m2,
  cons_kJ_m2        = cons_kJ_m2,
  prod_summary      = prod_summary,
  ratio_post        = ratio_post,
  # P regression objects, per variant (log scale)
  P_lm_gross        = P_reg_gross$P_lm,
  P_lm_struct       = P_reg_struct$P_lm,
  coef_draws_gross  = P_reg_gross$coef_draws,   # JAGS P (gross) posterior draws
  coef_draws_struct = P_reg_struct$coef_draws,  # JAGS P (structural) posterior draws
  dredge_gross      = P_reg_gross$dredge_P,
  dredge_struct     = P_reg_struct$dredge_P,
  # C regression objects (single fit)
  C_lm              = C_lm,         # log(C) ~ log(GPP) (lm)
  cC_draws          = cC_draws,     # JAGS C model posterior draws (sims.list)
  streams           = streams,
  stream_area       = stream_area,
  interval_days     = interval_days,
  CALORIC_DENSITY_KJ_PER_G = CALORIC_DENSITY_KJ_PER_G
), "production_results.RDS")

write_csv(prod_summary, "production_summary.csv")

cat("\nDone. Results saved to production_results.RDS and production_summary.csv\n")
