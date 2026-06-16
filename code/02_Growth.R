# 02_Growth.R
# Westslope cutthroat trout — length-weight and growth models
# Aug–Oct 2025, 8 isolated streams
#
# Models:
#   L-W:     log(W) ~ log(FL) + occasion + (1 | Location)  — hierarchical
#            Common allometric slope; occasion (Dep1/Dep2) FIXED-effect
#            intercept shift (Aug vs Oct condition); stream random intercept.
#            W (g) = exp(b0 + b1*log(FL) + b_occ*[occ==dep2] + u_stream)
#   Growth:  dFL_day ~ log(FL_aug) + (1 | Location)  — mm day-1
#            dW_day  ~ log(W_aug)  + (1 | Location)  — g day-1
#   Both fit by REML on PIT-tagged Dep1 fish recaptured at Dep2.
#   Using FL_aug as the size covariate for both (common reference; all
#   dep1 fish have FL, not all have W).
#
# NOTE: Production.R currently predicts total-interval dFL and divides by
#   days_s in the MC loop. Once you switch to per-day growth models here,
#   that loop will need updating: multiply predicted daily rate by days_s
#   before computing FL_oct and dW.
#
# Outputs: growth_results.RDS

library(tidyverse)
library(lubridate)
library(lme4)
library(patchwork)

# ============================================================
# 1. Load and clean data
# ============================================================

dat <- read.csv("Data/CTT.Data.2025.F.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(
    Date     = dmy(Date),
    Tag      = as.character(Tag),
    Location = trimws(Location)
  ) %>%
  filter(!Location %in% c("Jake Canyon", "Alkali", "Pintler"))

streams <- sort(unique(dat$Location))

dep1 <- dat %>% filter(Depletion == 1)
dep2 <- dat %>% filter(Depletion == 2)

# ============================================================
# 2. Length-weight model (hierarchical)
# ============================================================
# Fit on all fish with both W and FL recorded.
#   log(W) = b0 + b1*log(FL) + b_occ*[occasion==dep2] + u_stream
# - Common allometric slope b1 across streams/occasions.
# - occasion (Dep1 Aug vs Dep2 Oct) enters as a FIXED-effect intercept shift:
#   a 2-level grouping factor cannot support a reliable variance component,
#   so a random effect would be near-singular. The fixed effect still gives a
#   distinct Aug vs Oct relationship (occasion-accurate predicted weight).
# - Stream (8 levels) is a random intercept.

lw_dat <- dat %>%
  filter(!is.na(W), !is.na(FL), W > 0, FL > 0) %>%
  mutate(
    occasion = factor(if_else(Depletion == 1, "dep1", "dep2"),
                      levels = c("dep1", "dep2")),
    logFL    = log(FL)
  )

lw_mod <- lmer(log(W) ~ logFL + occasion + (1 | Location),
               data = lw_dat, REML = TRUE)

lw_fe    <- fixef(lw_mod)
lw_sigma <- sigma(lw_mod)
lw_re    <- ranef(lw_mod)$Location
lw_stream_re <- setNames(
  sapply(streams, function(s) if (s %in% rownames(lw_re)) lw_re[s, 1] else 0),
  streams
)

# Back-compat population-average scalars (dep1 reference, stream mean RE = 0).
# Used by 05_Growth_Uncertainty.R for its chain-rule sensitivity calc.
lw_a <- unname(lw_fe["(Intercept)"])
lw_b <- unname(lw_fe["logFL"])

# Portable occasion- and stream-aware predictor (median weight, g).
# Closes over fixed effects + stream REs only, so it serialises cleanly with
# growth_results.RDS; 03_Bioenergetics.R and 04_Production.R both call this so
# weights are guaranteed consistent across the pipeline.
pred_W <- local({
  b0  <- unname(lw_fe["(Intercept)"])
  b1  <- unname(lw_fe["logFL"])
  bo2 <- {
    nm <- grep("^occasion", names(lw_fe), value = TRUE)
    if (length(nm)) unname(lw_fe[nm[1]]) else 0
  }
  sre <- lw_stream_re
  function(FL, stream, occasion = "dep1") {
    re <- unname(sre[as.character(stream)]); re[is.na(re)] <- 0
    oc <- ifelse(as.character(occasion) == "dep2", bo2, 0)
    exp(b0 + b1 * log(FL) + oc + re)
  }
})

cat("\n--- Length-weight model (hierarchical: stream RE + occasion FE) ---\n")
print(summary(lw_mod))

# ============================================================
# 3. Growth data: match Dep1 tagged fish to Dep2 recaptures
# ============================================================

tagged_d1 <- dep1 %>%
  filter(!is.na(Tag), Tag != "") %>%
  select(Tag, Location, FL_aug = FL, W_aug = W, Date_aug = Date)

recaps_d2 <- dep2 %>%
  filter(Recap == 1) %>%
  select(Tag, Location, FL_oct = FL, W_oct = W, Date_oct = Date)

growth_dat <- tagged_d1 %>%
  inner_join(recaps_d2, by = c("Tag", "Location")) %>%
  mutate(
    days    = as.numeric(Date_oct - Date_aug),
    dFL     = FL_oct - FL_aug,
    dW      = W_oct  - W_aug,
    dFL_day = dFL / days,
    dW_day  = dW  / days
  ) %>%
  filter(
    !is.na(dFL), !is.na(FL_aug), FL_aug > 0,
    days > 0
  )

cat("\n--- Growth data summary (all matched pairs with FL) ---\n")
growth_dat %>%
  group_by(Location) %>%
  summarise(
    n          = n(),
    n_W        = sum(!is.na(dW) & W_aug > 0 & W_oct > 0),
    n_neg_dFL  = sum(dFL_day <= 0),
    n_neg_dW   = sum(dW_day  <= 0, na.rm = TRUE),
    FL_range   = paste0(round(min(FL_aug)), "\u2013", round(max(FL_aug))),
    dFL_day_mn = round(mean(dFL_day), 3),
    dFL_day_sd = round(sd(dFL_day),   3),
    dW_day_mn  = round(mean(dW_day, na.rm = TRUE), 4),
    dW_day_sd  = round(sd(dW_day,  na.rm = TRUE), 4)
  ) %>%
  print()

cat("\nFL range in Dep1 fish (all, not just tagged):",
    min(dep1$FL, na.rm = TRUE), "\u2013", max(dep1$FL, na.rm = TRUE), "mm\n")
cat("FL range in tagged/recaptured fish:",
    min(growth_dat$FL_aug), "\u2013", max(growth_dat$FL_aug), "mm\n")
cat("NOTE: growth predictions for fish below the tagging threshold are extrapolated.\n")
cat("NOTE: fish with dFL_day <= 0 or dW_day <= 0 excluded from log-scale models (see n_neg above).\n\n")

# Model-fitting subsets: log transform requires strictly positive growth.
# growth_dat kept whole for growth tracks plot (raw FL values).
growth_dat_FL <- growth_dat %>%
  filter(dFL_day > 0) %>%
  mutate(logFL_aug    = log(FL_aug),
         log_dFL_day  = log(dFL_day))

growth_dat_W <- growth_dat %>%
  filter(!is.na(dW), W_aug > 0, W_oct > 0, dW_day > 0) %>%
  mutate(logW_aug    = log(W_aug),
         log_dW_day  = log(dW_day))

# ============================================================
# 4. Growth models
# ============================================================

# --- 4a. log(dFL_day): log mm per day ---
growth_mod_FL <- lmer(log_dFL_day ~ logFL_aug + (1 | Location),
                      data = growth_dat_FL, REML = TRUE)
cat(sprintf("--- Growth model: log(dFL_day)  [n = %d; %d excluded, dFL <= 0] ---\n",
            nrow(growth_dat_FL),
            nrow(growth_dat) - nrow(growth_dat_FL)))
print(summary(growth_mod_FL))

fe_FL      <- fixef(growth_mod_FL)
re_FL      <- ranef(growth_mod_FL)$Location
sigma_FL   <- sigma(growth_mod_FL)

stream_re_FL <- setNames(
  sapply(streams, function(s) {
    if (s %in% rownames(re_FL)) re_FL[s, 1] else 0
  }),
  streams
)
cat("\n--- Stream random effects: log(dFL_day) intercept adjustments (log mm day-1) ---\n")
print(round(stream_re_FL, 4))

# --- 4b. log(dW_day): log g per day ---
growth_mod_W <- lmer(log_dW_day ~ logW_aug + (1 | Location),
                     data = growth_dat_W, REML = TRUE)
n_W_all <- sum(!is.na(growth_dat$dW) & growth_dat$W_aug > 0 & growth_dat$W_oct > 0)
cat(sprintf("\n--- Growth model: log(dW_day)  [n = %d; %d excluded, dW <= 0 or missing W] ---\n",
            nrow(growth_dat_W),
            n_W_all - nrow(growth_dat_W)))
print(summary(growth_mod_W))

fe_W       <- fixef(growth_mod_W)
re_W       <- ranef(growth_mod_W)$Location
sigma_W    <- sigma(growth_mod_W)

stream_re_W <- setNames(
  sapply(streams, function(s) {
    if (s %in% rownames(re_W)) re_W[s, 1] else 0
  }),
  streams
)
cat("\n--- Stream random effects: log(dW_day) intercept adjustments (log g day-1) ---\n")
print(round(stream_re_W, 5))

# ============================================================
# 5. Stream-level growth summary
# ============================================================

# --- 5a. Summary table (natural scale for interpretability) ---
summ_FL <- growth_dat_FL %>%
  group_by(Location) %>%
  summarise(
    n_FL       = n(),
    dFL_median = round(median(dFL_day), 3),
    dFL_mean   = round(mean(dFL_day),   3),
    dFL_sd     = round(sd(dFL_day),     3),
    .groups    = "drop"
  )

summ_W <- growth_dat_W %>%
  group_by(Location) %>%
  summarise(
    n_W       = n(),
    dW_median = round(median(dW_day), 4),
    dW_mean   = round(mean(dW_day),   4),
    dW_sd     = round(sd(dW_day),     4),
    .groups   = "drop"
  )

stream_growth_summ <- full_join(summ_FL, summ_W, by = "Location") %>%
  arrange(desc(dFL_median))

cat("\n--- Stream-level growth summary (natural scale, per day) ---\n")
cat("dFL: mm day-1  |  dW: g day-1\n\n")
print(stream_growth_summ)

# --- 5b. Per-stream violin + pointrange (natural scale) ---
# Streams ordered by median dFL_day for both panels so they align visually.
stream_order_g <- summ_FL %>% arrange(dFL_median) %>% pull(Location)

ps_FL <- growth_dat_FL %>%
  mutate(Location = factor(Location, levels = stream_order_g)) %>%
  group_by(Location) %>%
  mutate(med = median(dFL_day), lo = quantile(dFL_day, 0.25),
         hi  = quantile(dFL_day, 0.75)) %>%
  ungroup()

p_stream_FL <- ggplot(ps_FL, aes(x = dFL_day, y = Location)) +
  geom_violin(fill = "tomato", alpha = 0.45, colour = NA) +
  geom_pointrange(aes(x = med, xmin = lo, xmax = hi),
                  linewidth = 0.7, size = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = expression(Delta*"FL (mm day"^{-1}*")"), y = NULL,
       title = "Fork-length growth by stream") +
  theme_bw(base_size = 11)

ps_W <- growth_dat_W %>%
  mutate(Location = factor(Location, levels = stream_order_g)) %>%
  group_by(Location) %>%
  mutate(med = median(dW_day), lo = quantile(dW_day, 0.25),
         hi  = quantile(dW_day, 0.75)) %>%
  ungroup()

p_stream_W <- ggplot(ps_W, aes(x = dW_day, y = Location)) +
  geom_violin(fill = "steelblue", alpha = 0.45, colour = NA) +
  geom_pointrange(aes(x = med, xmin = lo, xmax = hi),
                  linewidth = 0.7, size = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = expression(Delta*"W (g day"^{-1}*")"), y = NULL,
       title = "Weight growth by stream") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

print(p_stream_FL | p_stream_W)

# --- 5c. Random effects forest plot (log scale — model estimates) ---
# condVar = TRUE gives approximate conditional SEs for each stream RE.
# Intervals are ±2 conditional SD (approx 95%).

re_FL_df <- as.data.frame(ranef(growth_mod_FL, condVar = TRUE)) %>%
  filter(term == "(Intercept)") %>%
  transmute(Location = grp, condval, condsd, Model = "log(\u0394FL day\u207B\u00B9)")

re_W_df <- as.data.frame(ranef(growth_mod_W, condVar = TRUE)) %>%
  filter(term == "(Intercept)") %>%
  transmute(Location = grp, condval, condsd, Model = "log(\u0394W day\u207B\u00B9)")

re_df <- bind_rows(re_FL_df, re_W_df) %>%
  mutate(
    lo95 = condval - 2 * condsd,
    hi95 = condval + 2 * condsd,
    Location = factor(Location, levels = stream_order_g),
    Model    = factor(Model, levels = c("log(\u0394FL day\u207B\u00B9)",
                                        "log(\u0394W day\u207B\u00B9)"))
  )

p_re <- ggplot(re_df, aes(x = condval, y = Location, colour = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo95, xmax = hi95),
                  position = position_dodge(width = 0.5),
                  linewidth = 0.6, size = 0.4) +
  scale_colour_manual(values = c("log(\u0394FL day\u207B\u00B9)" = "tomato",
                                 "log(\u0394W day\u207B\u00B9)"  = "steelblue")) +
  labs(
    x      = "Random effect (log scale, deviation from population mean)",
    y      = NULL, colour = NULL,
    title  = "Stream-level random effects \u00B1 2 conditional SD",
    subtitle = "Both models on log scale; streams ordered by median \u0394FL"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

print(p_re)

# ============================================================
# 6. Length-frequency diagnostics: Aug vs Oct
# ============================================================
# Uses all fish with valid FL — not just tagged/recaptured — to represent
# the full captured population at each occasion.

lf_dat <- bind_rows(
  dep1 %>% filter(!is.na(FL), FL > 0) %>%
    mutate(Occasion = "Dep1 (Aug)"),
  dep2 %>% filter(!is.na(FL), FL > 0) %>%
    mutate(Occasion = "Dep2 (Oct)")
) %>%
  mutate(Occasion = factor(Occasion, levels = c("Dep1 (Aug)", "Dep2 (Oct)")))

# --- 5a. Summary table ---
cat("\n--- Length-frequency summary by stream and occasion ---\n")
lf_dat %>%
  group_by(Location, Occasion) %>%
  summarise(
    n      = n(),
    FL_min = round(min(FL), 0),
    FL_med = round(median(FL), 0),
    FL_max = round(max(FL), 0),
    FL_mn  = round(mean(FL), 1),
    FL_sd  = round(sd(FL),   1),
    .groups = "drop"
  ) %>%
  arrange(Location, Occasion) %>%
  print()

# --- 5b. Faceted length-frequency histograms ---
# Bin width fixed at 10 mm so panels are directly comparable across streams.
p_lf <- ggplot(lf_dat, aes(x = FL, fill = Occasion, colour = Occasion)) +
  geom_histogram(binwidth = 10, position = "identity", alpha = 0.5,
                 boundary = 0) +
  geom_vline(xintercept = 80, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  facet_wrap(~ Location, ncol = 2, scales = "free_y") +
  scale_fill_manual(values  = c("Dep1 (Aug)" = "tomato",
                                "Dep2 (Oct)" = "steelblue")) +
  scale_colour_manual(values = c("Dep1 (Aug)" = "tomato",
                                 "Dep2 (Oct)" = "steelblue")) +
  labs(
    x      = "Fork length (mm)",
    y      = "Count",
    fill   = NULL, colour = NULL,
    title  = "Length-frequency distributions: Dep1 (Aug) vs Dep2 (Oct)",
    subtitle = "All fish with valid FL; bin width = 10 mm; y-axis free per stream"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

print(p_lf)

# --- 5c. Individual growth tracks for PIT-tagged recaptured fish ---
# Each line connects FL_aug to FL_oct for one fish; coloured by stream.
# Helps visualise whether the recaptured size range is representative and
# whether any fish show implausible shrinkage.

tracks <- growth_dat %>%
  select(Tag, Location, FL_aug, FL_oct, dFL_day) %>%
  pivot_longer(cols = c(FL_aug, FL_oct),
               names_to  = "Occasion",
               values_to = "FL") %>%
  mutate(Occasion = recode(Occasion,
                           FL_aug = "Dep1 (Aug)",
                           FL_oct = "Dep2 (Oct)"),
         Occasion = factor(Occasion, levels = c("Dep1 (Aug)", "Dep2 (Oct)")))

p_tracks <- ggplot(tracks, aes(x = Occasion, y = FL,
                                group = Tag, colour = Location)) +
  geom_line(alpha = 0.45, linewidth = 0.5) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ Location, ncol = 2) +
  labs(
    x      = NULL,
    y      = "Fork length (mm)",
    colour = NULL,
    title  = "Individual growth tracks: PIT-tagged recaptures",
    subtitle = "Each line = one fish; negative slopes = apparent shrinkage"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

print(p_tracks)

# ============================================================
# 7. Model diagnostics and fit plots
# ============================================================

growth_dat_FL$pred_log_dFL_day <- predict(growth_mod_FL)
growth_dat_W$pred_log_dW_day   <- predict(growth_mod_W)

# Helper: three-panel diagnostic for a fitted lmer
diag_plots <- function(df, obs_col, pred_col, loc_col = "Location",
                       ylab, title_prefix) {

  resid_col <- paste0(obs_col, "_resid")
  df[[resid_col]] <- df[[obs_col]] - df[[pred_col]]

  p_fit <- ggplot(df, aes(x = .data[[pred_col]], y = .data[[obs_col]],
                           colour = .data[[loc_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7) +
    labs(x = paste0("Predicted ", ylab), y = paste0("Observed ", ylab),
         title = paste0(title_prefix, ": fitted vs observed"),
         colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  p_res <- ggplot(df, aes(x = .data[[pred_col]], y = .data[[resid_col]],
                           colour = .data[[loc_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7) +
    labs(x = paste0("Predicted ", ylab), y = paste0("Residual ", ylab),
         title = paste0(title_prefix, ": residuals vs fitted"),
         colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  p_qq <- ggplot(df, aes(sample = .data[[resid_col]])) +
    stat_qq(alpha = 0.6) +
    stat_qq_line(colour = "steelblue") +
    labs(x = "Theoretical quantile", y = paste0("Sample quantile (", ylab, ")"),
         title = paste0(title_prefix, ": Q-Q residuals")) +
    theme_bw(base_size = 11)

  list(fit = p_fit, res = p_res, qq = p_qq)
}

diag_FL <- diag_plots(growth_dat_FL, "log_dFL_day", "pred_log_dFL_day",
                      ylab = "log(dFL) (log mm day\u207B\u00B9)", title_prefix = "\u0394FL")
diag_W  <- diag_plots(growth_dat_W,  "log_dW_day",  "pred_log_dW_day",
                      ylab = "log(dW) (log g day\u207B\u00B9)",   title_prefix = "\u0394W")

# Model fit plots: log growth ~ logFL coloured by stream
p_mod_FL <- ggplot(growth_dat_FL,
                   aes(x = logFL_aug, y = log_dFL_day, colour = Location)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = pred_log_dFL_day, group = Location), linewidth = 0.8) +
  labs(
    x     = "log(Fork length at tagging, mm)",
    y     = expression("log("*Delta*"FL) (log mm day"^{-1}*")"),
    title = "log(\u0394FL) model: log(dFL_day) ~ log(FL_aug) + (1|Location)",
    colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

p_mod_W <- ggplot(growth_dat_W,
                  aes(x = logW_aug, y = log_dW_day, colour = Location)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = pred_log_dW_day, group = Location), linewidth = 0.8) +
  labs(
    x     = "log(Weight at tagging, g)",
    y     = expression("log("*Delta*"W) (log g day"^{-1}*")"),
    title = "log(\u0394W) model: log(dW_day) ~ log(W_aug) + (1|Location)",
    colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

lw_occ_shift <- { nm <- grep("^occasion", names(lw_fe), value = TRUE)
                  if (length(nm)) unname(lw_fe[nm[1]]) else 0 }
p_lw <- ggplot(lw_dat, aes(x = log(FL), y = log(W), colour = occasion)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = c(dep1 = "tomato", dep2 = "steelblue")) +
  labs(
    x = "log(Fork length, mm)",
    y = "log(Weight, g)",
    colour = "Occasion",
    title = sprintf(
      "Length-weight (hierarchical): log(W) = %.3f + %.3f \u00D7 log(FL); Dep2 shift = %.3f",
      lw_a, lw_b, lw_occ_shift)
  ) +
  theme_bw(base_size = 11)

# Print figures
print(p_lw)
print(p_mod_FL)
print(p_mod_W)

# Parallel diagnostic panels
(diag_FL$fit | diag_FL$res | diag_FL$qq) /
(diag_W$fit  | diag_W$res  | diag_W$qq) +
  plot_annotation(title = "Growth model diagnostics: \u0394FL (top) and \u0394W (bottom)",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

# ============================================================
# 8. Save
# ============================================================

saveRDS(list(
  # Length-weight (hierarchical: stream RE + occasion FE)
  lw_mod       = lw_mod,
  lw_fe        = lw_fe,
  lw_sigma     = lw_sigma,
  lw_stream_re = lw_stream_re,
  pred_W       = pred_W,      # pred_W(FL, stream, occasion = "dep1"/"dep2") -> g
  lw_a         = lw_a,        # back-compat scalars (population-average, dep1)
  lw_b         = lw_b,
  # dFL model (mm day-1)
  growth_mod_FL = growth_mod_FL,
  fe_FL         = fe_FL,
  sigma_FL      = sigma_FL,
  stream_re_FL  = stream_re_FL,
  # dW model (g day-1)
  growth_mod_W  = growth_mod_W,
  fe_W          = fe_W,
  sigma_W       = sigma_W,
  stream_re_W   = stream_re_W,
  # Data (growth_dat = all matched pairs, used for tracks plot)
  growth_dat    = growth_dat,
  growth_dat_FL = growth_dat_FL,
  growth_dat_W  = growth_dat_W
), "growth_results.RDS")

cat("\nDone. Results saved to growth_results.RDS\n")

