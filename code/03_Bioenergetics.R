# 03_Bioenergetics.R
# Westslope cutthroat trout — individual-level bioenergetics using fb4package
# Fish Bioenergetics 4.0 (Deslauriers et al. 2017, fb4package)
# Aug–Oct 2025, 8 isolated streams
#
# Runs a separate fb4 simulation for every Dep1 fish with measured FL
# (n = 829), producing a per-fish distribution of consumption and respiration.
#
# Growth and length-weight models are fit once in 02_Growth.R (single source
# of truth) and reused here so bioenergetics and 04_Production.R stay
# consistent. The per-day log(dFL) growth model (PIT-tag recaptures) is applied
# uniformly to every Dep1 fish to predict its trajectory, then the hierarchical
# L-W model (stream RE + occasion FE) maps length to weight at each occasion:
#   dFL_day   = exp(b0 + u_stream + b1*log(FL_aug) + sigma^2/2)   [mm day-1]
#   FL_oct    = FL_aug + dFL_day * n_days
#   W_initial = pred_W(FL_aug, stream, "dep1")   Aug L-W relationship
#   W_final   = pred_W(FL_oct, stream, "dep2")   Oct L-W relationship
#
# Species: Oncorhynchus clarki (WCT) from fish4_parameters; falls back to RBT.
# Diet: 100% aquatic invertebrate, 2500 cal g-1 = 10460 J g-1 wet mass (KEY ASSUMPTION).
# oxycal: 14100 J g-1 O2 — passed to run_fb4(); daily_output columns are in
#          J/day (already converted), so post-processing uses /1000 for kJ.
#
# Outputs: bioenergetics_results.RDS
#   fish_bioen  — per-fish data frame: W_initial, W_final, p_value,
#                 cum_resp_kJ, cum_cons_kJ, cum_growth_kJ, W_final_source
#   stream_temp, temp_vec — mean temperature per stream over production interval

library(tidyverse)
library(lubridate)
library(lme4)
library(fb4package)

# ============================================================
# 1. Species lookup in fish4_parameters
# ============================================================

data(fish4_parameters)
sp_names <- names(fish4_parameters)

cat("=== Species available in fish4_parameters ===\n")
cat(paste(sp_names, collapse = "\n"), "\n\n")

wct_match <- grep("clarkii|clarki|westslope|cutthroat",
                  sp_names, ignore.case = TRUE, value = TRUE)
rbt_match <- grep("mykiss|rainbow",
                  sp_names, ignore.case = TRUE, value = TRUE)

if (length(wct_match) > 0) {
  sp_key    <- wct_match[1]
  sp_common <- paste0("cutthroat trout (WCT/YCT): ", sp_key)
  cat("Found cutthroat trout in database:", sp_key, "\n\n")
} else if (length(rbt_match) > 0) {
  sp_key    <- rbt_match[1]
  sp_common <- paste0("rainbow trout (RBT proxy for WCT): ", sp_key)
  message(
    "NOTE: WCT/YCT not found in fish4_parameters; falling back to RBT.\n",
    "  To use WCT-specific params (Hanson et al. 1997: RA=0.00264, RB=-0.217,\n",
    "  RQ=0.06818), override with set_parameter_value() after object creation."
  )
  cat("Using RBT as WCT proxy:", sp_key, "\n\n")
} else {
  stop("No suitable salmonid in fish4_parameters. Available: ",
       paste(sp_names, collapse = ", "))
}

sp_entry  <- fish4_parameters[[sp_key]]
sp_info   <- sp_entry[["species_info"]]
sp_params <- sp_entry[["life_stages"]][["juvenile_adult"]]

# ============================================================
# 2. Load and clean fish data
# ============================================================

dat <- read.csv("Data/CTT.Data.2025.F.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(
    Date     = dmy(Date),
    Tag      = as.character(Tag),
    Location = trimws(Location)
  ) %>%
  filter(!Location %in% c("Jake Canyon", "Alkali", "Pintler"))

streams <- sort(unique(dat$Location))
dep1    <- dat %>% filter(Depletion == 1)
dep2    <- dat %>% filter(Depletion == 2)

# ============================================================
# 3. Growth & length-weight models (centralised in 02_Growth.R)
# ============================================================
# Single source of truth: 02_Growth.R fits the hierarchical L-W model
# (stream RE + occasion fixed effect) and the per-day log(dFL) growth model and
# saves them to growth_results.RDS. Reusing them here guarantees that the
# length->weight mapping and growth projection match 04_Production.R exactly.

grw          <- readRDS("growth_results.RDS")
pred_W       <- grw$pred_W        # pred_W(FL, stream, occasion = "dep1"/"dep2") -> g
fe_FL        <- grw$fe_FL         # per-day log(dFL) fixed effects: (Intercept), logFL_aug
sigma_FL     <- grw$sigma_FL      # residual SD on log scale (log mm day-1)
stream_re_FL <- grw$stream_re_FL  # named vector: stream random intercepts

# ============================================================
# 5. Production interval dates per stream
# ============================================================

interval_days <- dep1 %>%
  group_by(Location) %>%
  summarise(date_d1 = min(Date, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    dep2 %>% group_by(Location) %>%
      summarise(date_d2 = min(Date, na.rm = TRUE), .groups = "drop"),
    by = "Location"
  ) %>%
  mutate(days = as.numeric(date_d2 - date_d1)) %>%
  rename(Stream = Location)

prod_dates <- interval_days %>% mutate(date_d2 = date_d1 + days)

cat("--- Production interval lengths (days) ---\n")
print(interval_days %>% select(Stream, date_d1, days))

# ============================================================
# 6. Temperature data — daily means over production interval
# ============================================================

temp_full <- read.csv("data/temp.csv") %>%
  mutate(
    Date   = as.Date(Date, "%m/%d/%Y"),
    Stream = recode(Stream, CC = "CCT")
  )

daily_temp_prod <- temp_full %>%
  left_join(prod_dates %>% select(Stream, date_d1, date_d2), by = "Stream") %>%
  filter(Date >= date_d1 & Date <= date_d2) %>%
  group_by(Stream, Date) %>%
  summarise(daily_T = mean(Temp, na.rm = TRUE), .groups = "drop")

# Pre-build per-stream temperature data frames (reused across all fish in that stream)
stream_temp_dfs <- lapply(setNames(streams, streams), function(s) {
  daily_temp_prod %>%
    filter(Stream == s) %>%
    arrange(Date) %>%
    mutate(Day = row_number()) %>%
    select(Day, Temperature = daily_T)
})

# ============================================================
# 7. Build per-fish dataset
# ============================================================
# Growth model applied to every Dep1 fish — recaptures inform the model
# but are not used as individual endpoints here.

fish_all <- dep1 %>%
  filter(!is.na(FL), FL > 0) %>%
  select(Tag, Location, FL_aug = FL, W_aug = W) %>%   # W_aug retained for reference
  left_join(interval_days %>% select(Stream, n_days = days),
            by = c("Location" = "Stream")) %>%
  mutate(
    fish_id   = row_number(),
    re_FL     = coalesce(unname(stream_re_FL[Location]), 0),
    dFL_day   = exp(fe_FL["(Intercept)"] + re_FL +
                    fe_FL["logFL_aug"] * log(FL_aug) + sigma_FL^2 / 2),  # mm day-1
    FL_oct    = pmax(FL_aug + dFL_day * n_days, 1),
    W_initial = pred_W(FL_aug, Location, "dep1"),   # Aug L-W relationship
    W_final   = pred_W(FL_oct, Location, "dep2")    # Oct L-W relationship
  )

cat("\n--- Fish dataset summary ---\n")
cat("Total fish:", nrow(fish_all), "\n")
fish_all %>%
  group_by(Location) %>%
  summarise(n = n(), mean_W_initial = round(mean(W_initial), 1),
            mean_W_final = round(mean(W_final), 1), .groups = "drop") %>%
  print()

# ============================================================
# 8. Per-fish fb4 simulations
# ============================================================

OXYCAL_J       <- 14100   # J g-1 O2
PREY_ENERGY_J  <- 2500 * 4.184   # J g-1 wet mass (2500 cal g-1 × 4.184 J cal-1)
                                  # fb4package expects J g-1 throughout

results_list <- vector("list", nrow(fish_all))

cat("\nRunning", nrow(fish_all), "individual fb4 simulations...\n")

for (i in seq_len(nrow(fish_all))) {

  fish_i  <- fish_all[i, ]
  sname   <- fish_i$Location
  n_days  <- fish_i$n_days
  t_df    <- stream_temp_dfs[[sname]]

  if (nrow(t_df) == 0 || is.na(n_days)) {
    results_list[[i]] <- NULL
    next
  }

  diet_prop   <- data.frame(Day = seq_len(n_days), invertebrate = 1.0)
  prey_energy <- data.frame(Day = seq_len(n_days), invertebrate = PREY_ENERGY_J)

  bio_obj <- Bioenergetic(species_params = sp_params, species_info = sp_info)
  bio_obj <- set_environment(bio_obj,        temperature_data = t_df)
  bio_obj <- set_diet(bio_obj,               diet_proportions = diet_prop,
                                              prey_energies    = prey_energy)
  bio_obj <- set_simulation_settings(bio_obj, initial_weight  = fish_i$W_initial,
                                               duration        = n_days)

  res_i <- tryCatch(
    run_fb4(bio_obj,
            fit_to    = "Weight",
            fit_value = fish_i$W_final,
            strategy  = "binary_search",
            oxycal    = OXYCAL_J,
            verbose   = FALSE),
    error = function(e) NULL
  )

  if (!is.null(res_i)) {
    daily <- res_i$daily_output
    results_list[[i]] <- list(
      fish_id        = fish_i$fish_id,
      p_value        = res_i$summary$p_value,
      converged      = res_i$summary$converged,
      cum_resp_kJ    = sum(daily$Respiration,        na.rm = TRUE) / 1000,
      cum_cons_kJ    = sum(daily$Consumption_energy, na.rm = TRUE) / 1000,
      cum_growth_kJ  = sum(daily$Net_energy,         na.rm = TRUE) / 1000
    )
  }

  if (i %% 100 == 0) cat("  Completed", i, "/", nrow(fish_all), "\n")
}

cat("  Done.\n\n")

# ============================================================
# 9. Assemble per-fish results data frame
# ============================================================

res_df <- bind_rows(lapply(results_list, as.data.frame))

fish_bioen <- fish_all %>%
  select(fish_id, Location, FL_aug, W_aug, W_initial, FL_oct, W_final, n_days) %>%
  left_join(res_df, by = "fish_id") %>%
  rename(Stream = Location)

cat("=== Per-fish simulation summary ===\n")
cat("Converged:", sum(fish_bioen$converged, na.rm = TRUE), "/",
    sum(!is.na(fish_bioen$converged)), "\n\n")

fish_bioen %>%
  group_by(Stream) %>%
  summarise(
    n           = n(),
    p_median    = round(median(p_value,     na.rm = TRUE), 3),
    p_lo95      = round(quantile(p_value,   0.025, na.rm = TRUE), 3),
    p_hi95      = round(quantile(p_value,   0.975, na.rm = TRUE), 3),
    cons_med_kJ = round(median(cum_cons_kJ, na.rm = TRUE), 2),
    resp_med_kJ = round(median(cum_resp_kJ, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

# ============================================================
# 10. Mean temperature per stream (used in downstream regression)
# ============================================================

stream_temp_sum <- daily_temp_prod %>%
  group_by(Stream) %>%
  summarise(mean_temp = mean(daily_T, na.rm = TRUE), .groups = "drop")

temp_vec <- stream_temp_sum %>%
  arrange(match(Stream, streams)) %>%
  pull(mean_temp)

cat("\n--- Mean temperature over production interval (°C) ---\n")
print(stream_temp_sum)

# ============================================================
# 11. Save
# ============================================================

saveRDS(list(
  # Per-fish results (main output)
  fish_bioen      = fish_bioen,
  # Species used
  species_key     = sp_key,
  species_common  = sp_common,
  stream_temp     = stream_temp_sum,
  temp_vec        = temp_vec,
  prod_dates      = prod_dates,
  daily_temp_prod = daily_temp_prod,
  interval_days   = interval_days,
  # Supporting models (centralised in 02_Growth.R; kept here for provenance)
  lw_mod          = grw$lw_mod,
  growth_mod      = grw$growth_mod_FL,
  # Constants
  OXYCAL_J      = OXYCAL_J,
  PREY_ENERGY_J = PREY_ENERGY_J
), "bioenergetics_results.RDS")

cat("\nDone. Results saved to bioenergetics_results.RDS\n")
