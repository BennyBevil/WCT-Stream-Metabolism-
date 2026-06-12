# 03_Bioenergetics.R
# Westslope cutthroat trout — individual-level bioenergetics using fb4package
# Fish Bioenergetics 4.0 (Deslauriers et al. 2017, fb4package)
# Aug–Oct 2025, 8 isolated streams
#
# Runs a separate fb4 simulation for every Dep1 fish with a measured weight
# (n = 829), producing a per-fish distribution of consumption and respiration.
#
# Final weight (W_final) source per fish:
#   empirical    — PIT-tag recaptures with W measured at Dep2 (n = 153)
#   growth_model — growth model (ΔFL ~ log(FL_aug) + (1|Location)) + L-W
#                  applied to remaining 676 non-recaptured fish
#
# Species: Oncorhynchus clarki (WCT) from fish4_parameters; falls back to RBT.
# Diet: 100% aquatic invertebrate, 2500 cal g-1 wet mass (KEY ASSUMPTION).
# oxycal: 14100 J g-1 O2 (= 14.1 kJ g-1; consistent with Production.R).
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
# 3. Length-weight model (for W_final of non-recaptured fish)
# ============================================================

lw_dat <- dat %>% filter(!is.na(W), !is.na(FL), W > 0, FL > 0)
lw_mod <- lm(log(W) ~ log(FL), data = lw_dat)
lw_a   <- coef(lw_mod)[1]
lw_b   <- coef(lw_mod)[2]
pred_W <- function(FL) exp(lw_a + lw_b * log(FL))

# ============================================================
# 4. Growth model (for W_final of non-recaptured fish)
# ============================================================

tagged_d1 <- dep1 %>%
  filter(!is.na(Tag), Tag != "") %>%
  select(Tag, Location, FL_aug = FL, Date_aug = Date)

recaps_d2 <- dep2 %>%
  filter(Recap == 1) %>%
  select(Tag, Location, FL_oct = FL, Date_oct = Date)

growth_dat <- tagged_d1 %>%
  inner_join(recaps_d2, by = c("Tag", "Location")) %>%
  mutate(dFL = FL_oct - FL_aug) %>%
  filter(!is.na(dFL), !is.na(FL_aug), FL_aug > 0) %>%
  mutate(logFL_aug = log(FL_aug))

growth_mod <- lmer(dFL ~ logFL_aug + (1 | Location),
                   data = growth_dat, REML = TRUE)
fe    <- fixef(growth_mod)
re_df <- ranef(growth_mod)$Location

stream_re <- setNames(
  sapply(streams, function(s) if (s %in% rownames(re_df)) re_df[s, 1] else 0),
  streams
)

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
# W_final: empirical for PIT-tag recaptures; growth model + L-W otherwise.

recaps_w <- dep2 %>%
  filter(Recap == 1, !is.na(W), W > 0) %>%
  select(Tag, Location, W_oct = W) %>%
  distinct(Tag, Location, .keep_all = TRUE)   # keep one recapture per fish

fish_all <- dep1 %>%
  filter(!is.na(W), W > 0, !is.na(FL), FL > 0) %>%
  select(Tag, Location, FL_aug = FL, W_aug = W) %>%
  mutate(fish_id = row_number()) %>%
  left_join(recaps_w, by = c("Tag", "Location")) %>%
  mutate(
    stream_re_val = stream_re[Location],
    mu_dFL        = fe["(Intercept)"] + stream_re_val +
                    fe["logFL_aug"] * log(FL_aug),
    FL_oct_pred   = pmax(FL_aug + mu_dFL, 1),
    W_final_pred  = pred_W(FL_oct_pred),
    W_final       = ifelse(!is.na(W_oct), W_oct, W_final_pred),
    W_final_source = ifelse(!is.na(W_oct), "empirical", "growth_model")
  ) %>%
  left_join(interval_days %>% select(Stream, n_days = days),
            by = c("Location" = "Stream"))

cat("\n--- Fish dataset summary ---\n")
cat("Total fish:", nrow(fish_all), "\n")
fish_all %>%
  count(Location, W_final_source) %>%
  pivot_wider(names_from = W_final_source, values_from = n, values_fill = 0) %>%
  print()

# ============================================================
# 8. Per-fish fb4 simulations
# ============================================================

OXYCAL_J        <- 14100   # J g-1 O2
PREY_ENERGY_CAL <- 2500    # cal g-1 wet mass; KEY ASSUMPTION

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
  prey_energy <- data.frame(Day = seq_len(n_days), invertebrate = PREY_ENERGY_CAL)

  bio_obj <- Bioenergetic(species_params = sp_params, species_info = sp_info)
  bio_obj <- set_environment(bio_obj,        temperature_data = t_df)
  bio_obj <- set_diet(bio_obj,               diet_proportions = diet_prop,
                                              prey_energies    = prey_energy)
  bio_obj <- set_simulation_settings(bio_obj, initial_weight  = fish_i$W_aug,
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
      cum_resp_kJ    = sum(daily$Respiration,       na.rm = TRUE) / 1000,
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
  select(fish_id, Location, FL_aug, W_aug, W_final, W_final_source, n_days) %>%
  left_join(res_df, by = "fish_id") %>%
  rename(Stream = Location)

cat("=== Per-fish simulation summary ===\n")
cat("Converged:", sum(fish_bioen$converged, na.rm = TRUE), "/",
    sum(!is.na(fish_bioen$converged)), "\n\n")

fish_bioen %>%
  group_by(Stream) %>%
  summarise(
    n           = n(),
    n_empirical = sum(W_final_source == "empirical"),
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
  # Supporting models (reused downstream)
  lw_mod          = lw_mod,
  growth_mod      = growth_mod,
  # Constants
  OXYCAL_J        = OXYCAL_J,
  PREY_ENERGY_CAL = PREY_ENERGY_CAL
), "bioenergetics_results.RDS")

cat("\nDone. Results saved to bioenergetics_results.RDS\n")
