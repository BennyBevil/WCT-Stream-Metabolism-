# 02_run_BaMM.R
# For each stream: reads cleaned O2 CSV, writes BaMM.dat, runs the BaMM ADMB
# executable, reads the MCMC output, and saves posteriors as an RDS.
#
# Usage:
#   Rscript 02_run_BaMM.R 1stage     # single-stage respiration model
#   Rscript 02_run_BaMM.R 2stage     # two-stage respiration model
# Defaults to 1stage if no argument is given. The two stages use separate,
# self-contained ADMB working directories (engine_1stage/ and engine_2stage/),
# each with its own BaMM executable, BaMM.cfg, and BaMM.pin — so both stages
# can be run concurrently without sharing/overwriting any config.
#
# Prerequisites:
#   - Run 01_clean_O2.R first to produce Data/O2_Files/cleaned/*_clean.csv
#   - engine_<stage>/ populated with BaMM executable, BaMM.cfg, BaMM.pin
#
# Output RDS naming:
#   1stage -> Data/GPP_EST/{stream}_IND_BaMM.RDS     (consumed downstream)
#   2stage -> Data/GPP_EST/{stream}_2stage_BaMM.RDS
#
# MCMC settings: 20,000 iterations, save every 10th → 2,000 posterior draws.
# Adjust mcmc_iter / mcmc_save below if needed.

library(R2admb)
library(readr)
library(readxl)
library(dplyr)
library(janitor)

# ---- Stage selection --------------------------------------------------------
stage <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(stage)) stage <- "1stage"
if (!stage %in% c("1stage", "2stage", "2stage_priors", "2stage_bounded"))
  stop("Stage must be '1stage', '2stage', '2stage_priors', or '2stage_bounded' ",
       "(got '", stage, "').")
# 2stage_priors  : 2-stage with weakly-informative priors on Eb/Ep/beta_prod
#                  (engine_2stage_priors/). Regularizes respiration only.
# 2stage_bounded : 2-stage with a RECOMPILED engine (engine_2stage_bounded/) in
#                  which log_Pmax/log_I_halfSat are bounded to physical ranges so
#                  flat-diel streams (CC, cow) cannot run the light-curve
#                  parameters off to ~1e23; identifies the Hessian without
#                  changing GPP (the Pmax/I_halfSat ratio is preserved).
# Both write the SAME _2stage_BaMM.RDS suffix so downstream consumers (production,
# model comparison) pick the fit up as the 2-stage result — intended for streams
# whose plain 2-stage fit fails to converge.
rds_suffix <- if (stage == "1stage") "_IND_BaMM.RDS" else "_2stage_BaMM.RDS"

# Optional 2nd arg: comma-separated stream names to restrict the run (else all).
# Useful for re-running only the streams that failed to converge.
only_streams <- commandArgs(trailingOnly = TRUE)[2]
message("=== BaMM run: ", stage, " ===")

# ---- Paths ------------------------------------------------------------------
proj_root      <- "~/Documents/GPP_FP"
bamm_dir       <- file.path(proj_root, "BAMM")
engine_dir     <- file.path(bamm_dir, paste0("engine_", stage))
clean_dir      <- file.path(proj_root, "Data/O2_Files/cleaned")
metrics_xl     <- file.path(proj_root, "Data/Stream.Metrics.xlsx")
good_days_file <- file.path(proj_root, "Data/GoodDays.csv")
out_dir        <- file.path(proj_root, "Data/GPP_EST")
dir.create(out_dir, showWarnings = FALSE)

# Per-day GPP extractor (parses the per-interval P series from BaMM_mcmc.dat
# while it is still fresh in the engine dir, before the next stream overwrites it).
source(file.path(bamm_dir, "extract_daily_gpp.R"))
daily_suffix <- sub("_BaMM\\.RDS$", "_dailyGPP.RDS", rds_suffix)

if (!file.exists(file.path(engine_dir, "BaMM")))
  stop("BaMM executable not found in ", engine_dir,
       " — set up engine_", stage, "/ first.")

# Anchor DOY for Henry (matches doy_start_overrides in 01_clean_O2.R)
henry_anchor <- 227L

# ---- Load GoodDays (selected modelling window per stream) -------------------
good_days <- if (file.exists(good_days_file)) {
  read_csv(good_days_file, show_col_types = FALSE)
} else {
  warning("GoodDays.csv not found — all cleaned data will be passed to BaMM.")
  NULL
}

# ---- MCMC settings ----------------------------------------------------------
mcmc_iter <- 20000
mcmc_save <- 10   # save every nth draw → mcmc_iter / mcmc_save posterior draws

# ---- Site parameters --------------------------------------------------------
# NOTE: BaMM expects slope in degrees (per documentation; default = 0.001).
# Values in Stream.Metrics.xlsx (0.23–3.96) are in degrees and used as-is.
# Constants below match the BaMM defaults from Write.datfile_TJC.R.
site_params <- read_excel(metrics_xl) %>%
  clean_names() %>%
  select(stream, depth_m, elevation_m, latitude, longitude, slope, aspect) %>%
  rename(
    Depth     = depth_m,
    altitude  = elevation_m,
    Latitude  = latitude,
    Longitude = longitude,
    Slope     = slope,
    aspect    = aspect
  ) %>%
  mutate(stream_upper = toupper(stream))   # for case-insensitive joining

# Fixed constants (same for all streams)
DST         <-  -1         # DST offset applied to local time
TimeZone    <-  -6         # Mountain time (UTC-6)
SolarConst  <- 1367
Transmiss   <-   0.8
Area        <-   1
Salinity    <-   0.01
alphaP      <-   1
alphaGK     <-   0.9972
SMeasFlag   <-   0         # 0 = light estimated from solar model (no PAR sensor)
LightSatFlag <-  1
SScale      <-   0.27626576
RKstep      <-   2

# ---- Helper: read MCMC output -----------------------------------------------
read_bamm_mcmc <- function(mcmc_file = "BaMM_mcmc.dat") {
  raw <- read_tsv(mcmc_file, col_names = FALSE,
                  show_col_types = FALSE, progress = FALSE)

  lab_idx     <- seq(3, ncol(raw), by = 2)
  val_idx     <- lab_idx + 1
  param_names <- raw[1, lab_idx] %>% unlist(use.names = FALSE) %>% as.character()
  param_names <- make.names(param_names, unique = TRUE)

  draw_mat <- raw[, val_idx] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  colnames(draw_mat) <- param_names

  as_tibble(draw_mat) %>%
    mutate(iter = as.integer(raw[[1]]),
           obj  = as.numeric(raw[[2]]),
           .before = 1)
}

# ---- Main loop --------------------------------------------------------------
clean_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)
if (length(clean_files) == 0) stop("No cleaned CSVs found in ", clean_dir,
                                    ". Run 01_clean_O2.R first.")

if (!is.na(only_streams)) {
  keep <- trimws(strsplit(only_streams, ",")[[1]])
  clean_files <- clean_files[sub("_clean\\.csv$", "", basename(clean_files)) %in% keep]
  if (length(clean_files) == 0)
    stop("No cleaned CSVs match requested streams: ", only_streams)
  message("Restricting run to: ",
          paste(sub("_clean\\.csv$", "", basename(clean_files)), collapse = ", "))
}

for (f in clean_files) {
  #f <- clean_files[3]
  stream <- sub("_clean\\.csv$", "", basename(f))
  message("\n===== ", stream, " =====")

  # --- Look up site parameters (case-insensitive) ---
  sp <- site_params %>% filter(stream_upper == toupper(stream))
  if (nrow(sp) == 0) {
    message("  No site parameters found for '", stream, "' — skipping.")
    next
  }
  sp <- sp[1, ]   # take first row if somehow duplicated

  # Replace NA site params with safe defaults
  if (is.na(sp$aspect)) sp$aspect <- 0.001
  if (is.na(sp$Slope))  sp$Slope  <- 0.001

  # --- Read cleaned data; drop rows with NA O2 or temperature ---
  dat <- read_csv(f, show_col_types = FALSE) %>%
    filter(!is.na(dissolved_oxygen), !is.na(temperature))

  # --- Filter to GoodDays modelling window ---
  if (!is.null(good_days)) {
    gd <- good_days %>% filter(toupper(Stream) == toupper(stream))
    if (nrow(gd) > 0) {
      udoys <- sort(unique(dat$DOY))
      if (toupper(stream) == "HENRY") {
        start_doy <- henry_anchor + gd$DayStart - 1L
        end_doy   <- henry_anchor + gd$DayEnd   - 1L
      } else {
        start_doy <- udoys[gd$DayStart]
        end_doy   <- udoys[gd$DayEnd]
      }
      dat <- dat %>% filter(DOY >= start_doy, DOY <= end_doy) %>%
        mutate(TimeSinceStart = as.numeric(
          difftime(mountain_standard_time, min(mountain_standard_time), units = "hours")
        ))
      message("  GoodDays window: DOY ", start_doy, "-", end_doy)
    }
  }

  n   <- nrow(dat)
  message("  n = ", n, " obs | DOY ", min(dat$DOY), "-", max(dat$DOY))

  # --- Assemble BaMM data vectors ---
  O2conc        <- dat$dissolved_oxygen
  tC            <- dat$temperature
  irr           <- dat$PAR            # all zeros (no PAR sensor)
  DOY_vec       <- dat$DOY
  TOD_vec       <- dat$DecimalTime
  Interval_vec  <- dat$TimeSinceStart
  n_obs         <- n
  n_d18O        <- 2
  delta_18O_O2  <- c(16, 16)
  DOY_18O       <- c(head(DOY_vec, 1), tail(DOY_vec, 1))
  TOD_18O       <- c(head(TOD_vec, 1), tail(TOD_vec, 1))
  Int_18O       <- c(head(Interval_vec, 1), tail(Interval_vec, 1))
  endofdata     <- -999

  # --- Switch to engine directory, write .dat, run model, then restore wd ---
  orig_wd <- getwd()
  setwd(engine_dir)

  tryCatch({
    # Delete all ADMB-generated files so each stream starts fresh from BaMM.pin.
    # Without this, ADMB reads the previous run's .par file as starting values
    # and BaMM_mcmc.dat is not overwritten by -mcmc (only by -mceval).
    admb_generated <- c(
      "BaMM.par", "BaMM.bar", "BaMM.std", "BaMM.cor", "BaMM.rep",
      "BaMM.psv", "BaMM.mcm", "BaMM.mc2", "BaMM_mcmc.dat",
      "BaMM.b01", "BaMM.b02", "BaMM.b03", "BaMM.b04",
      "BaMM.p01", "BaMM.p02", "BaMM.p03", "BaMM.p04",
      "BaMM.r01", "BaMM.r02", "BaMM.r03", "BaMM.r04",
      "admodel.cov", "admodel.hes", "admodel.dep",
      "admodel.ecm", "admodel.eva", "fmin.log", "BaMM.log"
    )
    file.remove(admb_generated[file.exists(admb_generated)])

    # Write BaMM.dat
    write_dat("BaMM", list(
      altitude     = sp$altitude,
      aspect       = sp$aspect,
      DST          = DST,
      Latitude     = sp$Latitude,
      Longitude    = sp$Longitude,
      Slope        = sp$Slope,
      SolarConst   = SolarConst,
      TimeZone     = TimeZone,
      Transmiss    = Transmiss,
      Area         = Area,
      Depth        = sp$Depth,
      Salinity     = Salinity,
      alphaP       = alphaP,
      alphaGK      = alphaGK,
      SMeasFlag    = SMeasFlag,
      LightSatFlag = LightSatFlag,
      SScale       = SScale,
      RKstep       = RKstep,
      n_irr        = n_obs,
      n_tC         = n_obs,
      n_O2conc     = n_obs,
      n_d18O       = n_d18O,
      irr          = irr,        DOY = DOY_vec, TOD = TOD_vec, Interval = Interval_vec,
      temp         = tC,         DOY = DOY_vec, TOD = TOD_vec, Interval = Interval_vec,
      O2data       = O2conc,     DOY = DOY_vec, TOD = TOD_vec, Interval = Interval_vec,
      O18data      = delta_18O_O2, DOY_18O = DOY_18O, TOD_18O = TOD_18O, Int_18O = Int_18O,
      endofdata    = endofdata
    ))

    # Run MCMC. The default proposal uses the inverse Hessian. For weakly
    # identified 2-stage fits (e.g. CC, cow) the Hessian is not positive
    # definite, so -mcmc aborts before writing BaMM_mcmc.dat. In that case fall
    # back to -mcdiag (identity proposal covariance), which lets MCMC run
    # without inverting the Hessian. Streams that fall back have weakly
    # identified respiration parameters (Eb/Ep/beta_prod) and wide posteriors.
    message("  Running BaMM MCMC (", mcmc_iter, " iters, saving every ", mcmc_save, ")...")
    mcmc_method <- "hessian"
    rc_mcmc <- system(paste0("./BaMM -mcmc ", mcmc_iter, " -mcsave ", mcmc_save))
    rc_eval <- system("./BaMM -mceval")

    if (!file.exists("BaMM_mcmc.dat")) {
      message("  Default (Hessian) MCMC failed — retrying with -mcdiag ",
              "(identity proposal covariance)...")
      mcmc_method <- "mcdiag"
      rc_mcmc <- system(paste0("./BaMM -mcmc ", mcmc_iter,
                               " -mcsave ", mcmc_save, " -mcdiag"))
      rc_eval <- system("./BaMM -mceval")
    }

    if (!file.exists("BaMM_mcmc.dat")) {
      # Capture BaMM output for diagnosis
      bamm_out <- system("./BaMM", intern = TRUE)
      message("  BaMM failed (mcmc rc=", rc_mcmc, " eval rc=", rc_eval,
              "). BaMM stdout:\n", paste(bamm_out, collapse = "\n"))
      stop("BaMM did not produce BaMM_mcmc.dat")
    }

    # Read and save posteriors (record which proposal covariance was used)
    post <- read_bamm_mcmc("BaMM_mcmc.dat")
    attr(post, "mcmc_method") <- mcmc_method
    message("  Posterior draws: ", nrow(post), " (proposal: ", mcmc_method, ")")

    rds_path <- file.path(out_dir, paste0(stream, rds_suffix))
    saveRDS(post, rds_path)
    message("  Saved: ", rds_path)

    # Per-day GPP posterior from the (still-fresh) BaMM_mcmc.dat. `dat` holds the
    # exact time grid (DOY, TimeSinceStart) aligned to the per-interval P series.
    daily <- tryCatch(
      extract_daily_gpp("BaMM_mcmc.dat", dat),
      error = function(e) { message("  daily-GPP extract failed: ", e$message); NULL })
    if (!is.null(daily)) {
      daily_path <- file.path(out_dir, paste0(stream, daily_suffix))
      saveRDS(daily, daily_path)
      message("  Saved daily GPP: ", daily_path,
              " (", nrow(daily$summary), " days; recon-vs-stored max diff ",
              signif(daily$check$max_abs_diff, 2), ")")
    }

  }, error = function(e) {
    message("  ERROR for ", stream, ": ", e$message)
  }, finally = {
    setwd(orig_wd)
  })
}

message("\nAll streams complete (", stage, ").")
