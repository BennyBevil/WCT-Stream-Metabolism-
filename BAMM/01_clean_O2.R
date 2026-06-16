# 01_clean_O2.R
# Reads each MiniDOT .Cat.TXT file, cleans it, and writes a cleaned CSV.
# Output: Data/O2_Files/cleaned/{stream}_clean.csv
#
# Cleaning steps:
#   1. Parse Mountain Standard Time timestamps
#   2. Compute DOY, DecimalTime, DayFrac, TimeSinceStart, PAR (= 0, no sensor)
#   3. Retain first 28 unique DOYs (full context window; GoodDays selection applied in 02_run_BaMM.R)
#   4. Q score filter: dissolved_oxygen / saturation → NA where q < q_threshold
#   5. Per-stream DOY exclusions (manual)
#   6. Per-stream point exclusions (manual)
#   7. Automated QC:
#        a. Physical bounds (DO, saturation, temperature)
#        b. Flatline detection (consecutive identical DO values)
#        c. Hampel spike filter (rolling median ± k·MAD)
#        d. Step-change filter (large 10-min DO jump)
#      Flagged DO and saturation values are set to NA; temperature is only
#      flagged by physical bounds. BaMM run script drops rows with NA dissolved_oxygen.

library(readr)
library(dplyr)
library(lubridate)
library(janitor)

# ---- Paths ------------------------------------------------------------------
proj_root <- "~/Documents/GPP_FP"
raw_dir   <- file.path(proj_root, "Data/O2_Files")
out_dir   <- file.path(raw_dir, "cleaned")
dir.create(out_dir, showWarnings = FALSE)

txt_files <- list.files(raw_dir, pattern = "\\.TXT$", full.names = TRUE)
if (length(txt_files) == 0) stop("No .TXT files found in ", raw_dir)

# ---- Q score threshold -------------------------------------------------------
# PME MiniDOT Q reflects optical foil condition; PME threshold is 0.7 (serious
# problem below). Using 0.8 as a conservative threshold. Observations below this
# have dissolved_oxygen and dissolved_oxygen_saturation set to NA; temperature
# is retained (thermistor is independent). BaMM run script drops NA O2 rows.
# No effect on 2025 dataset (all streams q >= 0.857) but retained as a safeguard.
q_threshold <- 0.8

# ---- Automated QC thresholds ------------------------------------------------
# 1. Physical bounds — values outside these ranges are physically implausible
do_bounds   <- c(0, 20)     # mg/L
sat_bounds  <- c(0, 150)    # %
temp_bounds <- c(0, 30)     # °C

# 2. Flatline detection — runs of >= this many consecutive identical DO values
flatline_n  <- 3

# 3. Hampel spike filter — rolling window of ±hampel_half_w observations (each
#    side); flag if |value - window median| > hampel_k × window MAD
hampel_half_w <- 6     # 6 obs × 10 min = ±1 hr window on each side
hampel_k      <- 3

# 4. Step-change threshold — flag 10-min DO steps larger than this
step_threshold <- 1.0  # mg/L

# ---- Per-stream DOY overrides -----------------------------------------------
# Named list: stream name → first DOY to keep (14 days will be retained from there).
# Streams not listed use the default: unique DOYs[2:15] (skip deployment day 1).
doy_start_overrides <- list(
  Henry = 227   # Aug 15; long time series, anchor to deployment window
)

# ---- Per-stream DOY exclusions ----------------------------------------------
# Named list: stream name → vector of DOYs to drop entirely after initial selection.
doy_exclude <- list(
  CC = 226   # Rapid DO oscillations throughout the day; likely sensor disturbance
)

# ---- Per-stream point exclusions --------------------------------------------
# Named list: stream name → data frame of DOY/Hour/Minute rows to drop.
point_exclude <- list(
  Buf = data.frame(DOY = 223, Hour = 5, Minute = 32)  # 1.73 mg/L step; anomalous single reading
)

# ---- Automated QC function --------------------------------------------------
# Applies checks 1-4 to dissolved_oxygen (and saturation). Temperature is only
# checked against physical bounds. Returns the data with flagged values as NA
# and prints a per-check summary.
apply_auto_qc <- function(dat, stream_name) {

  do   <- dat$dissolved_oxygen
  sat  <- dat$dissolved_oxygen_saturation
  temp <- dat$temperature
  n    <- nrow(dat)

  flag_do <- rep(FALSE, n)   # combined DO flag
  counts  <- integer(4)      # flags per check (for reporting)

  # --- 1. Physical bounds ---
  b_do   <- !is.na(do)   & (do   < do_bounds[1]   | do   > do_bounds[2])
  b_sat  <- !is.na(sat)  & (sat  < sat_bounds[1]  | sat  > sat_bounds[2])
  b_temp <- !is.na(temp) & (temp < temp_bounds[1] | temp > temp_bounds[2])
  flag_do   <- flag_do | b_do | b_sat
  counts[1] <- sum(b_do | b_sat | b_temp, na.rm = TRUE)

  # --- 2. Flatline detection ---
  do_round    <- round(do, 3)
  run_lengths <- rep(rle(do_round)$lengths, rle(do_round)$lengths)
  fl          <- !is.na(do) & run_lengths >= flatline_n
  flag_do     <- flag_do | fl
  counts[2]   <- sum(fl, na.rm = TRUE)

  # --- 3. Hampel spike filter ---
  hampel_flag <- vapply(seq_len(n), function(i) {
    if (is.na(do[i])) return(FALSE)
    idx   <- max(1L, i - hampel_half_w):min(n, i + hampel_half_w)
    x_win <- do[idx][!is.na(do[idx])]
    if (length(x_win) < 3L) return(FALSE)
    med <- median(x_win)
    s   <- mad(x_win, constant = 1)
    if (s == 0) return(FALSE)
    abs(do[i] - med) > hampel_k * s
  }, logical(1))
  flag_do   <- flag_do | hampel_flag
  counts[3] <- sum(hampel_flag, na.rm = TRUE)

  # --- 4. Step-change filter ---
  steps     <- c(NA_real_, abs(diff(do)))
  step_flag <- !is.na(steps) & !is.na(do) & steps > step_threshold
  flag_do   <- flag_do | step_flag
  counts[4] <- sum(step_flag, na.rm = TRUE)

  # --- Report ------------------------------------------------------------------
  labels <- c("physical bounds", "flatline", "Hampel spike", "step-change")
  for (i in seq_along(counts)) {
    if (counts[i] > 0)
      message(sprintf("  QC [%s]: %d obs flagged", labels[i], counts[i]))
  }
  total_new <- sum(flag_do & is.na(dat$dissolved_oxygen) == FALSE, na.rm = TRUE)
  if (sum(counts) == 0) {
    message("  QC: no flags")
  } else {
    message(sprintf("  QC total: %d obs set to NA (%.1f%%)",
                    sum(flag_do), 100 * sum(flag_do) / n))
  }

  # --- Apply flags -------------------------------------------------------------
  dat %>% mutate(
    dissolved_oxygen            = if_else(flag_do, NA_real_, dissolved_oxygen),
    dissolved_oxygen_saturation = if_else(flag_do, NA_real_, dissolved_oxygen_saturation),
    temperature                 = if_else(b_temp,  NA_real_, temperature)
  )
}

# ---- Cleaning function ------------------------------------------------------
clean_minidot <- function(file_path, start_doy = NULL) {
  dat <- read_csv(
    file_path,
    skip = 6,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    clean_names() %>%
    slice(-1) %>%           # drop units row
    mutate(
      mountain_standard_time      = ymd_hms(mountain_standard_time, tz = "America/Denver"),
      temperature                 = as.numeric(temperature),
      dissolved_oxygen            = as.numeric(dissolved_oxygen),
      dissolved_oxygen_saturation = as.numeric(dissolved_oxygen_saturation),
      q                           = as.numeric(q)
    ) %>%
    # Preserve raw sensor values before any QC flagging
    mutate(
      do_raw  = dissolved_oxygen,
      sat_raw = dissolved_oxygen_saturation,
      temp_raw = temperature
    ) %>%
    # Flag low-quality O2 readings — q reflects the optical O2 sensor only;
    # temperature (thermistor) is unaffected and retained.
    mutate(
      dissolved_oxygen            = if_else(q < q_threshold, NA_real_, dissolved_oxygen),
      dissolved_oxygen_saturation = if_else(q < q_threshold, NA_real_, dissolved_oxygen_saturation)
    ) %>%
    mutate(
      DOY         = yday(mountain_standard_time),
      Month       = month(mountain_standard_time),
      Year        = year(mountain_standard_time),
      Day         = day(mountain_standard_time),
      Hour        = hour(mountain_standard_time),
      Minute      = minute(mountain_standard_time),
      DayFrac     = DOY + (Hour * 60 + Minute) / (24 * 60),
      DecimalTime = Hour + Minute / 60
    )

  # Determine which 28 DOYs to keep
  unique_doys <- sort(unique(dat$DOY))
  if (!is.null(start_doy)) {
    keep_doys <- seq(start_doy, start_doy + 27)
  } else {
    keep_doys <- head(unique_doys, 28)
  }

  dat %>%
    filter(DOY %in% keep_doys) %>%
    mutate(
      TimeSinceStart = as.numeric(
        difftime(mountain_standard_time, min(mountain_standard_time), units = "hours")
      ),
      PAR = 0
    )
}

# ---- Loop over files --------------------------------------------------------
for (f in txt_files) {
  stream <- sub("\\.Cat\\.TXT$", "", basename(f), ignore.case = TRUE)
  message("\nCleaning: ", stream)

  cleaned <- tryCatch(
    clean_minidot(f, start_doy = doy_start_overrides[[stream]]),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  if (is.null(cleaned)) next

  # Apply DOY exclusions
  if (!is.null(doy_exclude[[stream]])) {
    n_before <- nrow(cleaned)
    cleaned  <- cleaned %>% filter(!DOY %in% doy_exclude[[stream]])
    message("  DOY exclusion: removed ", n_before - nrow(cleaned),
            " rows (DOY ", paste(doy_exclude[[stream]], collapse = ", "), ")")
  }

  # Apply point exclusions
  if (!is.null(point_exclude[[stream]])) {
    n_before <- nrow(cleaned)
    cleaned  <- cleaned %>%
      anti_join(point_exclude[[stream]], by = c("DOY", "Hour", "Minute"))
    message("  Point exclusion: removed ", n_before - nrow(cleaned), " row(s)")
  }

  # Apply automated QC (flags set to NA, rows retained for record)
  cleaned <- apply_auto_qc(cleaned, stream)

  # Recompute TimeSinceStart after all exclusions
  cleaned <- cleaned %>%
    mutate(TimeSinceStart = as.numeric(
      difftime(mountain_standard_time, min(mountain_standard_time), units = "hours")
    ))

  out_path <- file.path(out_dir, paste0(stream, "_clean.csv"))
  write_csv(cleaned, out_path)

  n_valid <- sum(!is.na(cleaned$dissolved_oxygen))
  message(sprintf("  n_total=%d  n_valid_DO=%d  DOY %d-%d  saved: %s",
                  nrow(cleaned), n_valid,
                  min(cleaned$DOY), max(cleaned$DOY), out_path))
}

message("\nDone. Cleaned files in: ", out_dir)
