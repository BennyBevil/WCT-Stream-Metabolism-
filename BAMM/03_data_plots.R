# 03_data_plots.R
# Per-stream dissolved oxygen and temperature time series plots.
#
# Output: BAMM/BaMM_data_plots.pdf
#
# Each page shows one stream over a 28-day context window:
#   - Grey points  : raw sensor measurements
#   - Red points   : QC-flagged observations (Hampel, flatline, step-change)
#   - Coloured line: retained values passed to BaMM (breaks at NA)
#   - Green lines  : selected modelling window from Data/GoodDays.csv

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(lubridate)
library(janitor)

# ---- Config -----------------------------------------------------------------
proj_root      <- "~/Documents/GPP_FP"
clean_dir      <- file.path(proj_root, "Data/O2_Files/cleaned")
raw_dir        <- file.path(proj_root, "Data/O2_Files")
good_days_file <- file.path(proj_root, "Data/GoodDays.csv")
pdf_out        <- file.path(proj_root, "BAMM/BaMM_data_plots.pdf")

# Anchor DOY for Henry (matches doy_start_overrides in 01_clean_O2.R)
henry_anchor <- 227L

# ---- Load cleaned O2 data ---------------------------------------------------
clean_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)
if (length(clean_files) == 0) stop("No cleaned O2 CSVs found in ", clean_dir)

all_o2 <- bind_rows(lapply(clean_files, function(f) {
  stream <- sub("_clean\\.csv$", "", basename(f))
  read_csv(f, show_col_types = FALSE) %>% mutate(stream = stream, .before = 1)
}))

# ---- Load GoodDays ----------------------------------------------------------
good_days <- if (file.exists(good_days_file)) {
  read_csv(good_days_file, show_col_types = FALSE)
} else {
  warning("GoodDays.csv not found — selected-range lines will be omitted.")
  NULL
}

# ---- Load raw O2: 28-day context window per stream --------------------------
load_raw_28 <- function(stream_name) {
  f <- list.files(raw_dir, pattern = paste0("^", stream_name, "\\.Cat\\.TXT$"),
                  full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) return(NULL)
  dat <- read_csv(f[1], skip = 6, show_col_types = FALSE,
                  col_types = cols(.default = col_character())) %>%
    clean_names() %>% slice(-1) %>%
    mutate(
      mountain_standard_time      = ymd_hms(mountain_standard_time, tz = "America/Denver"),
      temperature                 = as.numeric(temperature),
      dissolved_oxygen            = as.numeric(dissolved_oxygen),
      dissolved_oxygen_saturation = as.numeric(dissolved_oxygen_saturation),
      DOY     = yday(mountain_standard_time),
      Hour    = hour(mountain_standard_time),
      Minute  = minute(mountain_standard_time),
      DayFrac = DOY + (Hour * 60 + Minute) / (24 * 60)
    )
  if (toupper(stream_name) == "HENRY") {
    dat %>% filter(DOY >= henry_anchor, DOY <= henry_anchor + 27L)
  } else {
    keep <- head(sort(unique(dat$DOY)), 28)
    dat %>% filter(DOY %in% keep)
  }
}

all_raw <- bind_rows(lapply(clean_files, function(f) {
  stream <- sub("_clean\\.csv$", "", basename(f))
  rd <- load_raw_28(stream)
  if (!is.null(rd)) mutate(rd, stream = stream, .before = 1) else NULL
}))

# ---- Compute selected DOY range per stream ----------------------------------
get_selected_doys <- function(stream_name) {
  if (is.null(good_days)) return(NULL)
  gd <- good_days %>% filter(toupper(Stream) == toupper(stream_name))
  if (nrow(gd) == 0) return(NULL)
  if (toupper(stream_name) == "HENRY") {
    list(start_doy = henry_anchor + gd$DayStart - 1L,
         end_doy   = henry_anchor + gd$DayEnd   - 1L)
  } else {
    udoys <- sort(unique(filter(all_raw, stream == stream_name)$DOY))
    list(start_doy = udoys[gd$DayStart],
         end_doy   = udoys[gd$DayEnd])
  }
}

# ---- Theme ------------------------------------------------------------------
theme_ts <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())

# ---- Plot -------------------------------------------------------------------
pdf(pdf_out, width = 11, height = 8.5)

for (s in sort(unique(all_o2$stream))) {
  raw_s   <- filter(all_raw, stream == s)
  clean_s <- filter(all_o2,  stream == s)

  has_raw <- all(c("do_raw", "sat_raw", "temp_raw") %in% names(clean_s))
  if (has_raw) {
    clean_s <- clean_s %>% mutate(
      do_flagged   = is.na(dissolved_oxygen) & !is.na(do_raw),
      sat_flagged  = is.na(dissolved_oxygen_saturation) & !is.na(sat_raw),
      temp_flagged = is.na(temperature) & !is.na(temp_raw)
    )
  }

  sel <- get_selected_doys(s)
  vline_layer <- if (!is.null(sel)) {
    list(
      geom_vline(xintercept = sel$start_doy, linetype = "solid",
                 linewidth = 0.6, colour = "forestgreen"),
      geom_vline(xintercept = sel$end_doy,   linetype = "solid",
                 linewidth = 0.6, colour = "forestgreen")
    )
  }

  p_do <- ggplot(mapping = aes(x = DayFrac)) +
    geom_point(data = raw_s, aes(y = dissolved_oxygen),
               size = 0.3, colour = "grey75", alpha = 0.5) +
    { if (has_raw) geom_point(data = filter(clean_s, do_flagged),
                              aes(y = do_raw), size = 1.2,
                              colour = "firebrick", alpha = 0.8) } +
    geom_line(data = clean_s, aes(y = dissolved_oxygen),
              linewidth = 0.3, colour = "steelblue4", na.rm = TRUE) +
    vline_layer +
    labs(x = NULL, y = "DO (mg/L)") +
    theme_ts

  p_sat <- ggplot(mapping = aes(x = DayFrac)) +
    geom_point(data = raw_s, aes(y = dissolved_oxygen_saturation),
               size = 0.3, colour = "grey75", alpha = 0.5) +
    { if (has_raw) geom_point(data = filter(clean_s, sat_flagged),
                              aes(y = sat_raw), size = 1.2,
                              colour = "firebrick", alpha = 0.8) } +
    geom_line(data = clean_s, aes(y = dissolved_oxygen_saturation),
              linewidth = 0.3, colour = "steelblue2", na.rm = TRUE) +
    geom_hline(yintercept = 100, linetype = "dashed",
               linewidth = 0.3, colour = "grey50") +
    vline_layer +
    labs(x = NULL, y = "DO saturation (%)") +
    theme_ts

  p_temp <- ggplot(mapping = aes(x = DayFrac)) +
    geom_point(data = raw_s, aes(y = temperature),
               size = 0.3, colour = "grey80", alpha = 0.5) +
    { if (has_raw) geom_point(data = filter(clean_s, temp_flagged),
                              aes(y = temp_raw), size = 1.2,
                              colour = "firebrick", alpha = 0.8) } +
    geom_line(data = clean_s, aes(y = temperature),
              linewidth = 0.3, colour = "tomato3", na.rm = TRUE) +
    vline_layer +
    labs(x = "Day of year", y = "Temperature (°C)") +
    theme_ts

  sel_label <- if (!is.null(sel))
    sprintf("  |  green = selected DOY %d-%d", sel$start_doy, sel$end_doy)
  else ""

  p_combined <- (p_do / p_sat / p_temp) +
    plot_annotation(
      title = paste0("Stream: ", s,
                     "  |  grey = raw  |  red = QC flagged  |  line = retained",
                     sel_label),
      theme = theme(plot.title = element_text(size = 8, face = "bold",
                                              hjust = 0.5))
    )
  print(p_combined)
}

dev.off()
cat("Data plots written to:", pdf_out, "\n")
