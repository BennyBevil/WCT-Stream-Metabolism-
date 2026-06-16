# extract_daily_gpp.R
# Reusable extractor for per-day GPP posteriors from a BaMM MCMC run.
#
# Why this exists:
#   BaMM emits, for EVERY MCMC draw, the full per-interval production series P(i)
#   (mg O2 m-2 hr-1) in BaMM_mcmc.dat, but read_bamm_mcmc() in 02_run_BaMM.R only
#   keeps the scalar window-aggregate `Integrated.PP` (the vector fields are
#   coerced to NA). This script parses the P series back out and integrates it
#   per calendar day to give a daily-GPP posterior (draws x days).
#
# Model facts this code relies on (BaMM.tpl):
#   - P(i) is the instantaneous production RATE * area (area = 1 here), so the
#     production AMOUNT in interval i is  P(i) * dt(i)  with dt in hours
#     (sum_vector(i) = rate * time_step; see total_primary_production()).
#   - Integrated.PP is the mean, over moving 24-h windows, of sum_i P(i)*dt(i)
#     (calculate_integrated_moving_window_average()). extract_daily_gpp()
#     reconstructs this exactly as a correctness check.
#
# IMPORTANT interpretation:
#   BaMM fits ONE metabolic parameter set for the whole window and (with no PAR
#   sensor, SMeasFlag = 0) drives production with a deterministic solar model.
#   The resulting daily-GPP trajectory therefore reflects seasonal irradiance,
#   NOT observed day-to-day metabolic variability. CIs come from posterior
#   uncertainty in the parameters, not from independent per-day fits.

suppressMessages(library(dplyr))

# ---- Parse one named vector field from a BaMM_mcmc.dat line ------------------
# Lines are tab-delimited as: iter \t obj \t (label \t value) ...  where a vector
# value occupies a single tab-field of space-separated numbers.
.bamm_field <- function(fields, label) {
  i <- which(fields == label)[1]
  if (is.na(i)) stop("label '", label, "' not found in BaMM_mcmc.dat")
  as.numeric(strsplit(trimws(fields[i + 1]), "[ ]+")[[1]])
}

# ---- Main: daily-GPP posterior from a fresh BaMM_mcmc.dat -------------------
# mcmc_file : path to BaMM_mcmc.dat (one line per saved draw)
# time_df   : data.frame aligned to P columns, with at least
#               DOY            (integer day-of-year)
#               TimeSinceStart (cumulative hours from window start)
#             i.e. exactly the `dat` rows 02_run_BaMM.R wrote to BaMM.dat.
# Returns a list:
#   draws    : matrix (n_draws x n_days) of per-calendar-day GPP (mg O2 m-2 d-1)
#   summary  : tidy per-day posterior summary (+ hours coverage, partial flag)
#   check    : data.frame(stored, recon, max_abs_diff) validating the parse
extract_daily_gpp <- function(mcmc_file, time_df) {
  stopifnot(all(c("DOY", "TimeSinceStart") %in% names(time_df)))
  n     <- nrow(time_df)
  lines <- readLines(mcmc_file)
  flds  <- lapply(lines, function(l) strsplit(l, "\t")[[1]])

  Pmat  <- t(vapply(flds, .bamm_field, numeric(n), label = "P"))
  if (ncol(Pmat) != n)
    stop("P series length (", ncol(Pmat), ") != nrow(time_df) (", n,
         "); time grid is misaligned with the MCMC output.")
  IntPP <- vapply(flds, function(f) .bamm_field(f, "Integrated PP"), numeric(1))

  # production AMOUNT per interval, per draw: sum_vector(i) = P(i) * dt(i)
  dt   <- c(0, diff(time_df$TimeSinceStart))   # interval 1 contributes 0 (as in tpl)
  amt  <- sweep(Pmat, 2, dt, `*`)              # draws x n

  # ---- correctness check: reconstruct moving-window Integrated.PP -----------
  idx_end_first_day <- which(time_df$TimeSinceStart >= 24.0)[1]
  num_windows <- n - idx_end_first_day + 1
  cs   <- t(apply(amt, 1, cumsum))
  recon <- rowMeans(
    cs[, idx_end_first_day:(idx_end_first_day + num_windows - 1), drop = FALSE] -
    cs[, 1:num_windows, drop = FALSE])
  check <- data.frame(stored = median(IntPP), recon = median(recon),
                      max_abs_diff = max(abs(IntPP - recon)))
  if (check$max_abs_diff > 1e-2)
    warning("Reconstructed Integrated.PP differs from stored by ",
            signif(check$max_abs_diff, 3), " — check time-grid alignment.")

  # ---- per-calendar-day GPP posterior --------------------------------------
  doy   <- time_df$DOY
  days  <- sort(unique(doy))
  cover <- tapply(dt, doy, sum)
  draws <- vapply(days, function(d) rowSums(amt[, doy == d, drop = FALSE]),
                  numeric(nrow(amt)))
  colnames(draws) <- days

  summary <- tibble(
    DOY      = days,
    hours    = round(as.numeric(cover[as.character(days)]), 2),
    partial  = as.numeric(cover[as.character(days)]) < 23,
    GPP_med  = apply(draws, 2, median),
    GPP_lo95 = apply(draws, 2, quantile, 0.025),
    GPP_hi95 = apply(draws, 2, quantile, 0.975),
    GPP_sd   = apply(draws, 2, sd))

  list(draws = draws, summary = summary, check = check)
}
