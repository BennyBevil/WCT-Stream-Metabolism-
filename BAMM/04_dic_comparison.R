# 04_dic_comparison.R
# Compares the 1-stage and 2-stage respiration models by DIC, using the
# posterior obj (negative log-posterior) draws saved by 02_run_BaMM.R.
#
# Prerequisites:
#   - 02_run_BaMM.R 1stage  -> Data/GPP_EST/{stream}_IND_BaMM.RDS
#   - 02_run_BaMM.R 2stage  -> Data/GPP_EST/{stream}_2stage_BaMM.RDS
#
# Output: Data/GPP_EST/DIC_comparison.csv
#
# DIC note:
#   D(theta) = 2 * obj  (obj = negative log-posterior ≈ negative log-likelihood
#                         when priors are uniform / only bounded)
#   p_D      = var(D) / 2  (Gelman 2004 alternative p_D)
#   DIC      = D_bar + p_D
#   Strictly, WAIC needs pointwise log-likelihoods, which ADMB MCMC output does
#   not expose. DIC is the practical substitute computable from the obj column.

library(readr)
library(dplyr)
library(tidyr)

# ---- Paths ------------------------------------------------------------------
proj_root <- "~/Documents/GPP_FP"
out_dir   <- file.path(proj_root, "Data/GPP_EST")

# ---- DIC helper -------------------------------------------------------------
compute_dic <- function(obj_draws) {
  D     <- 2 * obj_draws
  D_bar <- mean(D, na.rm = TRUE)
  p_D   <- var(D, na.rm = TRUE) / 2
  tibble(D_bar = round(D_bar, 2), p_D = round(p_D, 2), DIC = round(D_bar + p_D, 2))
}

# ---- Collect DIC for one model set ------------------------------------------
dic_for <- function(pattern, model_label) {
  files <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    warning("No files matching ", pattern, " in ", out_dir)
    return(NULL)
  }
  bind_rows(lapply(files, function(f) {
    s    <- sub(pattern, "", basename(f))
    post <- readRDS(f)
    compute_dic(post$obj) %>% mutate(stream = s, model = model_label, .before = 1)
  }))
}

dic_1stage <- dic_for("_IND_BaMM\\.RDS$",    "1stage")
dic_2stage <- dic_for("_2stage_BaMM\\.RDS$", "2stage")

dic_all <- bind_rows(dic_1stage, dic_2stage) %>% arrange(stream, model)
if (nrow(dic_all) == 0) stop("No posterior RDS files found in ", out_dir)

# ---- Wide comparison: delta-DIC (2-stage minus 1-stage) ---------------------
# delta_DIC_2vs1 > 0 -> 2-stage worse (favours 1-stage)
# delta_DIC_2vs1 < 0 -> 2-stage better
dic_wide <- dic_all %>%
  select(stream, model, DIC) %>%
  pivot_wider(names_from = model, values_from = DIC) %>%
  mutate(delta_DIC_2vs1 = if (all(c("1stage", "2stage") %in% names(.)))
    round(`2stage` - `1stage`, 2) else NA_real_)

cat("\n=== DIC comparison (delta > 0 favours 1-stage) ===\n")
print(dic_wide, n = Inf)

csv_path <- file.path(out_dir, "DIC_comparison.csv")
write_csv(dic_wide, csv_path)
cat("\nDIC comparison written to:", csv_path, "\n")
