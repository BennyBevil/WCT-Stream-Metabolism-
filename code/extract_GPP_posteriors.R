library(dplyr)

# Extract Integrated.PP posteriors from BaMM oxygen model RDS files
# Units: mg O2 m-2 d-1  (BaMM daily average; confirmed from BaMM.tpl report line)

gpp_dir <- "data/GPP_EST"
files   <- list.files(gpp_dir, pattern = "\\.RDS$", full.names = TRUE)

# Parse site name from filename (strip _IND_BaMM.RDS)
site_from_file <- function(f) sub("_IND_BaMM\\.RDS$", "", basename(f))

posteriors <- lapply(files, function(f) {
  post <- readRDS(f)
  data.frame(
    site        = site_from_file(f),
    iter        = post$iter,
    IntegratedPP = post$Integrated.PP
  )
})

gpp_post <- bind_rows(posteriors)

# Quick summary across sites
gpp_summary <- gpp_post |>
  group_by(site) |>
  summarise(
    n_draws = n(),
    mean    = mean(IntegratedPP),
    sd      = sd(IntegratedPP),
    q2.5    = quantile(IntegratedPP, 0.025),
    q50     = quantile(IntegratedPP, 0.50),
    q97.5   = quantile(IntegratedPP, 0.975),
    .groups = "drop"
  )

print(gpp_summary)

# Save outputs
saveRDS(gpp_post, "data/GPP_EST/IntegratedPP_posteriors.RDS")
write.csv(gpp_summary, "data/GPP_EST/IntegratedPP_summary.csv", row.names = FALSE)
