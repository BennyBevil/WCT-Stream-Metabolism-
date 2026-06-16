# 05_model_comparison.R
# Compares the 1-stage and 2-stage respiration models on (a) performance and
# (b) parameter / flux estimates, using the posterior draws from 02_run_BaMM.R.
#
# Prerequisites:
#   - 02_run_BaMM.R 1stage -> Data/GPP_EST/{stream}_IND_BaMM.RDS
#   - 02_run_BaMM.R 2stage -> Data/GPP_EST/{stream}_2stage_BaMM.RDS
#
# Outputs:
#   - BAMM/BaMM_model_comparison.pdf
#   - Data/GPP_EST/model_performance.csv       (per stream x model)
#   - Data/GPP_EST/model_comparison_estimates.csv (param x stream x model summaries)
#
# Performance metrics:
#   D_bar = 2*mean(obj)      deviance (lower = better fit)
#   p_D   = var(2*obj)/2     effective number of parameters (model complexity)
#   DIC   = D_bar + p_D      (lower = better; trades fit vs complexity)
#   sigma_O2_conc (median)   O2 process/observation error SD (lower = tighter fit)
#   ess_min                  smallest effective sample size across estimated params

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(coda)
library(grid)
library(gridExtra)

# ---- Config -----------------------------------------------------------------
proj_root <- "~/Documents/GPP_FP"
out_dir   <- file.path(proj_root, "Data/GPP_EST")
pdf_out   <- file.path(proj_root, "BAMM/BaMM_model_comparison.pdf")

models <- c("1stage" = "_IND_BaMM\\.RDS$", "2stage" = "_2stage_BaMM\\.RDS$")

# Flux / parameter estimands shared by both models
estimands <- c("Integrated.PP", "Integrated.CR", "Integrated.G",
               "Pmax", "I_halfSat", "Rref", "k20", "sigma_O2_conc")
# Respiration-structure parameters estimated only in the 2-stage model
twostage_only <- c("Eb", "Ep", "beta_prod")

# Pretty labels
lab <- c(Integrated.PP = "GPP (Integrated.PP)",
         Integrated.CR = "ER (Integrated.CR)",
         Integrated.G  = "Integrated.G",
         Pmax = "Pmax", I_halfSat = "I_halfSat", Rref = "Rref",
         k20 = "k20 (/hr)", sigma_O2_conc = "sigma O2 (mg/L)",
         Eb = "Eb (eV)", Ep = "Ep (eV)", beta_prod = "beta_prod")

# ---- Load posteriors --------------------------------------------------------
load_model <- function(suffix, label) {
  files <- list.files(out_dir, pattern = suffix, full.names = TRUE)
  if (length(files) == 0) { warning("No files for ", label); return(NULL) }
  setNames(lapply(files, readRDS),
           sub(suffix, "", basename(files)))
}

post <- lapply(names(models), function(m) load_model(models[[m]], m))
names(post) <- names(models)

# ---- Summaries: median + 95% CI per estimand/stream/model ------------------
summarise_draws <- function(x) {
  x <- x[is.finite(x)]
  tibble(median = median(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975),
         sd = sd(x))
}

est_tbl <- bind_rows(lapply(names(post), function(m) {
  bind_rows(lapply(names(post[[m]]), function(s) {
    p <- post[[m]][[s]]
    params <- intersect(c(estimands, twostage_only), names(p))
    bind_rows(lapply(params, function(k)
      summarise_draws(p[[k]]) %>% mutate(param = k, stream = s, model = m, .before = 1)))
  }))
}))

# ---- Performance metrics per stream/model -----------------------------------
perf_tbl <- bind_rows(lapply(names(post), function(m) {
  bind_rows(lapply(names(post[[m]]), function(s) {
    p     <- post[[m]][[s]]
    D     <- 2 * p$obj
    D_bar <- mean(D, na.rm = TRUE)
    p_D   <- var(D, na.rm = TRUE) / 2
    # ESS over estimated estimands (those that actually vary)
    varying <- intersect(c(estimands, twostage_only), names(p))
    ess <- sapply(varying, function(k) {
      x <- p[[k]][is.finite(p[[k]])]
      if (length(x) < 10 || sd(x) == 0) return(NA_real_)
      as.numeric(effectiveSize(mcmc(x)))
    })
    tibble(stream = s, model = m,
           D_bar = round(D_bar, 1),
           p_D   = round(p_D, 1),
           DIC   = round(D_bar + p_D, 1),
           sigma_O2 = round(median(p$sigma_O2_conc, na.rm = TRUE), 4),
           ess_min  = round(min(ess, na.rm = TRUE)))
  }))
}))

# Wide performance with delta-DIC and sigma ratio
perf_wide <- perf_tbl %>%
  select(stream, model, DIC, sigma_O2) %>%
  pivot_wider(names_from = model, values_from = c(DIC, sigma_O2)) %>%
  mutate(delta_DIC = round(DIC_2stage - DIC_1stage, 1),
         better    = case_when(is.na(delta_DIC) ~ "1-stage only",
                               delta_DIC < 0 ~ "2-stage",
                               TRUE ~ "1-stage"))

# ---- Write CSVs -------------------------------------------------------------
write_csv(perf_tbl, file.path(out_dir, "model_performance.csv"))
write_csv(est_tbl %>% mutate(across(where(is.numeric), ~round(.x, 4))),
          file.path(out_dir, "model_comparison_estimates.csv"))

cat("=== Performance (wide) ===\n"); print(perf_wide, n = Inf)

# ---- Plot helpers -----------------------------------------------------------
theme_cmp <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
mcols <- c("1stage" = "steelblue4", "2stage" = "darkorange3")

pointrange_panel <- function(par_name) {
  d <- est_tbl %>% filter(param == par_name)
  if (nrow(d) == 0) return(NULL)
  d <- d %>% mutate(stream = reorder(stream, median, mean, na.rm = TRUE))
  ggplot(d, aes(x = median, y = stream, colour = model)) +
    geom_pointrange(aes(xmin = q2.5, xmax = q97.5),
                    position = position_dodge(width = 0.5),
                    size = 0.35, linewidth = 0.5) +
    scale_colour_manual(values = mcols, name = NULL) +
    labs(x = lab[[par_name]], y = NULL, title = lab[[par_name]]) +
    theme_cmp
}

# ---- Open PDF ---------------------------------------------------------------
pdf(pdf_out, width = 11, height = 8.5)

# Page 1: performance table
perf_disp <- perf_tbl %>% arrange(stream, model)
grid.arrange(
  textGrob("BaMM model comparison — performance",
           gp = gpar(fontsize = 13, fontface = "bold")),
  tableGrob(perf_disp, rows = NULL, theme = ttheme_minimal(base_size = 8)),
  textGrob(paste0(
    "D_bar = deviance (lower = better fit) | p_D = effective # params | ",
    "DIC = D_bar + p_D (lower = better)\n",
    "sigma_O2 = median O2 error SD (mg/L; lower = tighter fit) | ",
    "ess_min = smallest ESS across estimated params"),
    gp = gpar(fontsize = 8), x = 0.02, just = "left"),
  heights = c(0.08, 0.80, 0.12), ncol = 1)

# Page 2: delta-DIC + sigma summary table
grid.newpage()
grid.arrange(
  textGrob("DIC & O2-error comparison (delta_DIC = 2stage - 1stage; <0 favours 2-stage)",
           gp = gpar(fontsize = 12, fontface = "bold")),
  tableGrob(perf_wide, rows = NULL, theme = ttheme_minimal(base_size = 8)),
  heights = c(0.1, 0.9), ncol = 1)

# Page 3: GPP and ER side by side (primary scientific estimands)
p_pp <- pointrange_panel("Integrated.PP")
p_cr <- pointrange_panel("Integrated.CR")
print((p_pp | p_cr) +
  plot_annotation(title = "Flux estimates: GPP and ER (median, 95% CI)",
                  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))))

# Page 4: rate parameters
rate_panels <- Filter(Negate(is.null),
  lapply(c("k20", "Rref", "Pmax", "sigma_O2_conc"), pointrange_panel))
print(wrap_plots(rate_panels, ncol = 2) +
  plot_annotation(title = "Parameter estimates by model (median, 95% CI)",
                  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))))

# Page 5: estimate shift — does respiration model change GPP/ER?
shift <- est_tbl %>%
  filter(param %in% c("Integrated.PP", "Integrated.CR")) %>%
  select(stream, model, param, median) %>%
  pivot_wider(names_from = model, values_from = median) %>%
  filter(!is.na(`1stage`), !is.na(`2stage`))
if (nrow(shift) > 0) {
  p_shift <- ggplot(shift, aes(`1stage`, `2stage`, colour = param)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(label = stream), size = 2.5, show.legend = FALSE) +
    facet_wrap(~param, scales = "free") +
    labs(x = "1-stage posterior median", y = "2-stage posterior median",
         colour = NULL,
         title = "Estimate shift between models (points on dashed line = no change)") +
    theme_cmp + theme(legend.position = "none")
  print(p_shift)
}

# Page 6: 2-stage-only respiration parameters (identifiability check)
ts_panels <- Filter(Negate(is.null),
  lapply(twostage_only, function(par_name) {
    d <- est_tbl %>% filter(param == par_name, model == "2stage")
    if (nrow(d) == 0) return(NULL)
    d <- d %>% mutate(stream = reorder(stream, median))
    ggplot(d, aes(x = median, y = stream)) +
      geom_pointrange(aes(xmin = q2.5, xmax = q97.5),
                      size = 0.35, linewidth = 0.5, colour = "darkorange3") +
      labs(x = lab[[par_name]], y = NULL, title = lab[[par_name]]) +
      theme_cmp
  }))
if (length(ts_panels) > 0)
  print(wrap_plots(ts_panels, ncol = 3) +
    plot_annotation(
      title = "2-stage respiration parameters (wide CI = weakly identified)",
      theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))))

dev.off()
cat("\nModel comparison written to:", pdf_out, "\n")
