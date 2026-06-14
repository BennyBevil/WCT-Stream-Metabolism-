# 03_diagnostics.R
# MCMC convergence and parameter value diagnostics for BaMM posteriors.
#
# Output: BAMM/BaMM_diagnostics.pdf
#
# Sections:
#   1. ESS + flags summary table
#   2. Per-stream trace & running-mean plots (obj, Pmax, Rref, k20, Integrated.PP)
#   3. Per-stream ACF plots for the same parameters
#   4. Per-stream posterior density plots
#   5. Cross-stream posterior density comparison + Integrated.PP violin

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(coda)
library(grid)
library(gridExtra)
library(readr)

# ---- Config -----------------------------------------------------------------
proj_root  <- "~/Documents/GPP_FP"
rds_dir    <- file.path(proj_root, "Data/GPP_EST")
pdf_out    <- file.path(proj_root, "BAMM/BaMM_diagnostics.pdf")

# Parameters to diagnose (must exist as columns in the RDS)
key_params <- c("obj", "Pmax", "Rref", "k20", "Integrated.PP")

# Parameter bounds (back-transformed from BaMM.cfg log bounds) — used for
# near-bound flagging. Parameters not listed are unbounded in normal space.
param_bounds <- list(
  k20          = c(exp(-7.6009), exp( 0.6931)),   # ~0.0005 – 2.0 /hr
  init_O2_conc = c(exp( 1.6094), exp( 2.7081)),   # ~5 – 15 mg/L
  sigma_O2_conc= c(exp(-9.2103), exp( 0.0000))    # ~0.0001 – 1.0
)

ess_threshold   <- 200    # minimum acceptable ESS
near_bound_pct  <- 0.05   # flag if posterior mean within 5% of bound range
beta_extreme    <- 1e6    # flag if any beta_prod draw exceeds this in abs value

# ---- Load posteriors --------------------------------------------------------
rds_files <- list.files(rds_dir, pattern = "_IND_BaMM\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in ", rds_dir)

post_list <- lapply(rds_files, function(f) {
  stream <- sub("_IND_BaMM\\.RDS$", "", basename(f))
  readRDS(f) %>% mutate(stream = stream, .before = 1)
})
names(post_list) <- sapply(post_list, function(x) x$stream[1])

all_post <- bind_rows(post_list)

streams <- names(post_list)
cat("Streams loaded:", paste(streams, collapse = ", "), "\n")
cat("Draws per stream:", nrow(post_list[[1]]), "\n\n")

# ---- 1. ESS and flags -------------------------------------------------------
compute_ess <- function(post) {
  params <- intersect(key_params, names(post))
  ess_vals <- sapply(params, function(p) {
    x <- post[[p]]
    x <- x[is.finite(x)]
    if (length(x) < 10) return(NA_real_)
    as.numeric(effectiveSize(mcmc(x)))
  })
  tibble(param = params, ESS = round(ess_vals))
}

ess_table <- bind_rows(
  lapply(streams, function(s) {
    compute_ess(post_list[[s]]) %>% mutate(stream = s)
  })
) %>%
  pivot_wider(names_from = stream, values_from = ESS) %>%
  arrange(match(param, key_params))

# Near-bound flags
bound_flags <- bind_rows(lapply(streams, function(s) {
  post <- post_list[[s]]
  flags <- lapply(names(param_bounds), function(p) {
    if (!p %in% names(post)) return(NULL)
    x <- post[[p]]
    lo <- param_bounds[[p]][1]; hi <- param_bounds[[p]][2]
    rng <- hi - lo
    mn  <- mean(x, na.rm = TRUE)
    flag <- (mn - lo) / rng < near_bound_pct | (hi - mn) / rng < near_bound_pct
    if (flag) tibble(stream = s, param = p,
                     mean = round(mn, 4),
                     lower_bound = round(lo, 4),
                     upper_bound = round(hi, 4)) else NULL
  })
  bind_rows(flags)
}))

# Extreme beta_prod flags
beta_flags <- bind_rows(lapply(streams, function(s) {
  post <- post_list[[s]]
  if (!"beta_prod" %in% names(post)) return(NULL)
  n_extreme <- sum(abs(post$beta_prod) > beta_extreme, na.rm = TRUE)
  pct <- round(100 * n_extreme / nrow(post), 1)
  if (n_extreme > 0) tibble(stream = s, n_extreme = n_extreme,
                             pct_extreme = pct) else NULL
}))

# Print summaries to console
cat("=== ESS Table ===\n"); print(ess_table, n = Inf)
cat("\nESS < ", ess_threshold, ":\n")
low_ess <- ess_table %>%
  pivot_longer(-param, names_to = "stream", values_to = "ESS") %>%
  filter(!is.na(ESS), ESS < ess_threshold)
if (nrow(low_ess) == 0) cat("  None\n") else print(low_ess)

cat("\nNear-bound parameters:\n")
if (nrow(bound_flags) == 0) cat("  None\n") else print(bound_flags)

cat("\nExtreme beta_prod draws:\n")
if (nrow(beta_flags) == 0) cat("  None\n") else print(beta_flags)

# ---- Helpers for plotting ---------------------------------------------------
running_mean <- function(x) cumsum(x) / seq_along(x)

theme_diag <- theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ---- Open PDF ---------------------------------------------------------------
pdf(pdf_out, width = 11, height = 8.5)

# --- Page 1: ESS table -------------------------------------------------------
ess_df <- as.data.frame(ess_table)
# Flag low ESS cells with asterisk
for (s in streams) {
  col <- as.character(ess_df[[s]])
  ess_num <- as.numeric(col)
  col[!is.na(ess_num) & ess_num < ess_threshold] <-
    paste0(col[!is.na(ess_num) & ess_num < ess_threshold], " *")
  ess_df[[s]] <- col
}

title_grob <- textGrob(
  "BaMM MCMC Diagnostics — Effective Sample Size",
  gp = gpar(fontsize = 13, fontface = "bold")
)
note_grob <- textGrob(
  paste0("* = ESS < ", ess_threshold, " (consider longer chain or more thinning).\n",
         "Key params: obj = neg-log-likelihood; Pmax, Rref in model units; ",
         "Eb in eV; k20 in /hr; Integrated.PP in g O2 m-2 period-1."),
  gp = gpar(fontsize = 8), just = "left", x = 0.02
)

tbl_grob <- tableGrob(ess_df, rows = NULL,
                      theme = ttheme_minimal(base_size = 8))

grid.arrange(title_grob, tbl_grob, note_grob,
             heights = c(0.08, 0.78, 0.14), ncol = 1)

# Extreme beta_prod note (if any)
if (nrow(beta_flags) > 0) {
  beta_note <- paste0(
    "WARNING — extreme beta_prod draws (|value| > ", beta_extreme, "):\n",
    paste(sprintf("  %s: %d draws (%.1f%%)", beta_flags$stream,
                  beta_flags$n_extreme, beta_flags$pct_extreme), collapse = "\n")
  )
  grid.newpage()
  grid.text(beta_note, x = 0.05, y = 0.95, just = c("left","top"),
            gp = gpar(fontsize = 10, fontface = "bold", col = "firebrick"))
}

# --- Pages 2-N: Per-stream trace + running mean ------------------------------
for (s in streams) {
  post <- post_list[[s]]
  params_present <- intersect(key_params, names(post))

  plot_list <- lapply(params_present, function(p) {
    x   <- post[[p]]
    idx <- seq_along(x)
    df  <- tibble(draw = idx, value = x,
                  running_mean = running_mean(x))

    p_trace <- ggplot(df, aes(draw, value)) +
      geom_line(linewidth = 0.2, colour = "steelblue4", alpha = 0.8) +
      labs(x = NULL, y = p, title = NULL) +
      theme_diag

    p_rm <- ggplot(df, aes(draw, running_mean)) +
      geom_line(linewidth = 0.5, colour = "darkorange3") +
      labs(x = "draw", y = paste0("running mean\n(", p, ")")) +
      theme_diag

    p_trace / p_rm
  })

  title_text <- paste0("Stream: ", s, "  |  Trace (blue) & Running Mean (orange)")
  wrap_plots(plot_list, ncol = 3) +
    plot_annotation(title = title_text,
                    theme = theme(plot.title = element_text(
                      size = 11, face = "bold", hjust = 0.5))) %>%
    print()
}

# --- Pages N+1: Per-stream ACF plots -----------------------------------------
for (s in streams) {
  post <- post_list[[s]]
  params_present <- intersect(key_params, names(post))

  acf_plots <- lapply(params_present, function(p) {
    x <- post[[p]]
    x <- x[is.finite(x)]
    ac  <- acf(x, lag.max = 50, plot = FALSE)
    df  <- tibble(lag = ac$lag[-1], acf = ac$acf[-1])
    ci  <- qnorm(0.975) / sqrt(length(x))

    ggplot(df, aes(lag, acf)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                 colour = "firebrick", linewidth = 0.4) +
      geom_segment(aes(xend = lag, yend = 0), linewidth = 0.4,
                   colour = "steelblue4") +
      labs(x = "lag", y = "ACF", title = p) +
      theme_diag
  })

  wrap_plots(acf_plots, ncol = 3) +
    plot_annotation(
      title = paste0("Stream: ", s, "  |  Autocorrelation (dashed = 95% CI)"),
      theme = theme(plot.title = element_text(size = 11, face = "bold",
                                              hjust = 0.5))
    ) %>% print()
}

# --- Per-stream posterior density pages --------------------------------------
for (s in streams) {
  post <- post_list[[s]]
  params_present <- intersect(key_params, names(post))

  dens_plots <- lapply(params_present, function(p) {
    x <- post[[p]]
    x <- x[is.finite(x)]
    df <- tibble(value = x)

    # Median and 95% CI for subtitle
    med  <- round(median(x), 4)
    lo95 <- round(quantile(x, 0.025), 4)
    hi95 <- round(quantile(x, 0.975), 4)

    ggplot(df, aes(value)) +
      geom_histogram(aes(y = after_stat(density)), bins = 60,
                     fill = "steelblue4", alpha = 0.5, colour = NA) +
      geom_density(linewidth = 0.6, colour = "steelblue4") +
      geom_vline(xintercept = med, linetype = "solid",
                 linewidth = 0.5, colour = "darkorange3") +
      geom_vline(xintercept = c(lo95, hi95), linetype = "dashed",
                 linewidth = 0.4, colour = "darkorange3") +
      labs(x = p, y = "density",
           title = p,
           subtitle = paste0("median=", med,
                             "  95% CI [", lo95, ", ", hi95, "]")) +
      theme_diag +
      theme(plot.subtitle = element_text(size = 7, colour = "grey40"))
  })

  p_combined <- wrap_plots(dens_plots, ncol = 3) +
    plot_annotation(
      title = paste0("Stream: ", s, "  |  Posterior distributions"),
      theme = theme(plot.title = element_text(size = 11, face = "bold",
                                              hjust = 0.5))
    )
  print(p_combined)
}

# --- Final pages: Cross-stream posterior densities ---------------------------
cross_params <- c("Pmax", "Rref", "Eb", "k20", "Integrated.PP", "Integrated.CR")
cross_params  <- intersect(cross_params, names(all_post))

density_plots <- lapply(cross_params, function(p) {
  df <- all_post %>%
    select(stream, value = all_of(p)) %>%
    filter(is.finite(value))

  # Trim x-axis to 1st-99th percentile to keep bulk of distributions visible
  xlims <- quantile(df$value, c(0.05, 0.95), na.rm = TRUE)

  ggplot(df, aes(value, colour = stream, fill = stream)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    coord_cartesian(xlim = xlims) +
    labs(x = p, y = "density", colour = NULL, fill = NULL) +
    theme_diag +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7))
})

wrap_plots(density_plots, ncol = 2) +
  plot_annotation(
    title = "Cross-stream posterior densities",
    theme = theme(plot.title = element_text(size = 12, face = "bold",
                                            hjust = 0.5))
  ) %>% print()

# --- Cross-stream Integrated.PP comparison -----------------------------------
if ("Integrated.PP" %in% names(all_post)) {
  pp_df <- all_post %>%
    select(stream, Integrated.PP) %>%
    filter(is.finite(Integrated.PP)) %>%
    group_by(stream) %>%
    mutate(med = median(Integrated.PP)) %>%
    ungroup() %>%
    mutate(stream = reorder(stream, med))   # order by median

  pp_ci <- pp_df %>%
    group_by(stream) %>%
    summarise(med  = median(Integrated.PP),
              lo95 = quantile(Integrated.PP, 0.025),
              hi95 = quantile(Integrated.PP, 0.975),
              .groups = "drop")

  p_pp <- ggplot(pp_df, aes(x = Integrated.PP, y = stream)) +
    geom_violin(aes(fill = stream), alpha = 0.3, colour = NA,
                scale = "width", trim = TRUE) +
    geom_pointrange(data = pp_ci,
                    aes(x = med, xmin = lo95, xmax = hi95, y = stream),
                    linewidth = 0.6, size = 0.4, colour = "grey20") +
    labs(x = expression("Integrated PP (g O"[2]*" m"^{-2}*" period"^{-1}*")"),
         y = NULL,
         title = "Cross-stream comparison: Integrated Primary Production",
         subtitle = "Violin = posterior distribution; point = median; range = 95% CI") +
    theme_diag +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 8, colour = "grey40"))

  print(p_pp)
}

# Cross-stream summary table
summary_tbl <- all_post %>%
  select(stream, all_of(cross_params)) %>%
  pivot_longer(-stream, names_to = "param", values_to = "value") %>%
  filter(is.finite(value)) %>%
  group_by(stream, param) %>%
  summarise(mean = round(mean(value), 3),
            sd   = round(sd(value), 3),
            q2.5 = round(quantile(value, 0.025), 3),
            q50  = round(median(value), 3),
            q97.5= round(quantile(value, 0.975), 3),
            .groups = "drop") %>%
  arrange(param, stream)

grid.newpage()
grid.text("Cross-stream posterior summaries",
          x = 0.5, y = 0.97, gp = gpar(fontsize = 12, fontface = "bold"))
grid.draw(tableGrob(summary_tbl, rows = NULL,
                    theme = ttheme_minimal(base_size = 7),
                    vp = viewport(y = 0.47, height = 0.9)))

dev.off()
cat("\nDiagnostics written to:", pdf_out, "\n")
