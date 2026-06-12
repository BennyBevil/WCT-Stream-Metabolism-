# 01_Abundance.R
# Westslope cutthroat trout — depletion abundance estimation
# Aug–Oct 2025, 8 isolated streams
#
# Two models run in parallel for comparison:
#
#   Simple:    single q per stream, hierarchical across streams
#              logit(q[s]) ~ Normal(lqm, lqsd^2)
#
#   Size-bin:  q varies continuously with fork length via 20 mm bins
#              logit(q[s,b]) = a[s] + b_FL * (logFL_mid[b] - cenFL)
#              Gives N[s,b] posteriors — needed for size-explicit production MC.
#
# Both model files are written to disk; both sets of posteriors are saved.
# 04_Production.R loads N_post_d1 (simple) or N_post_bin_d1 (size-bin).
#
# NOTE: Jake Canyon (1 pass, n=4) and Alkali (sparse, no temp logger) excluded.
# NOTE: Cow Camp pass 2/3 correction: +15 to pass 2, -15 from pass 3
#       (misrecorded fish confirmed in field notes).
#       In the size-bin model this correction is distributed proportionally
#       across FL bins (size class of misrecorded fish unknown).

library(tidyverse)
library(lubridate)
library(R2jags)

# ============================================================
# 1. Load and clean data
# ============================================================

dat <- read.csv("Data/CTT.Data.2025.F.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(
    Date     = dmy(Date),
    Tag      = as.character(Tag),
    Location = trimws(Location)
  ) %>%
  filter(!Location %in% c("Jake Canyon", "Alkali", "Pintler"))

streams <- sort(unique(dat$Location))
Nstream <- length(streams)

dep1 <- dat %>% filter(Depletion == 1)
dep2 <- dat %>% filter(Depletion == 2)

# ============================================================
# 2. Production interval (Dep1 start -> Dep2 start per stream)
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

# ============================================================
# 3. Catch matrix helpers
# ============================================================

# --- 3a. Aggregate catch matrix [Nstream × 3 passes] — simple model ---

build_catch_matrix <- function(dep_df, stream_vec) {
  passes <- dep_df %>%
    group_by(Location, Pass) %>%
    summarise(Count = sum(Count), RL = median(RL), StreamW = median(StreamW),
              .groups = "drop") %>%
    mutate(
      Count = case_when(
        Location == "Cow" & Pass == 2 ~ Count + 15,
        Location == "Cow" & Pass == 3 ~ Count - 15,
        TRUE ~ Count
      )
    )

  full_grid <- expand.grid(Location = stream_vec, Pass = 1:3,
                           stringsAsFactors = FALSE)
  passes <- full_grid %>%
    left_join(passes, by = c("Location", "Pass")) %>%
    mutate(Count   = replace_na(Count, 0),
           RL      = ave(RL,      Location, FUN = function(x) max(x, na.rm = TRUE)),
           StreamW = ave(StreamW, Location, FUN = function(x) max(x, na.rm = TRUE)))

  mat <- passes %>%
    arrange(Location, Pass) %>%
    pivot_wider(names_from = Pass, values_from = Count, id_cols = Location) %>%
    column_to_rownames("Location") %>%
    as.matrix()
  mat <- mat[stream_vec, , drop = FALSE]

  RL_vec <- passes %>% group_by(Location) %>%
    summarise(RL = max(RL, na.rm = TRUE)) %>%
    arrange(match(Location, stream_vec)) %>% pull(RL)
  SW_vec <- passes %>% group_by(Location) %>%
    summarise(SW = max(StreamW, na.rm = TRUE)) %>%
    arrange(match(Location, stream_vec)) %>% pull(SW)

  list(mat = mat, RL = RL_vec, StreamW = SW_vec, area = RL_vec * SW_vec)
}

# --- 3b. Bin-width diagnostics: dep1 total catch per bin for candidate widths ---

cat("=== Dep1 catch per FL bin — candidate bin widths ===\n")
for (bw in c(20, 25, 30, 35, 40)) {
  b_cands <- seq(40, ceiling(max(dep1$FL[!is.na(dep1$FL)]) / bw) * bw, by = bw)
  b_df <- dep1 %>%
    filter(!is.na(FL)) %>%
    mutate(bin = cut(FL, breaks = b_cands, include.lowest = TRUE, right = FALSE)) %>%
    count(bin, .drop = FALSE)
  cat(sprintf("\n%d mm bins:\n", bw))
  print(as.data.frame(b_df), row.names = FALSE)
}
cat("\n")

# --- 3c. Size-bin catch array [Nstream × Nbin × 3 passes] ---

# FL bin width — increase to pool more fish per bin (try 25, 30, 35, 40 mm)
BIN_WIDTH <- 30

FL_all <- dat$FL[!is.na(dat$FL)]
bins   <- seq(40, ceiling(max(FL_all) / BIN_WIDTH) * BIN_WIDTH, by = BIN_WIDTH)
Nbin      <- length(bins) - 1L
bin_mids  <- (bins[-length(bins)] + bins[-1]) / 2
logFL_mid <- log(bin_mids)
cenFL     <- mean(logFL_mid)

cat(BIN_WIDTH, "mm FL bins (", Nbin, "total):\n")
cat(paste0("  [", bins[-length(bins)], ", ", bins[-1], ")"), sep = "\n")
cat("\n")

# Proportional rounding to integer target (Hamilton / largest-remainder method)
prop_round <- function(counts, target) {
  n_target <- as.integer(max(target, 0L))
  if (sum(counts) == 0L) return(integer(length(counts)))
  scaled  <- counts * n_target / sum(counts)
  floored <- floor(scaled)
  deficit <- n_target - sum(floored)
  if (deficit > 0) {
    topup <- order(scaled - floored, decreasing = TRUE)[seq_len(deficit)]
    floored[topup] <- floored[topup] + 1L
  }
  as.integer(pmax(floored, 0L))
}

build_catch_array <- function(dep_df, stream_vec, bins) {
  Nb <- length(bins) - 1L

  binned <- dep_df %>%
    filter(!is.na(FL)) %>%
    mutate(bin = as.integer(cut(FL, breaks = bins,
                                include.lowest = TRUE, right = FALSE))) %>%
    filter(!is.na(bin)) %>%
    group_by(Location, bin, Pass) %>%
    summarise(Count = sum(Count), .groups = "drop")

  full_grid <- expand.grid(Location = stream_vec, bin = seq_len(Nb), Pass = 1:3,
                           stringsAsFactors = FALSE)

  arr_df <- full_grid %>%
    left_join(binned, by = c("Location", "bin", "Pass")) %>%
    mutate(Count = replace_na(Count, 0L))

  # Cow Camp correction distributed proportionally across FL bins
  for (pass_num in c(2L, 3L)) {
    delta <- if (pass_num == 2L) 15L else -15L
    mask  <- arr_df$Location == "Cow" & arr_df$Pass == pass_num
    arr_df$Count[mask] <- prop_round(arr_df$Count[mask],
                                      sum(arr_df$Count[mask]) + delta)
  }

  arr <- array(0L, dim = c(length(stream_vec), Nb, 3L),
               dimnames = list(stream_vec,
                               paste0(bin_mids, "mm"),
                               paste0("P", 1:3)))
  for (row in seq_len(nrow(arr_df))) {
    s <- match(arr_df$Location[row], stream_vec)
    arr[s, arr_df$bin[row], arr_df$Pass[row]] <- arr_df$Count[row]
  }

  area_df <- dep_df %>%
    group_by(Location) %>%
    summarise(RL = median(RL, na.rm = TRUE),
              SW = median(StreamW, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(match(Location, stream_vec))

  list(arr = arr, area = area_df$RL * area_df$SW)
}

# Build both formats for each occasion
d1_mat  <- build_catch_matrix(dep1, streams)
d2_mat  <- build_catch_matrix(dep2, streams)
d1_arr  <- build_catch_array(dep1, streams, bins)
d2_arr  <- build_catch_array(dep2, streams, bins)

stream_area <- d1_mat$area   # same reach both occasions

# Drop bins with zero dep1 catch across all streams — these are uninformative
# and cause flat likelihoods / poor mixing in JAGS
active_bins <- which(apply(d1_arr$arr, 2, sum) > 0)

d1_arr$arr  <- d1_arr$arr[, active_bins, , drop = FALSE]
d2_arr$arr  <- d2_arr$arr[, active_bins, , drop = FALSE]
bins        <- bins[c(active_bins, max(active_bins) + 1L)]
Nbin        <- length(bins) - 1L
bin_mids    <- bin_mids[active_bins]
logFL_mid   <- logFL_mid[active_bins]
cenFL       <- mean(logFL_mid)

minN_d1 <- apply(d1_arr$arr, c(1, 2), sum)   # per stream × bin
minN_d2 <- apply(d2_arr$arr, c(1, 2), sum)

cat("--- Dep1 total catch per stream ---\n")
print(data.frame(
  Stream = streams,
  P1 = rowSums(d1_arr$arr[, , 1]),
  P2 = rowSums(d1_arr$arr[, , 2]),
  P3 = rowSums(d1_arr$arr[, , 3])
))

cat("\n--- Dep1 catch by FL bin (pooled across streams) ---\n")
print(data.frame(
  FL_range = paste0(bins[-length(bins)], "\u2013", bins[-1], " mm"),
  P1 = apply(d1_arr$arr[, , 1], 2, sum),
  P2 = apply(d1_arr$arr[, , 2], 2, sum),
  P3 = apply(d1_arr$arr[, , 3], 2, sum)
))

# ============================================================
# 4. Simple depletion model (single q per stream)
# ============================================================

cat("
model {
  lqm  ~ dunif(-3, 3)
  lqsd ~ dnorm(0, 3) T(0,)
  for (s in 1:Nstream) {
    lq[s] ~ dnorm(lqm, pow(lqsd, -2))
    q[s]  <- ilogit(lq[s])
    Np[s] ~ dpois(N[s])
    N[s]  ~ dunif(minN[s], 5000)
    catch[s,1] ~ dbin(q[s], Np[s])
    catch[s,2] ~ dbin(q[s], Np[s] - catch[s,1])
    catch[s,3] ~ dbin(q[s], Np[s] - catch[s,1] - catch[s,2])
  }
}
", file = "Depletion3pass.txt")

run_depletion_simple <- function(catch_mat, stream_vec, label) {
  cat(sprintf("\nRunning simple depletion model — %s...\n", label))
  fit <- jags(
    data               = list(catch   = catch_mat,
                               Nstream = length(stream_vec),
                               minN    = rowSums(catch_mat)),
    inits              = NULL,
    parameters.to.save = c("N", "q", "lqm", "lqsd"),
    model.file         = "Depletion3pass.txt",
    n.chains           = 3,
    n.thin             = 10,
    n.iter             = 100000,
    n.burnin           = 50000,
    DIC                = FALSE
  )
  N_post <- fit$BUGSoutput$sims.list$N
  q_post <- fit$BUGSoutput$sims.list$q
  colnames(N_post) <- stream_vec
  colnames(q_post) <- stream_vec

  cat(sprintf("\n--- Abundance posteriors (%s) ---\n", label))
  N_sum <- apply(N_post, 2, function(x)
    c(median = median(x), lo95 = quantile(x, 0.025), hi95 = quantile(x, 0.975)))
  print(round(t(N_sum)))

  cat(sprintf("\n--- Detection probability posteriors (%s) ---\n", label))
  q_sum <- apply(q_post, 2, function(x)
    c(median = median(x), lo95 = quantile(x, 0.025), hi95 = quantile(x, 0.975)))
  print(round(t(q_sum), 3))

  list(fit = fit, N_post = N_post, q_post = q_post)
}

res_simple_d1 <- run_depletion_simple(d1_mat$mat, streams, "Depletion 1 (August)")
res_simple_d2 <- run_depletion_simple(d2_mat$mat, streams, "Depletion 2 (October)")

# ============================================================
# 5. Size-bin depletion model (q ~ FL regression)
# ============================================================
# logit(q[s,b]) = a[s] + b_FL * (logFL_mid[b] - cenFL)
#
# b_FL: shared FL slope (positive => larger fish detected more easily)
# a[s]: stream-specific intercept, anchored by the simple model posterior.
#       a_prior[s] = logit(median q from simple model) — this fixes the
#       baseline detection per stream; the size-bin model only estimates
#       how q varies with FL around that known baseline.
# sigma_a: allowed deviation of a[s] from the simple-model anchor.

cat("
model {
  b_FL    ~ dnorm(0, 1)
  sigma_a ~ dnorm(0, pow(0.1, -2)) T(0,)   # half-normal; deviation from simple-model anchor

  for (s in 1:Nstream) {
    a[s] ~ dnorm(a_prior[s], pow(sigma_a, -2))

    for (b in 1:Nbin) {
      # cenFL centres the regression so a[s] = logit(q) at mean log(FL),
      # not at FL=1mm; this is purely computational (improves MCMC mixing)
      # and does not change the biological interpretation.
      # q_min provides a floor: even the smallest fish are detected at
      # some non-negligible rate, preventing q -> 0 and unbounded N.
      lq[s,b] <- a[s] + b_FL * (logFL_mid[b] - cenFL)
      q[s,b]  <- q_min + (1 - q_min) * ilogit(lq[s,b])

      N[s,b]  ~ dunif(minN[s,b], maxN)
      Np[s,b] ~ dpois(N[s,b])

      catch[s,b,1] ~ dbin(q[s,b], Np[s,b])
      catch[s,b,2] ~ dbin(q[s,b], Np[s,b] - catch[s,b,1])
      catch[s,b,3] ~ dbin(q[s,b], Np[s,b] - catch[s,b,1] - catch[s,b,2])
    }
  }
}
", file = "DepletionSizeBin.txt")

Q_MIN <- 0.05  # minimum per-pass detection probability (biologically defensible floor)

# Stream-specific prior centres for a[s]: logit of simple-model posterior median q
a_prior_d1 <- qlogis(apply(res_simple_d1$q_post, 2, median))
a_prior_d2 <- qlogis(apply(res_simple_d2$q_post, 2, median))

cat("\n--- Simple-model q (median) used as a_prior for size-bin model ---\n")
print(data.frame(Stream  = streams,
                 q_simple = round(plogis(a_prior_d1), 3),
                 a_prior  = round(a_prior_d1, 3)))

run_depletion_bins <- function(catch_arr, minN_mat, a_prior, label) {
  cat(sprintf("\nRunning size-bin depletion model — %s...\n", label))
  fit <- jags(
    data = list(
      catch     = catch_arr,
      Nstream   = dim(catch_arr)[1],
      Nbin      = dim(catch_arr)[2],
      minN      = minN_mat,
      maxN      = 500,
      logFL_mid = logFL_mid,
      cenFL     = cenFL,
      q_min     = Q_MIN,
      a_prior   = a_prior
    ),
    inits              = NULL,
    parameters.to.save = c("N", "a", "b_FL", "sigma_a"),
    model.file         = "DepletionSizeBin.txt",
    n.chains           = 3,
    n.thin             = 10,
    n.iter             = 100000,
    n.burnin           = 50000,
    DIC                = FALSE
  )

  N_post <- fit$BUGSoutput$sims.list$N   # [n_iter, Nstream, Nbin]
  dimnames(N_post)[[2]] <- streams
  dimnames(N_post)[[3]] <- paste0(bin_mids, "mm")

  N_total <- apply(N_post, c(1, 2), sum)   # [n_iter, Nstream]

  cat(sprintf("\n--- Total abundance posteriors, size-bin model (%s) ---\n", label))
  N_sum <- apply(N_total, 2, function(x)
    c(median = median(x), lo95 = quantile(x, 0.025), hi95 = quantile(x, 0.975)))
  print(round(t(N_sum)))

  b_FL_post <- fit$BUGSoutput$sims.list$b_FL
  cat(sprintf("\n--- FL detection slope b_FL (%s): median %.3f, 95%% CRI [%.3f, %.3f]\n",
              label, median(b_FL_post),
              quantile(b_FL_post, 0.025), quantile(b_FL_post, 0.975)))

  rhat_max <- max(fit$BUGSoutput$summary[, "Rhat"], na.rm = TRUE)
  n_high   <- sum(fit$BUGSoutput$summary[, "Rhat"] > 1.1, na.rm = TRUE)
  cat(sprintf("  Max Rhat = %.3f  (%d parameters > 1.1)\n", rhat_max, n_high))

  list(fit = fit, N_post = N_post, N_total = N_total, b_FL_post = b_FL_post)
}

res_bin_d1 <- run_depletion_bins(d1_arr$arr, minN_d1, a_prior_d1, "Depletion 1 (August)")
res_bin_d2 <- run_depletion_bins(d2_arr$arr, minN_d2, a_prior_d2, "Depletion 2 (October)")

# ============================================================
# 6. Comparison: simple vs size-bin total N
# ============================================================

cat("\n=== Abundance comparison: simple vs size-bin model (Dep1, median [95% CRI]) ===\n")
fmt <- function(mat) apply(mat, 2, function(x)
  sprintf("%d [%d\u2013%d]", round(median(x)),
          round(quantile(x, 0.025)), round(quantile(x, 0.975))))

print(data.frame(
  Stream    = streams,
  Simple    = fmt(res_simple_d1$N_post),
  `Size-bin` = fmt(res_bin_d1$N_total),
  check.names = FALSE
))

# Posterior median q by FL bin for each stream (size-bin model, Dep1)
b_FL_med <- median(res_bin_d1$b_FL_post)
a_med    <- apply(res_bin_d1$fit$BUGSoutput$sims.list$a, 2, median)
q_by_bin <- outer(a_med, logFL_mid - cenFL,
                  function(a, l) Q_MIN + (1 - Q_MIN) * plogis(a + b_FL_med * l))
rownames(q_by_bin) <- streams
colnames(q_by_bin) <- paste0(bin_mids, "mm")

cat("\n--- Posterior median q by FL bin and stream (size-bin model, Dep1) ---\n")
print(round(q_by_bin, 3))

# --- sigma_a posterior diagnostic ---
sigma_a_post <- res_bin_d1$fit$BUGSoutput$sims.list$sigma_a
sigma_a_q    <- quantile(sigma_a_post, c(0.025, 0.5, 0.975))
cat(sprintf("\n--- sigma_a posterior (Dep1): median %.3f, 95%% CRI [%.3f, %.3f]\n",
            sigma_a_q[2], sigma_a_q[1], sigma_a_q[3]))

hist(sigma_a_post, breaks = 60,
     main = "sigma_a posterior — size-bin model (Dep1)",
     xlab = "sigma_a (logit scale)",
     col  = "steelblue", border = "white")
abline(v = sigma_a_q, col = c("navy","red","navy"), lty = c(2,1,2), lwd = 2)

# --- Detection probability curves: overall + per stream ---

sims       <- res_bin_d1$fit$BUGSoutput$sims.list
a_sims     <- sims$a                  # [n_iter, Nstream]
b_FL_sims  <- sims$b_FL               # [n_iter]

FL_seq     <- seq(min(bin_mids), max(bin_mids), length.out = 200)
lFL_seq    <- log(FL_seq) - cenFL

# Overall curve: use mean of stream a_prior as the intercept
q_overall <- sapply(lFL_seq, function(l)
  Q_MIN + (1 - Q_MIN) * plogis(mean(a_prior_d1) + b_FL_sims * l))
# q_overall is [n_iter, 200]; get median + 95% CRI
q_ov_med  <- apply(q_overall, 2, median)
q_ov_lo   <- apply(q_overall, 2, quantile, 0.025)
q_ov_hi   <- apply(q_overall, 2, quantile, 0.975)

overall_df <- data.frame(FL = FL_seq, q = q_ov_med,
                         lo = q_ov_lo, hi = q_ov_hi,
                         Stream = "Overall (mu_a)")

# Per-stream curves: posterior median a[s]
stream_df <- do.call(rbind, lapply(seq_along(streams), function(s) {
  q_s <- sapply(lFL_seq, function(l)
    Q_MIN + (1 - Q_MIN) * plogis(a_sims[, s] + b_FL_sims * l))
  data.frame(
    FL     = FL_seq,
    q      = apply(q_s, 2, median),
    lo     = apply(q_s, 2, quantile, 0.025),
    hi     = apply(q_s, 2, quantile, 0.975),
    Stream = streams[s]
  )
}))

plot_df <- rbind(stream_df, overall_df)
plot_df$type <- ifelse(plot_df$Stream == "Overall (mu_a)", "Overall", "Stream")

p_detect <- ggplot(plot_df, aes(x = FL, y = q, colour = Stream,
                                linetype = type, fill = Stream)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = data.frame(FL = bin_mids), aes(x = FL),
           inherit.aes = FALSE, sides = "b", colour = "grey40") +
  scale_linetype_manual(values = c(Overall = "dashed", Stream = "solid"),
                        guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Fork length (mm)", y = "Detection probability (q)",
       title = "Per-pass detection probability by fork length",
       subtitle = "Size-bin model, Dep1 — median + 95% CRI",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p_detect)
ggsave("detection_curves_dep1.png", p_detect, width = 8, height = 5, dpi = 150)

# ============================================================
# 7. Save
# ============================================================

saveRDS(res_simple_d1$fit, "depletion_dep1_fit.RDS")
saveRDS(res_simple_d2$fit, "depletion_dep2_fit.RDS")

saveRDS(list(
  # --- Simple model posteriors ---
  N_post_d1     = res_simple_d1$N_post,   # [n_iter, Nstream]
  N_post_d2     = res_simple_d2$N_post,
  q_post_d1     = res_simple_d1$q_post,
  q_post_d2     = res_simple_d2$q_post,

  # --- Size-bin model posteriors ---
  N_post_bin_d1 = res_bin_d1$N_post,      # [n_iter, Nstream, Nbin]
  N_post_bin_d2 = res_bin_d2$N_post,
  N_total_bin_d1 = res_bin_d1$N_total,    # [n_iter, Nstream]
  N_total_bin_d2 = res_bin_d2$N_total,
  b_FL_post     = res_bin_d1$b_FL_post,
  a_post_d1     = res_bin_d1$fit$BUGSoutput$sims.list$a,

  # --- Bin definitions (needed by 04_Production.R) ---
  bins      = bins,
  Nbin      = Nbin,
  bin_mids  = bin_mids,
  logFL_mid = logFL_mid,
  cenFL     = cenFL,

  # --- Study metadata ---
  streams       = streams,
  Nstream       = Nstream,
  stream_area   = stream_area,
  interval_days = interval_days,

  # --- Raw data (dep1 FL distributions needed in production MC) ---
  dep1 = dep1,
  dep2 = dep2
), "abundance_results.RDS")

cat("\nDone. Results saved to abundance_results.RDS\n")
