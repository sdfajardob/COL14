#
# BIN-LEVEL MODEL to test whether ENSO modulates the demographic intensity with XRF PC1 Pallcacocha (Mark et al. 2022)
#

library(here)  # locate repo root automatically (see R/README.md)
input_csv <- here::here("data", "COL14_v1.0.0.csv")
enso_csv  <- here::here("data", "paleoclimate", "xrf_pallcacochaPC1_mark_etal_2022.csv")
local_dir <- here::here("outputs", "xrfpallcacocha")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
setwd(local_dir)
cat("Working directory:", getwd(), "\n")

log_file <- file.path(local_dir, "mcmc_binlevel_enso_spd.log")

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", timestamp, msg)
  cat(line)
  cat(line, file = log_file, append = TRUE)
  flush.console()
}

library(rcarbon)
library(nimble)
library(nimbleCarbon)
library(parallel)
library(coda)
library(dplyr)

# MCMC CONFIGURATION
TOTAL_ITERS <- 100000
BURNIN      <- 10000
THIN        <- 2
NCHAINS     <- 4
SEEDS       <- c(123, 456, 789, 325)
stopifnot(length(SEEDS) == NCHAINS)
NCORES <- 4

# ANALYSIS WINDOW. Older bound BP of the changing point 95% HPDI = 3811.333 from double exponential model mcmc 100k
WIN_OLD <- 3811   
WIN_YNG <- 500

# LOAD C14 DATA
log_msg("START: data loading")

x <- read.csv(input_csv, stringsAsFactors = FALSE)
x$Material_Dated[is.na(x$Material_Dated)] <- "intcal20"
x <- subset(x, !is.na(C14Age) & !is.na(C14SD) &
              C14Age >= WIN_YNG & C14Age <= WIN_OLD)

# Treat rows with missing Site_name as distinct unknown-origin sites
na_site <- is.na(x$Site_name) | x$Site_name == ""
if (any(na_site)) {
  x$Site_name[na_site] <- sprintf("UNK_%04d", seq_len(sum(na_site)))
  message(sprintf("Assigned unique placeholders to %d rows with missing Site_name",
                  sum(na_site)))
}

site_levels <- sort(unique(x$Site_name))
x$SiteCode  <- factor(x$Site_name, levels = site_levels)
x$SiteCode  <- paste0("S", as.integer(x$SiteCode))
x$SiteCode  <- factor(x$SiteCode)
x$IsMarine  <- x$Material_Dated == "Shell"
curves      <- ifelse(x$IsMarine, "marine20", "intcal20")

log_msg(sprintf("END: data loading | N=%d dates", nrow(x)))

# CALIBRATION AND SITE BINNING
log_msg("START: calibration")

calDates <- calibrate(x = x$C14Age, errors = x$C14SD, calCurves = curves)
index    <- which.CalDates(calDates, BP < WIN_OLD & BP > WIN_YNG, p = 0.5)
calDates <- calDates[index]
x        <- x[index, ]
curves   <- curves[index]

# Site-level binning = 100 yr
bins <- binPrep(sites = x$SiteCode, ages = x$C14Age, h = 100)

medDates <- medCal(calDates)
medDates[medDates > WIN_OLD] <- WIN_OLD
medDates[medDates < WIN_YNG] <- WIN_YNG

log_msg(sprintf("END: calibration | %d dates after filtering | %d unique bins (events)",
                nrow(x), length(unique(bins))))

# LOAD ENSO DATA 
log_msg("START: ENSO data preparation")

enso_raw <- read.csv(enso_csv,
                     stringsAsFactors = FALSE)
					 
#enso binning = 100 yr bin represented by lower bound
enso_binned <- enso_raw |>
  rename(yearBP = year_BP, enso_pca = enso_pca) |>
  filter(yearBP <= WIN_OLD & yearBP >= WIN_YNG) |>
  mutate(bin_100yr = floor(yearBP / 100) * 100) |>
  group_by(bin_100yr) |>
  summarise(enso_pca = median(enso_pca), .groups = "drop") |> #I checked the distributions of median is better than the mean here.
  arrange(desc(bin_100yr))

enso_lookup <- enso_binned$enso_pca
enso_bins   <- enso_binned$bin_100yr

# Standardize ENSO
enso_mean <- mean(enso_lookup)
enso_sd   <- sd(enso_lookup)

enso_scaled <- (enso_lookup - enso_mean) / enso_sd
enso_pos    <- ifelse(enso_scaled > 0, enso_scaled, 0)
enso_neg    <- ifelse(enso_scaled < 0, abs(enso_scaled), 0)

log_msg(sprintf("END: ENSO data | %d bins | mean=%.4f sd=%.4f",
                length(enso_bins), enso_mean, enso_sd))

# SPD-DERIVED BIN COUNTS. it distributes each date calibrated probability across bins.

log_msg("START: computing SPD-derived bin counts")

bin_lower <- enso_bins
bin_upper <- enso_bins + 99

# One representative date per site-bin group
bin_rep     <- !duplicated(bins)
rep_indices <- which(bin_rep)
n_events    <- length(rep_indices)

log_msg(sprintf("  %d representative events (after site binning)", n_events))

# For each representative, extract calibrated probability and distribute across 100-yr time bins
spd_counts <- rep(0, length(enso_bins))

for (idx in rep_indices) {
  # Extract calibrated probability grid for this date
  cal_grid <- calDates$grids[[idx]]
  cal_bp   <- as.numeric(row.names(cal_grid))
  cal_prob <- cal_grid$PrDens

  # Normalize to sum to 1 (within our window)
  in_window <- cal_bp >= WIN_YNG & cal_bp <= WIN_OLD
  if (sum(cal_prob[in_window]) > 0) {
    cal_prob_win <- cal_prob * in_window
    cal_prob_win <- cal_prob_win / sum(cal_prob_win)
  } else {
    # Fallback: use median date if no probability in window
    med_bin <- floor(medDates[idx] / 100) * 100
    nearest <- which.min(abs(enso_bins - med_bin))
    spd_counts[nearest] <- spd_counts[nearest] + 1
    next
  }

  # Sum probability within each 100-yr bin
  for (b in seq_along(enso_bins)) {
    in_bin <- cal_bp >= bin_lower[b] & cal_bp <= bin_upper[b]
    spd_counts[b] <- spd_counts[b] + sum(cal_prob_win[in_bin])
  }
}

#  median-based counts for comparison
med_rep <- medDates[bin_rep]
med_bin_assign <- floor(med_rep / 100) * 100
med_bin_assign <- sapply(med_bin_assign, function(b) {
  enso_bins[which.min(abs(enso_bins - b))]
})
median_counts <- sapply(enso_bins, function(b) sum(med_bin_assign == b))

# raw counts (all dates, not site-binned) for reference
time_bin_all <- floor(medDates / 100) * 100
time_bin_all <- sapply(time_bin_all, function(b) {
  enso_bins[which.min(abs(enso_bins - b))]
})
raw_counts <- sapply(enso_bins, function(b) sum(time_bin_all == b))

# Summary
count_table <- data.frame(
  bin           = enso_bins,
  spd_count     = round(spd_counts, 3),
  median_count  = median_counts,
  spd_minus_med = round(spd_counts - median_counts, 3),
  raw_count     = raw_counts,
  enso_pca      = round(enso_lookup, 4),
  enso_scaled   = round(enso_scaled, 4),
  enso_pos      = round(enso_pos, 4),
  enso_neg      = round(enso_neg, 4)
)
write.csv(count_table, "bin_count_table.csv", row.names = FALSE)

log_msg(sprintf("END: SPD bin counts | %.1f total events | %d bins",
                sum(spd_counts), length(enso_bins)))
log_msg(sprintf("  SPD count range:    [%.1f, %.1f]",
                min(spd_counts), max(spd_counts)))
log_msg(sprintf("  Median count range: [%d, %d]",
                min(median_counts), max(median_counts)))
log_msg(sprintf("  Mean |SPD - median| per bin: %.2f",
                mean(abs(spd_counts - median_counts))))

#  MODEL It does a normal approximation to Poisson using dnorm. 
log_msg("START: specify model")


m.binlevel <- nimbleCode({
  for (t in 1:N_bins) {
    log_lambda[t] <- alpha + r * (a - bin_age[t]) +
                     beta_pos * ENSO_pos[t] + beta_neg * ENSO_neg[t]
    lambda[t] <- exp(log_lambda[t])
    n[t] ~ dnorm(mean = lambda[t], sd = sqrt(lambda[t]))
  }
  alpha    ~ dnorm(0, sd = 10)
  r        ~ dexp(1 / 0.0004)
  beta_pos ~ dnorm(0, sd = 1)
  beta_neg ~ dnorm(0, sd = 1)
})

constants <- list(
  N_bins   = length(enso_bins),
  a        = WIN_OLD,
  bin_age  = enso_bins,
  ENSO_pos = enso_pos,
  ENSO_neg = enso_neg
)

data_list <- list(
  n = spd_counts
)

initsFunction <- function() list(
  alpha    = rnorm(1, log(mean(spd_counts) + 0.1), 0.5),
  r        = rexp(1, 1 / 0.0004),
  beta_pos = rnorm(1, 0, 0.1),
  beta_neg = rnorm(1, 0, 0.1)
)

log_msg("END: specify model")

# 6. COMPILE and pre-flight check
log_msg("START: build model check")

Rmodel <- nimbleModel(
  code      = m.binlevel,
  constants = constants,
  data      = data_list,
  inits     = initsFunction()
)

lp <- Rmodel$calculate()
log_msg(sprintf("Model log-prob at inits: %.2f", lp))
if (!is.finite(lp)) stop("Non-finite log-prob at initial values.")

Cmodel <- compileNimble(Rmodel)
log_msg("END: model build check")

log_msg("START: pre-flight MCMC compile")

conf <- configureMCMC(
  Rmodel,
  monitors = c("alpha", "r", "beta_pos", "beta_neg")
)

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

log_msg("END: pre-flight MCMC compile")

# PARALLELIZING
log_msg("START: cluster setup")

run_chain <- function(seed, code, constants, data_list,
                      TOTAL_ITERS, BURNIN, THIN, mean_count) {
  library(nimble)
  library(nimbleCarbon)
  library(coda)

  set.seed(seed)

  inits <- list(
    alpha    = rnorm(1, log(mean_count + 0.1), 0.5),
    r        = rexp(1, 1 / 0.0004),
    beta_pos = rnorm(1, 0, 0.1),
    beta_neg = rnorm(1, 0, 0.1)
  )

  Rmodel <- nimbleModel(
    code      = code,
    constants = constants,
    data      = data_list,
    inits     = inits
  )

  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(
    Rmodel,
    monitors = c("alpha", "r", "beta_pos", "beta_neg")
  )

  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  Cmcmc$run(
    niter   = TOTAL_ITERS,
    nburnin = BURNIN,
    thin    = THIN
  )

  return(as.matrix(Cmcmc$mvSamples))
}

cl <- makeCluster(NCORES)

mean_count <- mean(spd_counts)

clusterExport(
  cl,
  c("m.binlevel", "constants", "data_list",
    "TOTAL_ITERS", "BURNIN", "THIN", "run_chain", "mean_count"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(nimble)
  library(nimbleCarbon)
  library(coda)
})

log_msg("END: cluster setup")
log_msg("START: MCMC sampling")

start_mcmc <- Sys.time()

chain_output <- parLapply(
  cl,
  SEEDS,
  run_chain,
  code        = m.binlevel,
  constants   = constants,
  data_list   = data_list,
  TOTAL_ITERS = TOTAL_ITERS,
  BURNIN      = BURNIN,
  THIN        = THIN,
  mean_count  = mean_count
)

elapsed <- as.numeric(difftime(Sys.time(), start_mcmc, units = "mins"))
log_msg(sprintf("END: MCMC sampling | elapsed %.2f min", elapsed))

stopCluster(cl)

# OUTPUT
log_msg("START: processing output")

chain_coda <- mcmc.list(lapply(chain_output, mcmc))

# Combined samples across chains
all_samples <- do.call(rbind, chain_output)

# Parameter summaries
scalar_params <- c("alpha", "r", "beta_pos", "beta_neg")

param_summary <- data.frame(
  parameter  = scalar_params,
  median     = sapply(scalar_params, function(p) round(median(all_samples[, p]), 4)),
  mean       = sapply(scalar_params, function(p) round(mean(all_samples[, p]), 4)),
  sd         = sapply(scalar_params, function(p) round(sd(all_samples[, p]), 4)),
  hpdi_lower = sapply(scalar_params, function(p)
    round(HPDinterval(as.mcmc(all_samples[, p]), prob = 0.95)[1, "lower"], 4)),
  hpdi_upper = sapply(scalar_params, function(p)
    round(HPDinterval(as.mcmc(all_samples[, p]), prob = 0.95)[1, "upper"], 4)),
  row.names = NULL
)

cat("\n── Parameter estimates ──────────────────────────────\n")
print(param_summary)
write.csv(param_summary, "param_summary.csv", row.names = FALSE)

for (i in seq_len(nrow(param_summary))) {
  log_msg(sprintf("Estimate | %s | median=%.4f mean=%.4f sd=%.4f HPDI=[%.4f, %.4f]",
                  param_summary$parameter[i],
                  param_summary$median[i], param_summary$mean[i],
                  param_summary$sd[i],
                  param_summary$hpdi_lower[i], param_summary$hpdi_upper[i]))
}

# Interpretation
cat("\n── Interpretation ──────────────────────────────────\n")
bp_med <- param_summary$median[param_summary$parameter == "beta_pos"]
bn_med <- param_summary$median[param_summary$parameter == "beta_neg"]
bp_lo  <- param_summary$hpdi_lower[param_summary$parameter == "beta_pos"]
bp_hi  <- param_summary$hpdi_upper[param_summary$parameter == "beta_pos"]
bn_lo  <- param_summary$hpdi_lower[param_summary$parameter == "beta_neg"]
bn_hi  <- param_summary$hpdi_upper[param_summary$parameter == "beta_neg"]

cat(sprintf("  beta_pos: median=%.3f → exp(%.3f)=%.2f → %.0f%% more events per +1 SD ENSO\n",
            bp_med, bp_med, exp(bp_med), (exp(bp_med) - 1) * 100))
cat(sprintf("  beta_neg: median=%.3f → exp(%.3f)=%.2f → %.0f%% fewer events per +1 SD |ENSO|\n",
            bn_med, bn_med, exp(bn_med), (1 - exp(bn_med)) * 100))

if (bp_lo > 0) {
  cat("  → beta_pos: 95% HPDI excludes zero — credible El Nino effect\n")
} else {
  cat("  → beta_pos: 95% HPDI includes zero — no clear El Nino effect\n")
}

if (bn_hi < 0) {
  cat("  → beta_neg: 95% HPDI excludes zero — credible La Nina effect\n")
} else {
  cat("  → beta_neg: 95% HPDI includes zero — no clear La Nina effect\n")
}

log_msg("END: processing output")

# CONVERGENCE DIAGNOSTICS
log_msg("START: convergence diagnostics")

cat("\n── Convergence diagnostics ─────────────────────────\n")

for (p in scalar_params) {
  ess_combined <- effectiveSize(chain_coda[, p])
  cat(sprintf("  %s: ESS=%.0f\n", p, ess_combined))
}

gelman <- gelman.diag(chain_coda[, scalar_params], multivariate = FALSE)
cat("\n  Gelman-Rubin R-hat:\n")
print(gelman)

log_msg("END: convergence diagnostics")

# SAVE OUTPUT
log_msg("START: saving")

save(chain_coda, file = "col14_binlevel_chain_coda.RData")

save(all_samples, param_summary,
     count_table, spd_counts, median_counts, raw_counts,
     enso_bins, enso_lookup, enso_scaled, enso_pos, enso_neg,
     enso_mean, enso_sd,
     constants, data_list,
     WIN_OLD, WIN_YNG,
     file = "col14_binlevel_results.RData")

writeLines(capture.output(sessionInfo()), "sessioninfo.txt")

log_msg("END: saving")

# PLOTS
log_msg("START: plots")

# Plot 1: Data overview (SPD vs median vs raw)
pdf("data_overview.pdf", width = 14, height = 10)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

# SPD counts per bin
barplot(spd_counts, names.arg = enso_bins,
        col = ifelse(enso_scaled > 0, "coral", "steelblue"),
        main = "A. SPD-derived event counts per bin",
        xlab = "Bin (cal BP)", ylab = "Count (fractional)",
        las = 2, cex.names = 0.6)

# ENSO values
barplot(enso_scaled, names.arg = enso_bins,
        col = ifelse(enso_scaled > 0, "coral", "steelblue"),
        main = "B. ENSO covariate (scaled)",
        xlab = "Bin (cal BP)", ylab = "ENSO (z-score)",
        las = 2, cex.names = 0.6)
abline(h = 0, col = "grey30")

# SPD vs median counts
barplot(rbind(spd_counts, median_counts), beside = TRUE,
        names.arg = enso_bins,
        col = c("steelblue", "coral"),
        main = "C. SPD-derived vs median-based counts",
        xlab = "Bin (cal BP)", ylab = "Count",
        las = 2, cex.names = 0.6)
legend("topright", c("SPD (fractional)", "Median (integer)"),
       fill = c("steelblue", "coral"), bty = "n")

# Scatter — SPD count vs ENSO
plot(enso_scaled, spd_counts, pch = 16, col = "steelblue",
     xlab = "ENSO (scaled)", ylab = "SPD event count",
     main = "D. Event count vs ENSO")
points(enso_scaled[enso_scaled > 0], spd_counts[enso_scaled > 0],
       pch = 16, col = "coral")
legend("topright", c("ENSO > 0 (El Nino)", "ENSO < 0 (La Nina)"),
       col = c("coral", "steelblue"), pch = 16, bty = "n")

dev.off()

# Plot 2: Posterior histograms
pdf("posterior_histograms.pdf", width = 16, height = 4)
par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 1))

for (p in scalar_params) {
  hist(all_samples[, p], breaks = 50,
       main = p, xlab = p,
       col = "steelblue", border = "white", freq = FALSE)
  abline(v = 0, col = "grey40", lty = 3)
  hpdi <- HPDinterval(as.mcmc(all_samples[, p]), prob = 0.95)[1, ]
  abline(v = hpdi, col = "firebrick", lty = 2, lwd = 1.5)
}
dev.off()

# Plot 3: Traceplots (all chains)
pdf("traceplots.pdf", width = 12, height = 10)
par(mfrow = c(4, 1), mar = c(4, 4, 3, 1))

colors <- c("steelblue", "coral", "forestgreen", "purple")
for (p in scalar_params) {
  plot(NULL, xlim = c(1, nrow(chain_output[[1]])),
       ylim = range(sapply(chain_output, function(ch) range(ch[, p]))),
       main = paste("Trace:", p), xlab = "Iteration", ylab = p)
  for (ch in seq_along(chain_output)) {
    lines(chain_output[[ch]][, p], col = adjustcolor(colors[ch], 0.4), lwd = 0.5)
  }
}
dev.off()

# Plot 4: Posterior correlations
pdf("posterior_correlations.pdf", width = 10, height = 10)
par(mfrow = c(4, 4), mar = c(4, 4, 2, 1))

cor_mat <- cor(all_samples[, scalar_params])
thin_idx <- seq(1, nrow(all_samples), length.out = min(2000, nrow(all_samples)))

for (i in seq_along(scalar_params)) {
  for (j in seq_along(scalar_params)) {
    if (i == j) {
      dens <- density(all_samples[, scalar_params[i]])
      plot(dens, main = scalar_params[i], xlab = "", col = "steelblue", lwd = 2)
    } else {
      plot(all_samples[thin_idx, scalar_params[j]],
           all_samples[thin_idx, scalar_params[i]],
           pch = ".", col = rgb(0.2, 0.4, 0.8, 0.3),
           xlab = scalar_params[j], ylab = scalar_params[i],
           main = sprintf("cor = %.3f", cor_mat[i, j]))
    }
  }
}
dev.off()

# Plot 5: Fitted vs observed
pdf("fitted_vs_observed.pdf", width = 12, height = 5)
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))

alpha_samples <- all_samples[, "alpha"]
r_samples     <- all_samples[, "r"]
bp_samples    <- all_samples[, "beta_pos"]
bn_samples    <- all_samples[, "beta_neg"]

n_post <- min(1000, nrow(all_samples))
post_idx <- seq(1, nrow(all_samples), length.out = n_post)

fitted_lambda <- matrix(NA, n_post, length(enso_bins))
for (s in seq_along(post_idx)) {
  idx <- post_idx[s]
  fitted_lambda[s, ] <- exp(
    alpha_samples[idx] +
    r_samples[idx] * (WIN_OLD - enso_bins) +
    bp_samples[idx] * enso_pos +
    bn_samples[idx] * enso_neg
  )
}

lambda_median <- apply(fitted_lambda, 2, median)
lambda_lower  <- apply(fitted_lambda, 2, quantile, 0.025)
lambda_upper  <- apply(fitted_lambda, 2, quantile, 0.975)

plot(enso_bins, spd_counts, type = "h", lwd = 4, col = "steelblue",
     xlab = "Bin (cal BP)", ylab = "Event count",
     main = "Observed SPD counts vs fitted intensity",
     xlim = rev(range(enso_bins)),
     ylim = c(0, max(c(spd_counts, lambda_upper)) * 1.1))
lines(enso_bins, lambda_median, col = "firebrick", lwd = 2)
polygon(c(enso_bins, rev(enso_bins)),
        c(lambda_lower, rev(lambda_upper)),
        col = adjustcolor("firebrick", 0.2), border = NA)
legend("topright",
       c("Observed (SPD)", "Fitted (median)", "95% CI"),
       col = c("steelblue", "firebrick", adjustcolor("firebrick", 0.2)),
       lwd = c(4, 2, 10), bty = "n")

dev.off()

log_msg("END: plots")
log_msg("END: analysis complete")

cat("\n")
cat("================================================================\n")
cat("Analysis complete. Key outputs:\n")
cat("  param_summary.csv        — parameter estimates + HPDI\n")
cat("  bin_count_table.csv      — SPD, median, raw counts + ENSO\n")
cat("  posterior_histograms.pdf — marginal posteriors\n")
cat("  fitted_vs_observed.pdf  — model fit to data\n")
cat("  traceplots.pdf           — chain convergence\n")
cat("  data_overview.pdf        — SPD vs median counts + ENSO\n")
cat("================================================================\n")
