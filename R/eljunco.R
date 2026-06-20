#
# BIN-LEVEL MODELs — El Junco proxies (log-botryococcene & dD)
# Each proxy modeled separately. Bins with no proxy data are dropped.
# Single linear beta (no pos/neg split).
#

library(here)  # locate repo root automatically (see R/README.md)
input_csv    <- here::here("data", "COL14_v1.0.0.csv")
logbotry_csv <- here::here("data", "paleoclimate", "eljunco_log_botryococcene.csv")
dD_csv       <- here::here("data", "paleoclimate", "eljunco_dD_botryococcene_avg.csv")
local_dir <- here::here("outputs", "eljunco")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
setwd(local_dir)
cat("Working directory:", getwd(), "\n")

log_file <- file.path(local_dir, "mcmc_binlevel_eljunco_spd.log")

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

# ANALYSIS WINDOW
WIN_OLD <- 3808
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

bins <- binPrep(sites = x$SiteCode, ages = x$C14Age, h = 100)

medDates <- medCal(calDates)
medDates[medDates > WIN_OLD] <- WIN_OLD
medDates[medDates < WIN_YNG] <- WIN_YNG

log_msg(sprintf("END: calibration | %d dates after filtering | %d unique bins (events)",
                nrow(x), length(unique(bins))))

# LOAD AND BIN BOTH EL JUNCO PROXIES
log_msg("START: El Junco proxy preparation")

# Log-botryococcene
logbotry_raw <- read.csv(logbotry_csv,
                         stringsAsFactors = FALSE)
logbotry_binned <- logbotry_raw |>
  rename(yearBP = year_BP,
         log_botry = Log_botryococcene_concentration_mgg) |>
  mutate(log_botry = as.numeric(trimws(log_botry))) |>
  filter(!is.na(log_botry)) |>
  mutate(bin_100yr = floor(yearBP / 100) * 100) |>
  group_by(bin_100yr) |>
  summarise(log_botry = median(log_botry), n_obs = n(), .groups = "drop") |>
  filter(bin_100yr >= WIN_YNG & bin_100yr <= WIN_OLD) |>
  arrange(desc(bin_100yr))

# Standardize within window
logbotry_mean <- mean(logbotry_binned$log_botry)
logbotry_sd   <- sd(logbotry_binned$log_botry)
logbotry_binned$proxy_scaled <- (logbotry_binned$log_botry - logbotry_mean) / logbotry_sd

log_msg(sprintf("  Log-botryococcene: %d bins in window | mean=%.4f sd=%.4f",
                nrow(logbotry_binned), logbotry_mean, logbotry_sd))
log_msg(sprintf("  Bins with n=1 obs: %d / %d",
                sum(logbotry_binned$n_obs == 1), nrow(logbotry_binned)))

# dD of botryococcene
dD_raw <- read.csv(dD_csv,
                   stringsAsFactors = FALSE)
dD_binned <- dD_raw |>
  rename(yearBP = year_BP,
         dD_avg = dD_botryococcene_avg) |>
  mutate(dD_avg = as.numeric(trimws(dD_avg))) |>
  filter(!is.na(dD_avg)) |>
  mutate(bin_100yr = floor(yearBP / 100) * 100) |>
  group_by(bin_100yr) |>
  summarise(dD_avg = median(dD_avg), n_obs = n(), .groups = "drop") |>
  filter(bin_100yr >= WIN_YNG & bin_100yr <= WIN_OLD) |>
  arrange(desc(bin_100yr))

# Standardize within window
dD_mean <- mean(dD_binned$dD_avg)
dD_sd   <- sd(dD_binned$dD_avg)
dD_binned$proxy_scaled <- (dD_binned$dD_avg - dD_mean) / dD_sd

log_msg(sprintf("  dD-botryococcene: %d bins in window | mean=%.4f sd=%.4f",
                nrow(dD_binned), dD_mean, dD_sd))
log_msg(sprintf("  Bins with n=1 obs: %d / %d",
                sum(dD_binned$n_obs == 1), nrow(dD_binned)))

log_msg("END: El Junco proxy preparation")

# SPD-DERIVED BIN COUNTS

log_msg("START: computing SPD-derived bin counts")

# All possible 100-yr bins in window
all_bins_in_window <- seq(floor(WIN_YNG / 100) * 100,
                          floor(WIN_OLD / 100) * 100, by = 100)
bin_lower <- all_bins_in_window
bin_upper <- all_bins_in_window + 99

bin_rep     <- !duplicated(bins)
rep_indices <- which(bin_rep)
n_events    <- length(rep_indices)

log_msg(sprintf("  %d representative events (after site binning)", n_events))

spd_counts_all <- rep(0, length(all_bins_in_window))

for (idx in rep_indices) {
  cal_grid <- calDates$grids[[idx]]
  cal_bp   <- as.numeric(row.names(cal_grid))
  cal_prob <- cal_grid$PrDens

  in_window <- cal_bp >= WIN_YNG & cal_bp <= WIN_OLD
  if (sum(cal_prob[in_window]) > 0) {
    cal_prob_win <- cal_prob * in_window
    cal_prob_win <- cal_prob_win / sum(cal_prob_win)
  } else {
    med_bin <- floor(medDates[idx] / 100) * 100
    nearest <- which.min(abs(all_bins_in_window - med_bin))
    spd_counts_all[nearest] <- spd_counts_all[nearest] + 1
    next
  }

  for (b in seq_along(all_bins_in_window)) {
    in_bin <- cal_bp >= bin_lower[b] & cal_bp <= bin_upper[b]
    spd_counts_all[b] <- spd_counts_all[b] + sum(cal_prob_win[in_bin])
  }
}

# Create a full lookup: bin -> SPD count

spd_lookup <- data.frame(
  bin       = all_bins_in_window,
  spd_count = spd_counts_all
)

log_msg(sprintf("END: SPD bin counts | %.1f total events | %d bins",
                sum(spd_counts_all), length(all_bins_in_window)))


# BUILD MATCHED DATASETS. Drop bins with no proxy data

log_msg("START: matching bins")

# --- Log-botryococcene matched set ---
matched_logbotry <- merge(
  logbotry_binned[, c("bin_100yr", "proxy_scaled", "n_obs")],
  spd_lookup,
  by.x = "bin_100yr", by.y = "bin",
  all.x = TRUE
)
matched_logbotry <- matched_logbotry[order(-matched_logbotry$bin_100yr), ]

log_msg(sprintf("  Log-botryococcene: %d matched bins (dropped %d)",
                nrow(matched_logbotry),
                length(all_bins_in_window) - nrow(matched_logbotry)))

# dD matched set
matched_dD <- merge(
  dD_binned[, c("bin_100yr", "proxy_scaled", "n_obs")],
  spd_lookup,
  by.x = "bin_100yr", by.y = "bin",
  all.x = TRUE
)
matched_dD <- matched_dD[order(-matched_dD$bin_100yr), ]

log_msg(sprintf("  dD-botryococcene:  %d matched bins (dropped %d)",
                nrow(matched_dD),
                length(all_bins_in_window) - nrow(matched_dD)))

# Save matched tables
write.csv(matched_logbotry, "matched_logbotry.csv", row.names = FALSE)
write.csv(matched_dD, "matched_dD.csv", row.names = FALSE)

log_msg("END: matching bins")

# MODEL single linear beta (no pos/neg split)

log_msg("START: specify models")

m.binlevel_proxy <- nimbleCode({
  for (t in 1:N_bins) {
    log_lambda[t] <- alpha + r * (a - bin_age[t]) + beta * proxy[t]
    lambda[t] <- exp(log_lambda[t])
    n[t] ~ dnorm(mean = lambda[t], sd = sqrt(lambda[t]))
  }
  alpha ~ dnorm(0, sd = 10)
  r     ~ dexp(1 / 0.0004)
  beta  ~ dnorm(0, sd = 1)
})

log_msg("END: specify models")

# FIT ONE PROXY MODEL


fit_proxy_model <- function(proxy_name, matched_data, model_code,
                            WIN_OLD, TOTAL_ITERS, BURNIN, THIN,
                            NCHAINS, SEEDS, NCORES) {

  log_msg(sprintf("═══ FITTING: %s ═══", proxy_name))

  constants <- list(
    N_bins  = nrow(matched_data),
    a       = WIN_OLD,
    bin_age = matched_data$bin_100yr,
    proxy   = matched_data$proxy_scaled
  )

  data_list <- list(n = matched_data$spd_count)

  mean_count <- mean(matched_data$spd_count)

  # ── Pre-flight check ──
  log_msg(sprintf("  %s: pre-flight model check", proxy_name))

  inits_test <- list(
    alpha = rnorm(1, log(mean_count + 0.1), 0.5),
    r     = rexp(1, 1 / 0.0004),
    beta  = rnorm(1, 0, 0.1)
  )

  Rmodel <- nimbleModel(
    code      = model_code,
    constants = constants,
    data      = data_list,
    inits     = inits_test
  )

  lp <- Rmodel$calculate()
  log_msg(sprintf("  %s: log-prob at inits = %.2f", proxy_name, lp))
  if (!is.finite(lp)) stop(sprintf("Non-finite log-prob for %s", proxy_name))

  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, monitors = c("alpha", "r", "beta"))
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  log_msg(sprintf("  %s: pre-flight OK", proxy_name))

  # Parallel MCMC
  run_chain <- function(seed, code, constants, data_list,
                        TOTAL_ITERS, BURNIN, THIN, mean_count) {
    library(nimble)
    library(nimbleCarbon)
    library(coda)

    set.seed(seed)

    inits <- list(
      alpha = rnorm(1, log(mean_count + 0.1), 0.5),
      r     = rexp(1, 1 / 0.0004),
      beta  = rnorm(1, 0, 0.1)
    )

    Rmodel <- nimbleModel(
      code      = code,
      constants = constants,
      data      = data_list,
      inits     = inits
    )
    Cmodel <- compileNimble(Rmodel)

    conf  <- configureMCMC(Rmodel, monitors = c("alpha", "r", "beta"))
    Rmcmc <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

    Cmcmc$run(niter = TOTAL_ITERS, nburnin = BURNIN, thin = THIN)
    return(as.matrix(Cmcmc$mvSamples))
  }

  log_msg(sprintf("  %s: starting cluster (%d cores)", proxy_name, NCORES))
  cl <- makeCluster(NCORES)

  clusterExport(
    cl,
    c("model_code", "constants", "data_list",
      "TOTAL_ITERS", "BURNIN", "THIN", "run_chain", "mean_count"),
    envir = environment()
  )
  clusterEvalQ(cl, {
    library(nimble)
    library(nimbleCarbon)
    library(coda)
  })

  log_msg(sprintf("  %s: MCMC sampling", proxy_name))
  start_mcmc <- Sys.time()

  chain_output <- parLapply(
    cl, SEEDS, run_chain,
    code        = model_code,
    constants   = constants,
    data_list   = data_list,
    TOTAL_ITERS = TOTAL_ITERS,
    BURNIN      = BURNIN,
    THIN        = THIN,
    mean_count  = mean_count
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_mcmc, units = "mins"))
  log_msg(sprintf("  %s: MCMC done | %.2f min", proxy_name, elapsed))
  stopCluster(cl)

  # ── Results ──
  chain_coda  <- mcmc.list(lapply(chain_output, mcmc))
  all_samples <- do.call(rbind, chain_output)

  scalar_params <- c("alpha", "r", "beta")

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

  cat(sprintf("\n── %s: Parameter estimates ──\n", proxy_name))
  print(param_summary)

  # Interpretation
  b_med <- param_summary$median[param_summary$parameter == "beta"]
  b_lo  <- param_summary$hpdi_lower[param_summary$parameter == "beta"]
  b_hi  <- param_summary$hpdi_upper[param_summary$parameter == "beta"]

  cat(sprintf("  beta: median=%.3f → exp(%.3f)=%.2f → %.0f%% change per +1 SD proxy\n",
              b_med, b_med, exp(b_med), (exp(b_med) - 1) * 100))

  if (b_lo > 0) {
    cat("  → 95% HPDI excludes zero (positive effect)\n")
  } else if (b_hi < 0) {
    cat("  → 95% HPDI excludes zero (negative effect)\n")
  } else {
    cat("  → 95% HPDI includes zero — no clear effect\n")
  }

  # Convergence
  cat(sprintf("\n── %s: Convergence ──\n", proxy_name))
  for (p in scalar_params) {
    ess <- effectiveSize(chain_coda[, p])
    cat(sprintf("  %s: ESS=%.0f\n", p, ess))
  }
  gelman <- gelman.diag(chain_coda[, scalar_params], multivariate = FALSE)
  print(gelman)

  for (i in seq_len(nrow(param_summary))) {
    log_msg(sprintf("  %s | %s | median=%.4f mean=%.4f sd=%.4f HPDI=[%.4f, %.4f]",
                    proxy_name, param_summary$parameter[i],
                    param_summary$median[i], param_summary$mean[i],
                    param_summary$sd[i],
                    param_summary$hpdi_lower[i], param_summary$hpdi_upper[i]))
  }

  return(list(
    proxy_name    = proxy_name,
    chain_coda    = chain_coda,
    all_samples   = all_samples,
    param_summary = param_summary,
    matched_data  = matched_data,
    constants     = constants,
    data_list     = data_list
  ))
}

#  FIT BOTH MODELS

results_logbotry <- fit_proxy_model(
  proxy_name  = "log_botryococcene",
  matched_data = matched_logbotry,
  model_code  = m.binlevel_proxy,
  WIN_OLD     = WIN_OLD,
  TOTAL_ITERS = TOTAL_ITERS,
  BURNIN      = BURNIN,
  THIN        = THIN,
  NCHAINS     = NCHAINS,
  SEEDS       = SEEDS,
  NCORES      = NCORES
)

results_dD <- fit_proxy_model(
  proxy_name  = "dD_botryococcene",
  matched_data = matched_dD,
  model_code  = m.binlevel_proxy,
  WIN_OLD     = WIN_OLD,
  TOTAL_ITERS = TOTAL_ITERS,
  BURNIN      = BURNIN,
  THIN        = THIN,
  NCHAINS     = NCHAINS,
  SEEDS       = SEEDS,
  NCORES      = NCORES
)

# SAVE

log_msg("START: saving")

write.csv(results_logbotry$param_summary, "param_summary_logbotry.csv", row.names = FALSE)
write.csv(results_dD$param_summary, "param_summary_dD.csv", row.names = FALSE)

save(results_logbotry, results_dD,
     spd_lookup, all_bins_in_window,
     logbotry_binned, dD_binned,
     logbotry_mean, logbotry_sd,
     dD_mean, dD_sd,
     matched_logbotry, matched_dD,
     WIN_OLD, WIN_YNG,
     file = "col14_binlevel_eljunco_results.RData")

chain_coda_logbotry <- results_logbotry$chain_coda
chain_coda_dD      <- results_dD$chain_coda
save(chain_coda_logbotry, file = "chain_coda_logbotry.RData")
save(chain_coda_dD, file = "chain_coda_dD.RData")

writeLines(capture.output(sessionInfo()), "sessioninfo.txt")

log_msg("END: saving")

#  10. PLOTS

log_msg("START: plots")

plot_proxy_results <- function(res, proxy_label, file_prefix) {

  matched  <- res$matched_data
  samples  <- res$all_samples
  bins_vec <- matched$bin_100yr

  # Posterior histograms
  pdf(paste0(file_prefix, "_posteriors.pdf"), width = 12, height = 4)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  for (p in c("alpha", "r", "beta")) {
    hist(samples[, p], breaks = 50,
         main = paste0(proxy_label, ": ", p), xlab = p,
         col = "steelblue", border = "white", freq = FALSE)
    abline(v = 0, col = "grey40", lty = 3)
    hpdi <- HPDinterval(as.mcmc(samples[, p]), prob = 0.95)[1, ]
    abline(v = hpdi, col = "firebrick", lty = 2, lwd = 1.5)
  }
  dev.off()

  # Traceplots
  chain_list <- res$chain_coda
  pdf(paste0(file_prefix, "_traceplots.pdf"), width = 12, height = 9)
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
  colors <- c("steelblue", "coral", "forestgreen", "purple")
  for (p in c("alpha", "r", "beta")) {
    all_vals <- do.call(c, lapply(chain_list, function(ch) ch[, p]))
    n_iter   <- nrow(chain_list[[1]])
    plot(NULL, xlim = c(1, n_iter), ylim = range(all_vals),
         main = paste0(proxy_label, " — Trace: ", p),
         xlab = "Iteration", ylab = p)
    for (ch in seq_along(chain_list)) {
      lines(chain_list[[ch]][, p], col = adjustcolor(colors[ch], 0.4), lwd = 0.5)
    }
  }
  dev.off()

  # Fitted vs observed
  pdf(paste0(file_prefix, "_fitted_vs_observed.pdf"), width = 12, height = 5)
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))

  n_post   <- min(1000, nrow(samples))
  post_idx <- seq(1, nrow(samples), length.out = n_post)

  fitted_lambda <- matrix(NA, n_post, length(bins_vec))
  for (s in seq_along(post_idx)) {
    i <- post_idx[s]
    fitted_lambda[s, ] <- exp(
      samples[i, "alpha"] +
      samples[i, "r"] * (WIN_OLD - bins_vec) +
      samples[i, "beta"] * matched$proxy_scaled
    )
  }

  lambda_median <- apply(fitted_lambda, 2, median)
  lambda_lower  <- apply(fitted_lambda, 2, quantile, 0.025)
  lambda_upper  <- apply(fitted_lambda, 2, quantile, 0.975)

  plot(bins_vec, matched$spd_count, type = "h", lwd = 4, col = "steelblue",
       xlab = "Bin (cal BP)", ylab = "Event count",
       main = paste0(proxy_label, ": Observed vs fitted"),
       xlim = rev(range(bins_vec)),
       ylim = c(0, max(c(matched$spd_count, lambda_upper)) * 1.1))
  lines(bins_vec, lambda_median, col = "firebrick", lwd = 2)
  polygon(c(bins_vec, rev(bins_vec)),
          c(lambda_lower, rev(lambda_upper)),
          col = adjustcolor("firebrick", 0.2), border = NA)
  legend("topright",
         c("Observed (SPD)", "Fitted (median)", "95% CI"),
         col = c("steelblue", "firebrick", adjustcolor("firebrick", 0.2)),
         lwd = c(4, 2, 10), bty = "n")
  dev.off()

  # Data overview: SPD + proxy
  pdf(paste0(file_prefix, "_data_overview.pdf"), width = 14, height = 5)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

  barplot(matched$spd_count, names.arg = bins_vec,
          col = "steelblue",
          main = paste0(proxy_label, ": SPD event counts"),
          xlab = "Bin (cal BP)", ylab = "Count",
          las = 2, cex.names = 0.6)

  barplot(matched$proxy_scaled, names.arg = bins_vec,
          col = ifelse(matched$proxy_scaled > 0, "coral", "steelblue"),
          main = paste0(proxy_label, ": Proxy (scaled)"),
          xlab = "Bin (cal BP)", ylab = "Proxy (z-score)",
          las = 2, cex.names = 0.6)
  abline(h = 0, col = "grey30")

  dev.off()
}

plot_proxy_results(results_logbotry, "Log-botryococcene", "logbotry")
plot_proxy_results(results_dD, "dD-botryococcene", "dD")

log_msg("END: plots")
log_msg("END: analysis complete")

cat("\n")
cat("================================================================\n")
cat("Analysis complete. Key outputs:\n")
cat("  param_summary_logbotry.csv  — log-botryococcene estimates\n")
cat("  param_summary_dD.csv        — dD estimates\n")
cat("  matched_logbotry.csv        — bins used (log-botryococcene)\n")
cat("  matched_dD.csv              — bins used (dD)\n")
cat("  *_posteriors.pdf            — marginal posteriors\n")
cat("  *_fitted_vs_observed.pdf    — model fit\n")
cat("  *_traceplots.pdf            — convergence\n")
cat("  *_data_overview.pdf         — data + proxy\n")
cat("================================================================\n")
