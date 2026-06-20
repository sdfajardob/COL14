
library(here)  # locate repo root automatically (see R/README.md)
input_csv <- here::here("data", "COL14_v1.0.0.csv")
local_dir <- here::here("outputs", "col14_mcmc100k")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
setwd(local_dir)
cat("Working directory:", getwd(), "\n")

log_file <- file.path(local_dir, "mcmc_timing.log")

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg),
      file = log_file, append = TRUE)
  flush.console()
}

library(rcarbon)
library(nimble)
library(nimbleCarbon)
library(parallel)
library(coda)
library(truncnorm)

writeLines(capture.output(sessionInfo()),
           file.path(local_dir, "sessionInfo.txt"))

provenance <- c(
  sprintf("run_timestamp      : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("R_version          : %s", R.version.string),
  sprintf("rcarbon_version    : %s", as.character(packageVersion("rcarbon"))),
  sprintf("nimble_version     : %s", as.character(packageVersion("nimble"))),
  sprintf("nimbleCarbon_ver   : %s", as.character(packageVersion("nimbleCarbon"))),
  sprintf("coda_version       : %s", as.character(packageVersion("coda"))),
  sprintf("truncnorm_version  : %s", as.character(packageVersion("truncnorm"))),
  sprintf("RNGkind            : %s", paste(RNGkind(), collapse = " / ")),
  "calibration_curves : intcal20 (terrestrial), marine20 (Shell)"
)
writeLines(provenance, file.path(local_dir, "provenance.txt"))
cat(provenance, sep = "\n"); cat("\n")

# MCMC CONFIGURATION

TOTAL_ITERS <- 100000
BURNIN <- 10000
THIN <- 2
NCHAINS <- 4
NCORES <- 4
SEEDS <- c(123, 456, 789, 325)
stopifnot(length(SEEDS) == NCHAINS)
NCORES <- min(NCHAINS, parallel::detectCores())

# LOAD AND PREPARE DATA

log_msg("START: data loading")

# hash the exact input file so a reader can confirm identical files
# input_csv is defined at the top of the script (portable path via 'here')
input_md5 <- tools::md5sum(input_csv)
log_msg(sprintf("input file: %s | md5: %s", input_csv, input_md5))
cat(sprintf("input md5: %s\n", input_md5))

x <- read.csv(input_csv,
              stringsAsFactors = FALSE)

x$Material_Dated[is.na(x$Material_Dated)] <- "intcal20"

x <- subset(x,
  !is.na(C14Age) &
  !is.na(C14SD) &
  C14Age >= 500 &
  C14Age <= 12900
)

# Treat rows with missing Site_name as distinct unknown-origin sites
na_site <- is.na(x$Site_name) | x$Site_name == ""
if (any(na_site)) {
  x$Site_name[na_site] <- sprintf("UNK_%04d", seq_len(sum(na_site)))
  message(sprintf("Assigned unique placeholders to %d rows with missing Site_name",
                  sum(na_site)))
}

site_levels <- sort(unique(x$Site_name))
x$SiteCode <- factor(x$Site_name, levels = site_levels)
x$SiteCode <- paste0("S", as.integer(x$SiteCode))
x$SiteCode <- factor(x$SiteCode)

x$IsMarine <- x$Material_Dated == "Shell"
curves <- ifelse(x$IsMarine, "marine20", "intcal20")


log_msg(sprintf("END: data loading | N=%d", nrow(x)))
log_msg("START: calibration")

calDates <- calibrate(
  x = x$C14Age,
  errors = x$C14SD,
  calCurves = curves
)


index <- which.CalDates(calDates,
  BP < 12900 & BP > 499, 
  p = 0.5
)


calDates <- calDates[index]
x <- x[index, ]

bins <- binPrep(
  sites = x$SiteCode, 
  ages = x$C14Age, 
  h = 100
)

spdc = spd(
  x = calDates, 
	timeRange = c(12900,500), 
  bins = bins
)

pdf("SPD_12900_500_BP.pdf", width = 7, height = 5)
plot(spdc)
dev.off()

medDates <- medCal(calDates)
medDates[medDates > 12900] <- 12900
medDates[medDates < 500] <- 500

obs.data <- data.frame(
  LabCode = x$Lab_code,
  CRA = x$C14Age,
  Error = x$C14SD,
  MedCalDate = medDates,
  SiteID = x$SiteCode,
  IsMarine = x$IsMarine
)

#save the exact set of dates (analysed) that survived calibration + the which.CalDates(p = 0.5), bins, and the SPD object. 

save(x, index, bins, obs.data, spdc, medDates, input_md5,
     file = file.path(local_dir, "col14_included_data.RData"))

data(intcal20)

constants <- list(
  N = nrow(obs.data),
  calBP = intcal20$CalBP,
  C14BP = intcal20$C14Age,
  C14err = intcal20$C14Age.sigma
)

data_list <- list(
  X = obs.data$CRA,
  sigma = obs.data$Error
)

log_msg("END: calibration")


# Specify model


log_msg("START: specify model")

m.dblexp <- nimbleCode({
  for (i in 1:N) {
    theta[i] ~ dDoubleExponentialGrowth(
      a = 12900,
      b = 500,
      r1 = r1,
      r2 = r2,
      mu = chp
    )

    mu[i] <- interpLin(z = theta[i], x = calBP[], y = C14BP[])
    sigmaCurve[i] <- interpLin(z = theta[i], x = calBP[], y = C14err[])
    sd[i] <- sqrt(sigma[i]^2 + sigmaCurve[i]^2)

    X[i] ~ dnorm(mean = mu[i], sd = sd[i])
  }

  r1 ~ dexp(1 / 0.0004)
  r2 ~ dexp(1 / 0.0004)
  chp ~ T(dnorm(3500, sd = 200), 501, 12899)
})

initsFunction <- function(medDates) list(
  r1 = rexp(1, 1 / 0.0004),
  r2 = rexp(1, 1 / 0.0004),
  chp = truncnorm::rtruncnorm(
    1, mean = 3500, sd = 200, a = 501, b = 12899
  ),
  theta = medDates
)

log_msg("END: specify model")


# Build and compile model
log_msg("START: build model")
cat("Building and compiling model...\n")

Rmodel <- nimbleModel(
  code = m.dblexp,
  constants = constants,
  data = data_list,
  inits = initsFunction(medDates)
)

Cmodel <- compileNimble(Rmodel)

log_msg("END: model build")

log_msg("START: MCMC compile")


conf <- configureMCMC(
  Rmodel,
  monitors = c("r1", "r2", "chp", "theta")
)

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

cat("Compilation finished.\n")
log_msg("END: MCMC compile")

#parallelizing

log_msg("START: Cluster setup")

run_chain <- function(seed, code, constants, data_list, medDates) {
  library(nimble)
  library(nimbleCarbon)
  library(coda)

  set.seed(seed)

  Rmodel <- nimbleModel(
    code = code,
    constants = constants,
    data = data_list,
    inits = initsFunction(medDates)
  )

  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(
    Rmodel,
    monitors = c("r1", "r2", "chp")
  )

  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  Cmcmc$run(
    niter = TOTAL_ITERS,
    nburnin = BURNIN,
    thin = THIN
  )

  return(as.matrix(Cmcmc$mvSamples))
}

cl <- makeCluster(NCORES)

clusterExport(
  cl,
  c(
    "m.dblexp",
    "constants",
    "data_list",
    "TOTAL_ITERS",
    "BURNIN",
    "THIN",
    "initsFunction",
    "run_chain"
  ),
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

on.exit({
  log_msg(sprintf(
    "MCMC sampling aborted | elapsed %.2f min",
    as.numeric(difftime(Sys.time(), start_mcmc, units = "mins"))
  ))
}, add = TRUE)

chain_output <- parLapply(
  cl,
  SEEDS,
  run_chain,
  code = m.dblexp,
  constants = constants,
  data_list = data_list,
  medDates = medDates
)

log_msg(sprintf(
  "MCMC sampling end | elapsed %.2f min",
  as.numeric(difftime(Sys.time(), start_mcmc, units = "mins"))
))

stopCluster(cl)

library(coda)

chain_coda <- mcmc.list(
  lapply(chain_output, mcmc)
)
# save output

log_msg("START: saving data")
save(chain_coda, file = "col14_chain_coda.RData")

params <- list(
  r1 = unlist(lapply(chain_coda, function(x) as.matrix(x)[, "r1"])),
  r2 = unlist(lapply(chain_coda, function(x) as.matrix(x)[, "r2"])),
  chp = unlist(lapply(chain_coda, function(x) as.matrix(x)[, "chp"]))
)

save(params, file = "col14_mcmc_params.RData")

writeLines(capture.output(sessionInfo()),
           file.path(local_dir, "sessionInfo_end.txt"))

log_msg("END: saving data")
cat("DONE.\n")
