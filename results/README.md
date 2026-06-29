# `results/` — Precomputed analysis outputs

The results of the bayesian models in this repository are expensive to run. To let others **verify the results without re-running anything**, the key outputs are committed here. Re-running the scripts writes fresh copies to `outputs/` (which is ignored by Git); this folder holds the versions that accompany the publication.


##  layout

```
results/
├── col14_mcmc_v12/
│   ├── col14_chain_coda.RData       # posterior chains: r1, r2, chp (coda mcmc.list)
│   ├── col14_mcmc_params.RData      # pooled posterior samples (r1, r2, chp)
│   ├── col14_included_data.RData    # exact analysed dates, bins, SPD object, input MD5
│   ├── col14_ppc_intcal20.RData     # Posterior predictive check data 
│   ├── SPD_12900_500_BP.pdf         # summed probability distribution
│   └── provenance.txt               # R/package versions, RNG, input-file MD5
├── eljunco/
│   ├── col14_binlevel_eljunco_results.RData
│   ├── chain_coda_logbotry.RData
│   ├── chain_coda_dD.RData
│   ├── param_summary_logbotry.csv   # headline estimates (α, r, β + HPDIs)
│   └── param_summary_dD.csv         # headline estimates (α, r, β + HPDIs)
│   └── sessioninfo.txt              # R/package versions and OS
└── xrfpallcacocha/
    ├── col14_binlevel_results.RData
    ├── col14_binlevel_chain_coda.RData
    ├── param_summary.csv
    └── bin_count_table.csv
    └── sessioninfo.txt              # R/package versions and OS
```

Small CSV summaries are included to compare re-run against the published numbers at a glance. the `.RData` files allow full inspection of the posteriors.

## Loading the results

```r
library(here)

# Posterior of the double-exponential growth model
load(here::here("results", "col14_mcmc_v12", "col14_chain_coda.RData"))
summary(chain_coda)          # 'chain_coda' is a coda mcmc.list of r1, r2, chp

# The exact dataset, bins, and SPD that the analysis used
load(here::here("results", "col14_mcmc_v12", "col14_included_data.RData"))
plot(spdc)                   # the summed probability distribution

# A proxy model's estimates
read.csv(here::here("results", "eljunco", "param_summary_logbotry.csv"))
```

Each `.RData` restores objects under their original names (e.g. `chain_coda`,`spdc`, `results_logbotry`); load one and run `ls()` to see what it contains.

## Files note

 `*.RData` files are marked binary in `.gitattributes` so Git will not corrupt them.
