# `R/` — Analysis code

The three scripts reproduce the analysis in the accompanying publication. They are the **exact publication-run scripts** — the model specification, priors, random seeds (`123, 456, 789, 325`), and iteration counts are unchanged. The *only* adaptation for this repository is the file-path handling: the scripts now locate their inputs through the [`here`](https://here.r-lib.org/) package, so **no path editing is required** (see [Paths](#paths)).

## Scripts

| Script             | What it does                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `col14_mcmc_v12.R` | **Core demographic model — run this first.** Loads and calibrates the COL14 dates (IntCal20, or Marine20 for `Material_Dated == "Shell"`), applies rcarbon site-binning (`h = 100`), computes the summed probability distribution (SPD) over 12900–500 BP, and fits a Bayesian **double-exponential growth model with a change point** (NIMBLE, `dDoubleExponentialGrowth`). Estimates the growth rates before and after the change point (`r1`, `r2`) and the change point itself (`chp`) from 4 parallel MCMC chains (100k iterations). Records run provenance (package versions, input-file MD5, RNG settings). |
| `eljunco.R`        | Bayesian bin-level regression (NIMBLE) testing whether two El Junco biomarker proxies — log-botryococcene concentration and δD of botryococcene (Zhang et al. 2014) — covary with SPD-derived demographic intensity per 100-year bin, over a 3808–500 cal BP window. Runs 4 parallel MCMC chains (100k iterations) per proxy and writes parameter summaries, diagnostics, and plots.                                                                                                                                                                                                                               |
| `xrfpallcacocha.R` | The same bin-level modelling approach using the Laguna Pallcacocha XRF PC1 ENSO proxy (Mark et al. 2022), with a positive/negative split of the standardised proxy.                                                                                                                                                                                                                                                                                                                                                                                                                                                |

The 3808–500 cal BP window used in `eljunco.R` and `xrfpallcacocha.R` is defined by the interval between the upper (older) bound (95% HPDI) of the change point (`chp`) estimated in `col14_mcmc_v12.R` and the approximate date of European contact.

It is highly recommended to **run** `col14_mcmc_v12.R` **first** when reproducing the full pipeline with a new version of the dataset. This ensures that the analysis window is correctly recalculated.

If the window changes, the corresponding ranges in `eljunco.R` and `xrfpallcacocha.R` must be updated accordingly to reflect the new bounds.

## Requirements

- **R** ≥ 4.1 (`eljunco.R` and `xrfpallcacocha.R` use the native `|>` pipe).
- R packages: `here`, `rcarbon`, `nimble`, `nimbleCarbon`, `coda`,
  `truncnorm`, `dplyr`, and `parallel` (base).
- **NIMBLE requires a working C++ toolchain** to compile models (Rtools on Windows; Xcode command-line tools on macOS; `build-essential` on Linux).
- The MCMC scripts run 4 chains in parallel; lower `NCORES` if you have fewer cores.

```r
install.packages(c("here", "rcarbon", "nimble", "nimbleCarbon",
                   "coda", "truncnorm", "dplyr"))
```

## Paths

No path editing is needed. Each script calls `library(here)` and resolves its inputs relative to the repository root — for example `here::here("data", "COL14_v1.0.0.csv")`. The
[`here`](https://here.r-lib.org/) package finds the root automatically via the `.here` anchor file (and the `.git` folder once the repo is initialised), so the scripts work whether you launch R from the repository root or from a subfolder, on any operating system.

Each script writes its outputs to its own subfolder under `outputs/` (e.g. `outputs/col14_mcmc100k/`), created automatically and ignored by Git. Curated, citable copies of the expensive results live in `results/` (see [Reproducibility](#reproducibility)).

## Running

Run with the repository (or any folder inside it) as the working directory. To reproduce everything in order:

```r
source("run_all.R")   # from the repo root
```

or run the scripts individually — `col14_mcmc_v12.R` first with a new dataset version, since it estimates the change point that defines the window used by the other two:

```r
source("R/col14_mcmc_v12.R")   # core demographic model + SPD
source("R/eljunco.R")          # El Junco proxy models
source("R/xrfpallcacocha.R")   # Pallcacocha XRF PC1 model
```

From a shell you can also run `Rscript run_all.R`.

## Outputs

Running the scripts writes to `outputs/<script>/` (ignored by Git):

- `col14_mcmc_v12.R` → the SPD plot (`SPD_12900_500_BP.pdf`), posterior chains (`col14_chain_coda.RData`, `col14_mcmc_params.RData`), the exact analysed dataset that survived calibration and binning (`col14_included_data.RData`), and provenance/session logs (`provenance.txt`, `sessionInfo*.txt`).
- `eljunco.R` → `col14_binlevel_eljunco_results.RData`, per-proxy chains (`chain_coda_logbotry.RData`, `chain_coda_dD.RData`), parameter-summary CSVs,
  and diagnostic PDFs.
- `xrfpallcacocha.R` → `col14_binlevel_results.RData`, `col14_binlevel_chain_coda.RData`, `param_summary.csv`, `bin_count_table.csv`,
  and diagnostic PDFs.

## Reproducibility


- **Precomputed results.** The expensive MCMC outputs are committed under  [`results/`](../results) so you can load and inspect them without re-running anything. See [`results/README.md`](../results/README.md).
- **Pinned package versions (recommended).** Radiocarbon results depend on the versions of `rcarbon`, `nimble`, and the IntCal20/Marine20 curves. To capture  yours exactly, initialise [`renv`](https://rstudio.github.io/renv/) and commit the lockfile:

  ```r
  install.packages("renv")
  renv::init()       # create renv.lock from the packages in use
  renv::snapshot()   # update it after any change
  ```

- **Seeds.** All chains are seeded (`123, 456, 789, 325`).
- **Provenance.** `col14_mcmc_v12.R` records the input-file MD5, R/package versions, RNG kind, and calibration curves to `provenance.txt`; each script also writes `sessionInfo()`.
- **Binary integrity.** `.gitattributes` marks `*.RData` and `*.pdf` as binary.
