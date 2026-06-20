# COL14 — A compiled radiocarbon database for Colombia

COL14 is a version-controlled database of radiocarbon dates from Colombia, together with the R code used to analyse it in the accompanying publication. The repository is designed also to include new radiocarbon dates from contributors who do not use Git.

The database is a result from the Mapping the Archaeological Pre-Columbian Heritage of South America (_MAPHSA_) project.  We used the database to identify what social en environmental conditions enabled the production of technologies during ~12900-500BP in northern South American. The database includes current available data for directly dating of the emergence of four technologies: plant domestication, pottery, metallurgy, and raised-field structures.

- **Maintainer:** Sebastian Fajardo — Leiden University; Delft University of Technology
- **Database version:** 1.0.0 · see [`CHANGELOG.md`](CHANGELOG.md)
- **DOI:** <!-- TODO: add Zenodo concept DOI badge once archived -->
- **Associated publication:** <!-- TODO: citation + DOI when available -->

---

## What's in here

```
COL14/
├── data/
│   ├── COL14_v1.0.0.csv               # the radiocarbon database
│   ├── data_dictionary.md             # definition of every column in the database
│   ├── SUBMISSION_TEMPLATE.csv        # blank template for contributing dates
│   └── paleoclimate/                  # third-party proxy series + sources
├── R/                                 # analysis scripts + how to run them
├── results/                           # precomputed MCMC outputs
├── run_all.R                          # reproduce the whole pipeline in order
├── .github/ISSUE_TEMPLATE/            # web forms for contributing/correcting dates
├── CITATION.cff                       # how to cite COL14
├── CONTRIBUTING.md                    # how to add or correct dates
├── CHANGELOG.md                       # what changed in each version
├── LICENSE                            # MIT — applies to the code
└── DATA-LICENSE.md                    # CC BY 4.0 — applies to the database
```

## The database

Each row is one radiocarbon determination. Ages are stored **uncalibrated** (`C14Age` ± `C14SD`, in radiocarbon years BP); calibration happens in the analysis code. Every column is defined in
[`data/data_dictionary.md`](data/data_dictionary.md).

## Reproducing the analysis

The analysis is in R. See [`R/README.md`](R/README.md) for the pipeline overview and package requirements. The scripts locate their inputs automatically (via the `here` package).[`run_all.R`](run_all.R) from the repository root to reproduce everything in order, or source the scripts individually. In brief, the three scripts (1) calibrate the dates, compute the SPD, and fit a Bayesian double-exponential growth model with a change point, and (2–3) fit Bayesian bin-level models relating demographic intensity to two independent ENSO/palaeoclimate proxies. It is recommend to run `col14_mcmc_v12.R` first as stand alone because  it estimates the change point that defines the window used by the other two scripts.

The MCMC is slow to run, so precomputed results are committed under [`results/`](results) — you can load and inspect the posteriors without re-running anything. For an exactly reproducible software environment, the repository is set up to use [`renv`](https://rstudio.github.io/renv/) (see the Reproducibility section of `R/README.md`).

The palaeoclimate series are reproduced from Mark et al. (2022) and Zhang et al. (2014); see [`data/paleoclimate/README.md`](data/paleoclimate/README.md) for full citations and licensing.

## Contributing new dates

New radiocarbon dates and corrections are welcome. There are three ways to contribute, from no technical setup to full Git — all described in [`CONTRIBUTING.md`](CONTRIBUTING.md):

1. **Web form (no Git needed).** Open an issue using the *Submit new radiocarbon date(s)* form and fill in the fields. Requires only a free GitHub account.
2. **Spreadsheet (no GitHub needed).** Fill in [`data/SUBMISSION_TEMPLATE.csv`](data/SUBMISSION_TEMPLATE.csv) and send it to
   the maintainer; it will be reviewed and committed for you.
3. **Pull request.** Edit the database CSV directly and open a PR.

All contributions are reviewed against the [data dictionary](data/data_dictionary.md) before being merged, please read it before contributing. 

## Versioning and permanent archiving

The database follows semantic-style versioning (`Major.Minor.Patch` format): new batches of dates bump the minor version (1.0.0 → 1.1.0); corrections bump the patch version. Each release is tagged on GitHub and recorded in [`CHANGELOG.md`](CHANGELOG.md).

<!-- TODO: This repository is connected to  **[Zenodo](https://zenodo.org/)**: DOI: -->

## How to cite

If you use COL14, please cite both the database and the associated publication. Citation metadata is in [`CITATION.cff`](CITATION.cff); GitHub renders a *"Cite this repository"* button from it. Once archived, cite the Zenodo DOI.

## Licensing

This repository is **dual-licensed**:

- **Code** (everything in `R/`) — [MIT License](LICENSE).
- **Database** (`data/COL14_v*.csv` and the data dictionary) —
  [Creative Commons Attribution 4.0 International (CC BY 4.0)](DATA-LICENSE.md).
  You may share and adapt the data, including commercially, as long as you give
  appropriate credit.

The third-party palaeoclimate series in `data/paleoclimate/` are **not** covered by these licences and remain subject to the terms of their original publications.

## Contact

Sebastian Fajardo — Leiden University / Delft University of Technology.
Website: https://sdfajardob.github.io/site/contact/  
Email: s.d.fajardo.bernal@liacs.leidenuniv.nl 
