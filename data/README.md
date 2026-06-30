# `data/`

This folder holds the COL14 radiocarbon dataset and the third-party palaeoclimate series used by the analysis.

## Files

| File                        | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| `COL14_v1.0.0.csv` | **The compiled radiocarbon dataset.**                       |
| `data_dictionary.md`        | Definition of every column.                                  |
| `SUBMISSION_TEMPLATE.csv`   | Empty CSV with the correct header row, for contributors.     |
| `paleoclimate/`             | Third-party palaeoclimate proxy series (see its own README). |

## Adding the dataset

Place your cleaned dataset in this folder named with its version, e.g.:

```
data/COL14_v1.0.0.csv
```

The analysis scripts in `../R/` currently read a file called `COL14_v1.0.0.csv` from an absolute path. After adding the dataset here, update the `read.csv(...)` path at the top of each script to point to this file (see `../R/README.md`).

## Versioning the dataset

- Every released version of the dataset is committed under a version-stamped
  name (`COL14_dataset_vMAJOR.MINOR.PATCH.csv`) and recorded in
  `../CHANGELOG.md`.
- New batches of dates bump the **MINOR** version (e.g. v1.0.0 → v1.1.0);
  corrections to existing rows bump the **PATCH** version.
- Each release is tagged on GitHub and archived to Zenodo for a citable DOI
  (see the main `README.md`).
