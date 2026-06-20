# Contributing to COL14

Thank you for helping to grow and improve COL14. New radiocarbon dates, corrections to existing entries, and added references are all welcome.

All contributions are reviewed against the [data dictionary](data/data_dictionary.md) before being merged, please read it before contributing.

## Three ways to contribute

### 1. Web form — easiest, no Git required

1. Go to the repository's **Issues** tab → **New issue**.
2. Choose **"Submit new radiocarbon date(s)"** (or **"Report an error / correction"**).
3. Fill in the fields and submit.

This only requires a free GitHub account. The maintainer transfers the
validated entry into the database for you.

### 2. Spreadsheet — no GitHub account required

1. Download [`data/SUBMISSION_TEMPLATE.csv`](data/SUBMISSION_TEMPLATE.csv).
2. Fill in one row per date, following the [data dictionary](data/data_dictionary.md). Keep the header row unchanged.
3. Email the file to the maintainer (see contact in the main `README.md`).

The maintainer reviews the rows and commits them, crediting you in the commit and the changelog.

### 3. Pull request — for Git users

1. Fork the repository and create a branch.
2. Add your rows to the current `data/COL14_database_v*.csv`, following the
   data dictionary exactly (column names, controlled values, units).
3. Add an entry under **[Unreleased]** in `CHANGELOG.md`.
4. Open a pull request describing what you added and the source(s).

## What makes a good entry

- **One row per radiocarbon determination.**
- Provide at least `Lab_code`, ` C14Age`, `C14SD`, and an `APA_reference`. Geographic Coordinates are strongly preferred.
- Ages must be **uncalibrated** radiocarbon years BP — do not enter calibrated ages.
- Use the controlled values defined in the data dictionary for
  `Location_quality` and the `Shell` flag in `Material_Dated`.
- Always include the source where the date is published (`APA_reference`).
- Save spreadsheets as UTF-8, with a period as the decimal separator.

## What happens after you submit (maintainer workflow)

1. **Review** — the entry is checked for format, required fields, plausible
   ranges, and duplicates against existing `Lab_code`s.
2. **Merge** — accepted rows are added to the database CSV and the change is
   noted under **[Unreleased]** in `CHANGELOG.md`.
3. **Release** — periodically, the maintainer publishes a new version: bump the
   version in the CSV filename, `CITATION.cff`, and `CHANGELOG.md`; tag a GitHub
   release; Zenodo then mints a DOI for that version.
## Reporting errors

Found a wrong coordinate, a mis-keyed age, a missing or incorrect reference, or
a recalibration issue? Use the **"Report an error / correction"** issue form, or
email the maintainer. Please point to the specific row (e.g. by `Lab_code`) and
the corrected information with its source.
