# Changelog

All notable changes to the COL14 dataset and code are documented here.
Versions follow a semantic-style scheme:

- **MAJOR** — a change in structure (columns, schema) or analysis that breaks
  compatibility with earlier versions;
- **MINOR** — new radiocarbon dates added;
- **PATCH** — corrections to existing entries, metadata, or documentation.
<!--
## [Unreleased]
- Changes staged for the next release go here. -->

## [1.0.0] — 2026-06-30
### Added
- Initial public release of the COL14 radiocarbon dataset.
- Analysis code: a double-exponential growth model with a change point and the
  SPD (`col14_mcmc_v12.R`), and bin-level models against the El Junco (Zhang et
  al. 2014) and Pallcacocha XRF PC1 (Mark et al. 2022) ENSO proxies.
- Dual-licensed: code under the MIT License; the dataset under CC BY 4.0. Third-party palaeoclimate series  extracts in `data/paleoclimate/` remain under the terms of their original publications.
