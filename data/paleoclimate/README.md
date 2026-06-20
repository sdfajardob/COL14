# Palaeoclimate proxy data

The bin-level models in `../../R/` relate radiocarbon-inferred demographic
intensity to two independent palaeoclimate archives. **These series are not
original to COL14** — they are reproduced here, for reproducibility only, from
the publications cited below. They remain the intellectual property of their
original authors and publishers and are **not** covered by the COL14 data
licence (CC BY 4.0). Anyone reusing them must cite the original sources and
observe the terms of the original publications.

---

## 1. Laguna Pallcacocha — XRF PC1 (ENSO proxy)

Used by `xrfpallcacocha.R`.

**Source**

> Mark, S.Z., Abbott, M.B., Rodbell, D.T., Moy, C.M., 2022. XRF analysis of
> Laguna Pallcacocha sediments yields new insights into Holocene El Niño
> development. *Earth and Planetary Science Letters* 593, 117657.
> https://doi.org/10.1016/j.epsl.2022.117657

**Expected file:** `xrf_pallcacochaPC1_mark_etal_2022.csv`

**Expected columns**

| Column     | Description |
|------------|-------------|
| `year_BP`  | Age in calendar years before present (BP). |
| `enso_pca` | First principal component of the XRF series, used as the ENSO proxy. |

---

## 2. El Junco Crater Lake — botryococcene biomarkers (ENSO proxy)

Used by `eljunco.R`.

**Source**

> Zhang, Z., Leduc, G., Sachs, J.P., 2014. El Niño evolution during the
> Holocene revealed by a biomarker rain gauge in the Galápagos Islands.
> *Earth and Planetary Science Letters* 404, 420–434.
> https://doi.org/10.1016/j.epsl.2014.07.013

**Expected files and columns**

`eljunco_log_botryococcene.csv`

| Column                                | Description |
|---------------------------------------|-------------|
| `year_BP`                             | Age in calendar years before present (BP). |
| `Log_botryococcene_concentration_mgg` | Log botryococcene concentration (mg/g). |

`eljunco_dD_botryococcene_avg.csv`

| Column                   | Description |
|--------------------------|-------------|
| `year_BP`                | Age in calendar years before present (BP). |
| `dD_botryococcene_avg`   | Average δD of botryococcene. |

