# Palaeoclimate proxy data

The bin-level models in `../../R/` relate radiocarbon-inferred demographic
intensity to two independent palaeoclimate archives. 

**These series are not original to COL14.** The files in this folder are
**reduced derived extracts**  selected columns, reformatted and lightly
cleaned so the analysis scripts run efficiently. They are reproduced here for reproducibility only, from the publications cited below. They remain the intellectual property of their
original authors and publishers and are **not** covered by the COL14 data
licence (CC BY 4.0). When this repository is archived (e.g. on Zenodo), the deposit-level licence (CC BY 4.0) does **not** extend to the third-party proxy files in this folder.

Anyone reusing these series must obtain the complete data from the cited
original sources and observe the terms of those sources. To reproduce or extend
this analysis, please cite the original publications (and, where applicable, the
data-archive DOIs) listed below.

---

## 1. Laguna Pallcacocha — XRF PC1 (ENSO proxy)

Used by `xrfpallcacocha.R`.

**Original publication**
> Mark, S.Z., Abbott, M.B., Rodbell, D.T., Moy, C.M., 2022. XRF analysis of
> Laguna Pallcacocha sediments yields new insights into Holocene El Niño
> development. *Earth and Planetary Science Letters* 593, 117657.
> https://doi.org/10.1016/j.epsl.2022.117657

**Data archive (source of the extract)**
> NOAA National Centers for Environmental Information / World Data Service for
> Paleoclimatology. Laguna Pallcacocha Geochemistry, Sediment Property, and
> Flood Reconstruction Data. Dataset DOI: https://doi.org/10.25921/41b1-w758
> Landing page: https://www.ncdc.noaa.gov/access/paleo-search/study/36596
> Accessed: 2026-06-30.

The NOAA archive requests that reusers cite the original publication, the NOAA
landing page, the dataset and publication DOIs, and the date accessed.

**Expected file:** `xrf_pallcacochaPC1_mark_etal_2022.csv`

**Expected columns**

| Column     | Description |
|------------|-------------|
| `year_BP`  | Age in calendar years before present (BP). |
| `enso_pca` | First principal component of the XRF series, used as the ENSO proxy. |

---

## 2. El Junco Crater Lake — botryococcene biomarkers (ENSO proxy)

Used by `eljunco.R`.

**Original publication**

> Zhang, Z., Leduc, G., Sachs, J.P., 2014. El Niño evolution during the
> Holocene revealed by a biomarker rain gauge in the Galápagos Islands.
> *Earth and Planetary Science Letters* 404, 420–434.
> https://doi.org/10.1016/j.epsl.2014.07.013

The complete data are available only from the publisher (Elsevier) as part of
the article and its supplementary material. The extracts here were taken from
that source; reusers must obtain the full data from the publisher and observe
Elsevier's terms.

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

