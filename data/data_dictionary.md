# COL14 — Data Dictionary

Each radiocarbon determination occupies one row. Column names are case-sensitive and contain no spaces.

---

## Columns

### `Site_name`
Reported archaeological site name, or — for individually dated objects — the object's catalogue/accession code.

- Codes containing the string **`MO`** denote artefacts curated by the **Museo del Oro (Bogotá)** that were dated either from organic residues preserved in the interior/core of the object or from an associated organic sample from the object itself with organic composition (e.g. textiles, hair).
- The name can include a locality given in **parentheses** (e.g. a municipality or *corregimiento*) records the most commonly reported geographic location, or the place name used by the source dataset or publication to identify the site.
- Rows with a missing site name are treated by the analysis scripts as distinct, unknown-origin records (each assigned a temporary `UNK_####` placeholder at run time).

*Type:* text · *Required:* recommended (may be blank for unprovenanced objects)

### `C14Age`
Measured (uncalibrated) radiocarbon age, expressed in radiocarbon years **before present (BP)**, where "present" is defined as 1950 common era (CE)

*Type:* integer · *Required:* yes

### `C14SD`
One-sigma measurement uncertainty (± 1 standard deviation) on `C14Age`, in radiocarbon years.

*Type:* integer · *Required:* yes

### `Lab_code`
Laboratory identifier of the dated sample (e.g. `Beta-123456`, `OxA-1234`, `GrN-1234`).

- Where the **sample number is known but the laboratory is not**, only the number is recorded.
- Where the **laboratory is known but the sample number is not**, the entry uses the laboratory prefix followed by a sequential placeholder   (e.g. `Beta-unknown1`, `Beta-unknown2`).

*Type:* text · *Required:* yes

### `Lon`
Longitude in decimal degrees, **WGS84** datum. Western longitudes are negative

*Type:* numeric (decimal degrees) · *Required:* recommended

### `Lat`
Latitude in decimal degrees, **WGS84** datum. Southern latitudes are negative (Colombian sites fall roughly between −4° and +13°).

*Type:* numeric (decimal degrees) · *Required:* recommended

### `Location_quality`
Reliability of the reported coordinates. One of:

| Value       | Meaning                                                                                                                                                                                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Reported`  | The location is stated in the source reference, or is otherwise well established (e.g. Teyuna / Ciudad Perdida in the Sierra Nevada de Santa Marta).                                                                                                                                                                                            |
| `Estimated` | The location is not explicitly reported but was approximated from contextual information (maps or text) in the reference. Where that information was insufficient to place the site, the coordinate instead marks the nearest present-day nucleated urban centre indicated by the reference where the archaeological activities were conducted. |
| `Unknown`   | No provenance information is available — typically objects recovered through looting activities and later curated by the Museo del Oro or another museum.                                                                                                                                                                                       |

*Type:* categorical · *Required:* yes

### `Material_Dated`
Material from which the sample was obtained, where reported (e.g. `Charcoal`, `bone`, `textile`, `organic residue`, `seed`, `Shell`).

- The value **`Shell`** flags a marine sample: the analysis calibrates these with the **Marine20** curve, and all other materials with **IntCal20**.
- Blank/missing values are treated as terrestrial (IntCal20) by the analysis  scripts.

*Type:* text · *Required:* recommended

### `APA_reference`
Bibliographic source in which the date is published or reported, formatted in
**APA** style.

*Type:* text · *Required:* yes

### `Ref_Complete`
Indicates whether the cited reference is bibliographically complete — i.e.
whether full citation details were available for the source.

*Values:* `Yes` / `No`
· *Type:* categorical · *Required:* yes

### `Archaeological_context`
Short free-text description of the materials and archaeological context associated with the date. Populated for dates that informed the **technological introduction intervals** defined in the accompanying publication.

*Type:* text · *Required:* optional

### `Technology`
The technology — among the introduction intervals analysed in the publication — to which the date pertains.

*Type:* text / categorical · *Required:* optional

---

## Notes on use

- **Calibration is performed at analysis time, not stored.** The database holds only uncalibrated `C14Age` ± `C14SD`; the R scripts calibrate with IntCal20 (terrestrial) or Marine20 (`Material_Dated == "Shell"`).
- **Encoding:** save the file as UTF-8 to preserve accented site names and author names.
- **Decimal separator:** use a period (`.`), not a comma, for `Lon`, `Lat`, .
- **Missing values:** leave the cell empty rather than entering `NA`, `0`, or
  `-`.

