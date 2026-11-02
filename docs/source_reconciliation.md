# Source reconciliation: GDC Release 46.0, TCGA-CDR and cBioPortal (TCGA COAD/READ)

**Status:** checkpoint closed (29 October 2026)
**Scope:** TCGA-COAD and TCGA-READ clinical data only. No genomic or RNA-seq data.
**Purpose:** before building the ETL, establish how the three available TCGA clinical sources relate to each other, where they agree, where they differ, and why.

| Role | Source |
|---|---|
| Primary ETL source | NCI GDC, Data Release 46.0 |
| Curated survival reference | TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR; Liu et al., *Cell* 2018) |
| Reconciliation only (not an ETL input) | cBioPortal, *Colorectal Adenocarcinoma (TCGA, PanCancer Atlas)* |

Evidence labels used throughout:
- **[V] Verified:** demonstrated directly by the scripts listed in the appendix.
- **[H] Hypothesis:** a plausible explanation consistent with the evidence but not proven.
- **[D] Deferred decision:** an ETL design choice informed by this work but **not yet made**.

---

## 1. Sources and provenance

| Source | Version / identifier | Acquisition | Files used |
|---|---|---|---|
| GDC | Data Release 46.0, 10 Aug 2026 (API status: commit `8f7c2a51`, tag 8.5.0) | Manual GDC portal download, 29 Oct 2026; projects TCGA-COAD + TCGA-READ | Clinical JSON; clinical TSV archive; biospecimen TSV archive |
| TCGA-CDR | Supplementary Table S1, Liu et al. 2018 | GDC PanCanAtlas publication page | `TCGA-CDR-SupplementalTableS1.xlsx`, sheet `TCGA-CDR` |
| cBioPortal | Study `coadread_tcga_pan_can_atlas_2018`; datahub commit `e1b75cd25bf7436d16719b7d76117ad8f3f66e54` (27 Oct 2026) | Scripted download pinned to commit, SHA-256 verified | `data_clinical_patient.txt`, `data_clinical_sample.txt` |

SHA-256 checksums of every input file are listed in Appendix A.2.

Other cBioPortal TCGA CRC studies considered: `coadread_tcga` (Firehose Legacy, 2016 snapshot; a possible later independent comparator) and `coadread_tcga_pub` (Nature 2012, 276 patients; too limited).

Caveat [V/H]: cBioPortal PanCanAtlas clinical data was built from the PanCanAtlas resources, so it is **not fully independent** of TCGA-CDR. Agreement between those two partly reflects shared origin.

---

## 2. Patient counts and overlap

### 2.1 Patients per source [V]

| Source | All | COAD | READ | Duplicate IDs | Missing IDs |
|---|---|---|---|---|---|
| GDC R46 | 633 | 461 | 172 | 0 | 0 |
| TCGA-CDR | 629 | 459 | 170 | 0 | 0 |
| cBioPortal (patient table) | 594 | 439 | 155 | 0 | 0 |

cBioPortal holds exactly one sample per patient (594 samples, all primary tumour, `-01`). GDC holds 2,037 samples for the same cohort (1,263 primary tumour, 657 blood normal, 112 solid-tissue normal, 3 metastatic, 2 recurrent). [V]

### 2.2 Overlap [V]

| In GDC | In CDR | In cBioPortal | Patients |
|---|---|---|---|
| ✓ | ✓ | ✓ | 591 |
| ✓ | ✓ | – | 38 |
| ✓ | – | ✓ | 3 |
| ✓ | – | – | 1 |

- No patient exists outside GDC: GDC contains every CDR and cBioPortal patient. [V]
- **The sources are not fully nested.** CDR ⊂ GDC and cBioPortal ⊂ GDC, but cBioPortal ⊄ CDR: 3 cBioPortal patients (TCGA-5M-AAT5, TCGA-5M-AATA, TCGA-F5-6810) are absent from the CDR. Source sizes (633 > 629 > 594) had suggested a simple funnel; the explicit test disproved it. [V]
- COAD/READ project labels are 100% concordant across sources (629/629 GDC–CDR; 594/594 GDC–cBioPortal). [V]

### 2.3 Patients missing from a source

- **4 GDC patients absent from the CDR** (TCGA-5M-AAT5, TCGA-5M-AATA, TCGA-AF-3912, TCGA-F5-6810): all have vital status "Not Reported" and no follow-up records in GDC, and cBioPortal has no overall-survival data for the 3 it contains. **No source holds outcome data for them.** [V] Their exclusion from a survival resource is consistent with having nothing to curate. [H]
- **39 GDC patients absent from cBioPortal:** all 39 have a non-FFPE primary-tumour sample in GDC, so "no tumour sample" is **rejected** as the explanation. [V] Exclusion from the PanCanAtlas molecular freeze for sample- or data-quality reasons (most are early TCGA tissue-source sites, e.g. AA, AG) is the remaining explanation. [H] It is testable with the PanCanAtlas merged sample-quality annotations file (not done; not needed for the ETL).
- **TCGA-A6-2670** is marked "Redacted" in the CDR, is absent from cBioPortal, and has no GDC follow-up records. [V] The reason for redaction is unknown. [H]
- **TCGA-A6-2670 and TCGA-AF-3912** are the only 2 GDC cases with no diagnosis labelled `primary` (each has a single diagnosis labelled "not reported"). [V]

### 2.4 Project label vs anatomical site [V]

Within GDC, the COAD/READ project label does not strictly correspond to anatomy (ICD-O-3 topography of the `primary` diagnosis):
- TCGA-COAD includes 7 rectosigmoid-junction tumours.
- TCGA-READ includes 71 rectosigmoid-junction, 6 "Colon, NOS", 2 "connective tissue of abdomen" and 1 "unknown primary site".

---

## 3. Clinical-variable concordance

Fields were compared only where definitions are compatible, after harmonising coding and units. Each comparison is restricted to patients present in both sources and classified as agree / disagree / missing in one source / missing in both. GDC values were taken from the diagnosis labelled `primary` **for comparison only**.

### 3.1 Agreement where both sources report a value [V]

| Variable | GDC–CDR | GDC–cBioPortal | CDR–cBioPortal |
|---|---|---|---|
| Sex | 629/629 | 591/591 | — |
| Primary site, ICD-O-3 (GDC `tissue_or_organ_of_origin`) | — | 590/590 | — |
| Histology, ICD-O-3 (GDC `morphology`) | — | 590/590 | — |
| AJCC pathological sub-stage | 608/608 | 574/574 | 574/574 |
| AJCC stage group (I–IV) | 608/608 | 574/574 | — |
| AJCC staging edition | — | 523/523 | — |
| Pathological T / N / M | — | 591 / 591 / 584 (all agree) | — |
| Vital status | 627/629 | 589/591 | 591/591 |
| Age (years) | ±1 year only (see 3.3) | ±1 year only | 591/591 |
| Overall-survival time (days) | — | — | 590/590 (after unit conversion) |

**Staging:** where two sources both report a stage, they never disagree, at sub-stage level, across all pairs.

### 3.2 Missingness [V]

| Item | Finding |
|---|---|
| Missing-value markers | Each source differs: GDC uses `Not Reported`, `unknown` and null; CDR uses `[Not Available]` and `[Discrepancy]`; cBioPortal uses empty strings. TSV exports use `'--`. |
| AJCC stage | Missing for 25 GDC patients (primary diagnosis). CDR: 20 missing (12 `[Not Available]`, 8 `[Discrepancy]`). GDC lacks a stage that the CDR has for 1 patient, and one that cBioPortal has for 6. |
| AJCC edition | Missing in both GDC and cBioPortal for 70 patients. |
| `[Discrepancy]` | A CDR curation flag meaning the original records conflicted (stage 8, histology 3). It is not simple missingness. |
| Vital status | "Not Reported" for 4 GDC patients: the same 4 with no outcome data anywhere. |

### 3.3 Definitional differences (apparent disagreements that are not contradictions)

| Topic | Finding | Evidence |
|---|---|---|
| Survival units | cBioPortal `OS_MONTHS` = days ÷ (365/12) = 30.4167, for 568/568 patients with non-zero time; **not** 365.25/12. | [V] |
| Age in years | No conversion of GDC day counts reproduces CDR age: ÷365.25 gives 21 differences (all −1), ÷365 gives 12 (all +1), GDC `age_at_index` gives 21 (all −1). CDR and cBioPortal ages match exactly. Ages in whole years are reliable only to ±1 year. | [V] |
| Why ages differ | TCGA probably computed ages from real calendar dates that are no longer available; cBioPortal ages are probably derived from the CDR. | [H] |
| Birth offsets | 4 GDC–CDR differences of exactly 1 day (−32,872 vs −32,873); all 4 are among the 19 patients with `age_is_obfuscated = true` (ages ≥ 90 capped for privacy). | [V] |
| Two GDC site fields | `tissue_or_organ_of_origin` matches cBioPortal ICD-O-3 site 590/590; `site_of_resection_or_biopsy` differs for 200 patients. They represent tumour origin vs tissue-sampling site. | [V] |
| Histology coding | CDR uses 4 broad text labels; GDC and cBioPortal use ICD-O-3 codes. The CDR folds rarer types (8263, 8211, 8255, 8560, 8574, 8010) into "Adenocarcinoma"; only mucinous yes/no is comparable. | [V] |
| Last follow-up | CDR leaves last contact blank for deceased patients by design (114 of 115 blanks are deaths); CDR `OS.time` = death day for 128/128 deaths and last contact for 500/501 living. cBioPortal `DAYS_LAST_FOLLOWUP` for deceased patients is inconsistent (equal to death day for 62; blank or different for 58). GDC has two candidate fields: diagnosis-level `days_to_last_follow_up` and the maximum `days_to_follow_up` across follow-up records. | [V] |

---

## 4. Unresolved contradictions

Patient IDs are listed in `results/exploration/05_flagged_contradictions.csv`.

| Issue | n | Evidence [V] | Explanation |
|---|---|---|---|
| Days to death, GDC vs CDR | 7 | GDC values 0, 1, 91, 570 days vs CDR 92–2,047; in every case GDC death day equals GDC's own last follow-up day; cBioPortal agrees with CDR. | [H] GDC death day populated from a follow-up record or placeholder. |
| Vital status, GDC dead vs CDR/cBioPortal alive | 2 | TCGA-F5-6702: GDC death at day 869, after CDR last contact at day 452. TCGA-CL-4957: GDC "Dead" with no death day. | [H] F5-6702 reflects a GDC update after the 2018 CDR snapshot. CL-4957 is unexplained. |
| Last follow-up | 11–12 (GDC–CDR), 5 (CDR–cBioPortal) | Differences −3,262 to +417 days (GDC–CDR) and 3–365 days (CDR–cBioPortal). | Unexplained. |
| Mucinous histology, CDR text vs ICD-O-3 | 7 | e.g. "Colon Mucinous Adenocarcinoma" coded 8140/3 (3 patients). | Coding contradiction or coarse text; minor. |
| Stage present in another source but missing in GDC | 1 (vs CDR), 6 (vs cBioPortal) | See 3.2. | To examine at the staging stage. |

---

## 5. Implications for the GDC-based OMOP ETL

### 5.1 Supported by verified findings [V]

- GDC Release 46.0 is the only complete patient set (633). It contains every patient in the other two sources, so using it as the ETL source loses no one.
- cBioPortal membership must **not** be used as an inclusion criterion: it would remove 39 patients for molecular-QC reasons unrelated to clinical eligibility.
- Colon vs rectum must be defined from ICD-O-3 topography, **not** the COAD/READ project label.
- GDC `tissue_or_organ_of_origin` (not `site_of_resection_or_biopsy`) represents tumour site; GDC `morphology` provides ICD-O-3 histology.
- Exact day offsets are more reliable than derived ages in years (±1 year).
- Missing-value normalisation must be source-specific, and `[Discrepancy]` should be retained as a distinct flag rather than collapsed into missing.
- The hierarchical GDC JSON, not the flattened TSV, is the appropriate structural input. Clinical TSV: 2,443 rows for 633 patients. [V, Python exploration; R reproduction pending, see A.4]

### 5.2 Deferred decisions [D]: not made at this checkpoint

| Decision | Evidence from this checkpoint |
|---|---|
| **Index diagnosis selection** | 631/633 cases have exactly one `primary` diagnosis; 2 have a single unlabelled diagnosis; 182 cases have 2–5 diagnoses (prior primaries, metastases, recurrences, 128 "Unknown"). The proposed "Option B" (primary label + documented CRC-site fallback) is **not yet approved**. |
| **Observation-period design** | Follow-up definitions differ by source; 4 patients have no follow-up anywhere; 2 GDC follow-up fields exist. |
| **Date and birth-year derivation** | Ages in years are ±1; birth offsets differ by 1 day for obfuscated 90+ patients; the synthetic anchor date is not yet designed. |
| **Survival-source precedence** | CDR is the designated survival reference, but 7 GDC death-day contradictions and 1 apparently newer GDC death (F5-6702) are the test cases for any rule. |
| **Stage II/III eligibility** | Stage concordance is perfect where present; open issues are missing stage, CDR `[Discrepancy]`, 70 missing editions, mixed AJCC editions (5th/6th/7th), clinical vs pathological staging and possible neoadjuvant treatment. Staging variables have **not** yet been profiled in depth. |
| **Handling of flagged cases** | TCGA-A6-2670 (CDR-redacted, unlabelled diagnosis), TCGA-AF-3912 (GDC only, unlabelled diagnosis, no outcomes), the 6 cases with unusual primary sites, and the 9 survival contradictions. |

---

## Appendix A. Reproduction inventory

### A.1 Software (as tested)

R 4.3.3; packages `jsonlite` 1.8.8, `readxl` 1.4.3, `readr` 2.1.5, `dplyr` 1.1.4, `tidyr` 1.3.1, `openssl` 2.1.1. Run all scripts from the RStudio project root; every path is relative to it.

### A.2 Input files and expected locations (git-ignored `data/`)

| Local path | Origin | Size (bytes) | SHA-256 |
|---|---|---|---|
| `data/raw/gdc_r46/clinical/clinical.cohort.2026-10-29.json` | GDC portal, clinical JSON | 5,037,387 | `73dd7c78f39e5b4375ef2c516e81cfe034a9eb15c8567000afa1941e16c084f6` |
| `data/raw/gdc_r46/clinical/` ← extract of `clinical.cohort.2026-10-29.tar.gz` (clinical.tsv, follow_up.tsv, pathology_detail.tsv, family_history.tsv, exposure.tsv) | GDC portal, clinical TSV archive | 471,109 (archive) | `533a990f9c6017fa0bf40820faebdabc20627d43768a95300b1884f5f69bf90d` |
| `data/raw/gdc_r46/biospecimen/` ← extract of `biospecimen.cohort.2026-10-29.tar.gz` (sample.tsv, portion.tsv, analyte.tsv, aliquot.tsv, slide.tsv) | GDC portal, biospecimen TSV archive | 1,332,313 (archive) | `8c69ccce0fac660f4a84697aaeac1d07e8434fdc4592f26aa3adc74f891b7d65` |
| `data/raw/cdr/TCGA-CDR-SupplementalTableS1.xlsx` | GDC PanCanAtlas publication page | 2,945,129 | `ea594c0fbb6731477c7ac511fab449ca9c38b0d42d269591ed9f5c4090e75a5a` |
| `data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt` | Script 02 (below) | 175,194 | `3c5738f41b9047313af8f809c244bc66419991379b0c1978816dfc56c184c313` |
| `data/raw/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_sample.txt` | Script 02 (below) | 105,858 | `b9cbcc7cbc7471b44a3767aac48d4d0625001e3c07f8776c50f77364341f7817` |

GDC checksums identify the exact files analysed (the original downloads). A fresh portal download may not be byte-identical even within the same release. The `data_release` value (46.0) is the authoritative version identifier.

### A.3 Scripts (run in this order)

| # | Script | Inputs | Outputs (`results/exploration/`) |
|---|---|---|---|
| 1 | `R/01_acquisition/02_download_cbioportal_clinical.R` | Internet access to `media.githubusercontent.com` | The two cBioPortal files (A.2), verified by size and SHA-256; stops on mismatch |
| 2 | `R/02_exploration/03_source_patient_counts.R` | GDC JSON; CDR xlsx; cBioPortal patient + sample tables | `source_patient_counts.csv` (10 rows) |
| 3 | `R/02_exploration/04_source_patient_overlap.R` | GDC JSON; GDC `biospecimen/sample.tsv`; CDR xlsx; cBioPortal patient table | `04_overlap_patterns.csv` (4), `04_overlap_pairwise.csv` (3), `04_nesting_tests.csv` (3), `04_project_concordance.csv` (2), `04_check_cbio_missing_tumour_sample.csv` (2), `04_project_vs_primary_site.csv` (18), `04_discrepant_patients.csv` (42; IDs + membership flags only) |
| 4 | `R/02_exploration/05_source_variable_comparison.R` | GDC JSON; CDR xlsx; cBioPortal patient table | `05_variable_comparison_summary.csv` (30), `05_os_unit_check.csv` (1), `05_categorical_disagreements.csv` (1), `05_numeric_disagreements.csv` (8), `05_birth_vs_age_obfuscation.csv` (4), `05_histology_cdr_vs_gdc.csv` (22), `05_missing_value_markers.csv` (41), `05_flagged_contradictions.csv` (10; IDs + issue type only) |

Numbers in brackets are data rows (excluding the header). Outputs contain aggregate counts, categories and open-access TCGA barcodes only, with no patient-level clinical values.

**Quick verification after reproduction:** `source_patient_counts.csv` shows 633 / 629 / 594; `04_overlap_patterns.csv` shows 591 / 38 / 3 / 1; `04_nesting_tests.csv` shows 3 violations for "cBioPortal within TCGA-CDR"; `05_os_unit_check.csv` shows a median of 30.417 days per month.

### A.4 Not yet reproduced in R

The structural exploration of the GDC clinical files (Stage 1, Step 2a) was performed in Python during an exploratory session and has **not** yet been scripted in R:
- clinical TSV rows vs cases (2,443 rows for 633 patients; up to 32 rows per patient);
- diagnoses per case (1: 451; 2–5: 182), diagnosis classifications, primary-site and morphology distributions;
- treatments per case (0–32), treatment types and therapeutic agents;
- follow-ups per case (0–14), pathology details and molecular tests (MMR IHC, KRAS, BRAF).

An R script reproducing these (`R/02_exploration/01_gdc_structure.R` or similar) is the outstanding item from the original exploration plan.

### A.5 Bundle contents (`tcga-crc-omop_reconciliation_step1-3.zip`)

```
docs/source_reconciliation.md
R/01_acquisition/02_download_cbioportal_clinical.R
R/02_exploration/03_source_patient_counts.R
R/02_exploration/04_source_patient_overlap.R
R/02_exploration/05_source_variable_comparison.R
results/exploration/source_patient_counts.csv
results/exploration/04_overlap_patterns.csv
results/exploration/04_overlap_pairwise.csv
results/exploration/04_nesting_tests.csv
results/exploration/04_project_concordance.csv
results/exploration/04_check_cbio_missing_tumour_sample.csv
results/exploration/04_project_vs_primary_site.csv
results/exploration/04_discrepant_patients.csv
results/exploration/05_variable_comparison_summary.csv
results/exploration/05_os_unit_check.csv
results/exploration/05_categorical_disagreements.csv
results/exploration/05_numeric_disagreements.csv
results/exploration/05_birth_vs_age_obfuscation.csv
results/exploration/05_histology_cdr_vs_gdc.csv
results/exploration/05_missing_value_markers.csv
results/exploration/05_flagged_contradictions.csv
```

No raw data, vocabularies or credentials are included.
