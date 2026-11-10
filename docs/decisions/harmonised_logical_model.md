# TCGA COAD/READ — Harmonised Intermediate Data Model

**Project:** TCGA COAD/READ → OMOP  
**Phase:** 4 — Harmonised Intermediate Data Model  
**Status:** Conceptual design agreed; physical schemas and transformations pending  
**Initial decision log date:** 2026-11-09  
**Primary ETL source:** GDC Release 46.0  
**Supplementary comparison sources:** TCGA-CDR (curated survival reference); cBioPortal PanCancer Atlas (comparison only)

## 1. Purpose and boundaries

This document records the team's *intermediate-model design decisions*. The model reorganises GDC clinical and biospecimen records into explicit entities while preserving the grain, native identifiers, lineage and parent–child relationships. It is **not** an OMOP CDM mapping specification, and names here do not imply final OMOP target tables or standard vocabulary mappings.

The GDC release and extraction metadata are held at the dataset/extract level rather than repeated on each patient row. No raw patient-level source data should be committed to GitHub. The source files remain unchanged; harmonised outputs are derived from them.

## 2. Agreed principles

1. **Source-faithful grain:** Represent each original entity or event at its actual source-record level; do not flatten away repeated records.
2. **Relationship fidelity:** Preserve GDC's observed parent–child relations. Do not fabricate diagnosis links for case-level follow-ups.
3. **Identifiers and provenance:** Retain GDC UUIDs and relevant TCGA submitter barcodes when present. Generated keys may be introduced only where necessary, with a deterministic documented rule.
4. **No premature OMOP decisions:** OMOP target tables, concept mappings and event semantics are decided in a later phase.
5. **Separation of concerns:** Source discovery/profiling, intermediate modelling, implementation and OMOP mapping are distinct activities.
6. **Metadata normalisation:** Store GDC source and release metadata once per extraction, not as constant values in each patient row.

## 3. Decision register

| ID | Entity / relationship | Agreed conceptual grain and linkage | Status |
|---|---|---|---|
| HIM-001 | `patient` | One row per unique GDC case; incorporate one-to-one `demographic` attributes | Agreed |
| HIM-002 | `diagnosis` | One row per GDC diagnosis record; foreign key to `patient` | Agreed |
| HIM-003 | `treatment` | One row per source treatment record; foreign key to `diagnosis` | Agreed |
| HIM-004 | `follow_up` | One row per source follow-up record; foreign key directly to `patient`, not an inferred diagnosis | Agreed |
| HIM-005 | `family_history` | One row per source family-history record; foreign key to `patient` | Agreed |
| HIM-006 | `exposure` | One row per source exposure record; foreign key to `patient` | Agreed |
| HIM-007 | `pathology_detail` | One row per source pathology-detail record; foreign key to `diagnosis` | Agreed |
| HIM-008 | `molecular_test` | One row per source clinical molecular-test record; foreign key to `follow_up` | Agreed |
| HIM-009 | `other_clinical_attribute` | One row per source other-clinical-attribute record; foreign key to `follow_up`; preserve timepoint | Agreed |
| HIM-010 | `sample` | One row per unique GDC sample UUID; foreign key to `patient` | Agreed |
| HIM-011 | `portion` | One row per unique GDC portion UUID; foreign key to `sample` | Agreed |
| HIM-012 | `analyte` | One row per unique GDC analyte UUID; foreign key to `portion` | Agreed |
| HIM-013 | `aliquot` | One row per unique GDC aliquot UUID; foreign key to `analyte` | Agreed |
| HIM-014 | `slide` | One row per unique slide identifier | Agreed |
| HIM-015 | `slide_portion` | One row per distinct (`slide_id`, `portion_id`) association; composite key; no association deduplication | Agreed |
| HIM-016 | Dataset/extract metadata | Keep source system, release and extract metadata separately from patient rows | Agreed |

**Status wording:** “Agreed” means the conceptual choice was accepted in discussion, **not** that the table exists or its physical field schema has been validated.

## 4. Conceptual entity relationships

```text
patient  [incorporates demographic]
├── diagnosis
│   ├── treatment
│   └── pathology_detail
├── follow_up
│   ├── molecular_test
│   └── other_clinical_attribute
├── family_history
├── exposure
└── sample
    └── portion
        ├── analyte
        │   └── aliquot
        └── slide_portion ── slide
```

All branches represent candidate one-to-many parent/child relationships except the slide–portion linkage, which is modelled with a junction table. This diagram shows agreed conceptual links; actual permissible cardinalities and key constraints should be checked against existing profiling results when writing the physical schema.

## 5. Source evidence already established

The following counts and integrity summaries were verified in the exploration output files provided for this project:

| Finding | Value | Evidence file |
|---|---:|---|
| GDC cases | 633 | `results/exploration/01_gdc_structure_summary.csv` |
| GDC diagnosis records | 878 | `results/exploration/01_gdc_structure_summary.csv` |
| Cases with multiple diagnoses | 182 | `results/exploration/01_gdc_structure_summary.csv` |
| GDC treatment records | 2,354 | `results/exploration/01_gdc_structure_summary.csv` |
| Samples / distinct sample IDs | 2,037 | `results/exploration/02_biospecimen_entity_summary.csv` |
| Portions / distinct portion IDs | 3,217 | Same |
| Analytes / distinct analyte IDs | 4,436 | Same |
| Aliquots / distinct aliquot IDs | 12,574 | Same |
| Slide association rows | 2,035 | `results/exploration/02_biospecimen_slide_summary.csv` |
| Distinct slide IDs | 1,998 | Same |
| Slide IDs linked to multiple portions | 35 | `results/exploration/02_biospecimen_relationship_summary.csv` |
| Duplicate slide–portion pairs | 0 | Same |
| Clinical cases linked to biospecimen cases | 633 of 633 | `results/exploration/02_clinical_biospecimen_linkage.csv` |

### Slide–portion resolution

The 2,035 source slide rows reflect 2,035 unique slide–portion associations for 1,998 distinct slides. Thirty-five slide identifiers have multiple parent portions (maximum three). There were no duplicate slide–portion pairs and no reported attribute conflicts for repeated slide IDs. **Do not drop these valid associations to force a one-parent-per-slide design.**

## 6. Important semantic safeguards

- **Demographics:** Fold the single case-level demographic object into `patient`; retain useful source IDs and attributes after checking names and completeness.
- **Clinical molecular tests:** Preserve test records, including non-genetic tests. Do not infer that every result is somatic or germline; variant origin requires explicit supporting evidence.
- **Timepoints:** `other_clinical_attribute` may refer to periods such as initial or prior to diagnosis, even though the source objects are nested under `follow_up`. Preserve the source timepoint category.
- **Treatments:** Preserve source treatment records and original treatment-status indicators. A recorded treatment object does not automatically prove treatment administration; investigate meaning before OMOP mapping.
- **Survival:** Keep TCGA-CDR as a curated comparison/reference, not an untracked overwrite of GDC clinical fields.

## 7. Unresolved items before implementing tables

- Confirm the **exact source field paths** and physical column names for every entity using the existing JSON schema inventory and biospecimen outputs.
- Confirm UUID availability, nullability, uniqueness and potential generated-key requirements for nested records.
- Verify physical key and foreign-key constraints from existing exploration outputs; document any exceptions.
- Choose data types and missing-value treatment; distinguish absent objects from present-but-unknown fields.
- Define deterministic extraction logic and validation tests (source counts, uniqueness, referential integrity, preservation of arrays and associations).
- Determine which source attributes remain in parent records versus related child entities.
- Decide whether cohort-specific source fields need extension tables; do not force every source field into the initial physical schema.
- Specify output location and Git ignore rules for derived patient-level data.

## 8. Implementation and Git policy

- Exploration scripts: `R/02_exploration/` (existing work retained).
- Intermediate-model specification: this file under `docs/`.
- Source model and field-path catalogue: future `docs/gdc_source_model.md` (create only after deriving it from confirmed source evidence).
- Planned transformation scripts: a separate subsequent phase; do not write them until schemas and validations are agreed.
- Never commit `data/raw/`, extracts containing patient-level records, or sensitive-derived outputs. Verify `git status` and `git diff --cached` before committing.

## 9. Change log

| Date | Change |
|---|---|
| 2026-11-09 | Initial Phase 4 conceptual model and decision register assembled from agreed discussion; evidence cross-referenced to existing exploration outputs. |
