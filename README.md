# Code and data for: Phage–antibiotic synergy in ESKAPE pathogens — a scoping review of antibiotic-partner selection

This repository contains the extraction dataset, derived data tables and analysis code that reproduce every number, table and statistic reported in the manuscript.

**It does not contain verbatim article text.** Abstracts and full texts retrieved during screening are third-party copyrighted material and are therefore not redistributed here; the screening audit file records each record's bibliographic identifiers, relevance score and scoring rationale, which is sufficient to reproduce the screening decisions and re-retrieve any record from Europe PMC.

---

## Repository contents

```
data/
  eskape_pas_dataset_clean.csv      419 x 31  main extraction dataset (1 row = 1 study x antibiotic arm)
  crkp_pas_dataset_clean.csv         74 x 41  CRKP continuity-check dataset (35 studies)
  screening_audit_4402.csv         4,402 rows screening decisions with score + rationale, no abstracts
  eligibility_exclusions_97.csv        97 rows eligibility exclusions with reason
  table1_corrected.json                     Table 1 values, ICC, design effect, effective n, cluster CI, McNemar, GEE
  within_study_paired.csv              23 rows the 23 head-to-head paired studies
  final_references.json                30 rows verified references with PMID/DOI

scripts/
  06_eskape_search.py                       search of Europe PMC (6 organism queries, 2015-2026)
  07_eskape_triage.py                       weighted relevance scoring and threshold
  09_eskape_extract.py                      structured extraction (includes antibiotic label vocabulary)
  18_reconcile_and_reanalyse.py             Table 1 reconciliation, within-study design
  19_corrected_table1_and_design.py         Table 1, ICC, design effect, McNemar, GEE  <- reproduces all statistics
  validate_results.py                       NEW: asserts every headline number in the manuscript

figures/
  prisma_scr_flow.svg / .png                PRISMA-ScR flow diagram

environment_pip_freeze.txt                  Python environment
```

---

## Quick start — verify the headline numbers in one command

```bash
pip install -r requirements.txt
python scripts/validate_results.py
```

`validate_results.py` recomputes and asserts **every number reported in the manuscript** — the 272/419/157 counts, all nine Table 1 rows with Wilson CIs, the ICC and design effect, all Fisher/McNemar/GEE results, and the organism and assay breakdowns. It prints a PASS/FAIL line per check and exits non-zero on any mismatch. This is the fastest way for a reviewer or reader to confirm the paper's arithmetic.

Expected output ends with:

```
96 checks: 93 PASS, 0 FAIL, 3 WARN
```

(The three WARNs are documented discrepancies — see `docs/KNOWN_DISCREPANCIES.md`.)

## Full reproduction

```bash
# screening (requires network access to Europe PMC)
python scripts/06_eskape_search.py          # regenerates the record pool
python scripts/07_eskape_triage.py          # regenerates the score and the >=14 selection

# statistics (offline, uses data/ only)
python scripts/19_corrected_table1_and_design.py
python scripts/18_reconcile_and_reanalyse.py
```

Note that re-running the search will return a different number of records unless run on the same date, because Europe PMC is continuously updated. The search date used for the manuscript is recorded in `docs/METHODS_SPECIFICATION.md`.

---

## Dataset structure

`eskape_pas_dataset_clean.csv` — one row per antibiotic arm within one study.

| Column group | Columns | Notes |
|---|---|---|
| Identity | `PMID_or_DOI`, `DOI`, `Publication_Year`, `Study_Title`, `Journal` | `PMID_or_DOI` is the clustering unit ("study") throughout |
| Organism | `Organism` | one of the six ESKAPE groups |
| Antibiotic | `Antibiotic_Tested`, `Antibiotic_Standardised`, `ChEMBL_ID`, `Antibiotic_Class`, `Antibiotic_Mechanism` | class/mechanism assigned from a fixed vocabulary in `scripts/09_eskape_extract.py` |
| Phage | `Phage_Name`, `Phage_Family_or_Genus` | |
| Genotype | `Carbapenemase_Genotype*`, `CRKP_Strain_ST_and_Capsule` | |
| Exposure | `Concentration_or_MIC_Fraction`, `Dosing_Sequence` | |
| Assay | `Assay_Type`, `Assay_Types_All` | `Assay_Type` lists the primary assay first |
| Outcome | `Quantitative_Outcome`, `Outcome_Log10_Reduction`, `Outcome_FIC`, `Outcome_Percent`, `Interaction_Classification`, `Synergy_Outcome` | **`Synergy_Outcome`** is 1/0 for determinate calls, blank for indeterminate |
| Provenance | `Evidence_Snippet`, `Classification_Basis`, `Evidence_Level`, `Data_Quality_Flags` | `Evidence_Snippet` is a short quoted phrase supporting the classification |

**Analysis set:** rows with non-null `Synergy_Outcome` (n = 272 of 419; 105 study–organism entries).

---

## How the clustering and the design-based analysis work

Three subset definitions are used in the paper. They are easy to confuse, so they are stated precisely here; all three are implemented in `scripts/19`.

| Subset | Definition | n |
|---|---|---|
| **Pooled** | all rows with a determinate `Synergy_Outcome` | 272 arms / 105 study–organism entries |
| **Design-based** | studies with `nunique(Antibiotic_Mechanism) >= 2` — selected on *design*, never on outcome | 209 arms / 50 studies |
| **Outcome-selected** | studies with `nunique(Antibiotic_Tested) >= 2` **and** `0 < synergy_count < arm_count`. **Exploratory only; must not be used to estimate class performance** | 17 studies |
| **Head-to-head** | design-based studies measuring *both* cell-wall PBP and DNA-gyrase mechanisms | 23 studies |

The outcome-selected subset is retained in the paper purely to demonstrate the size of the selection effect (it inflates the fluoroquinolone–β-lactam gap from 8.2 to 40.7 percentage points). It is not an estimate of comparative class performance.

---

## Caveats we ask you to read before citing this work

1. **Screening was single-reviewer.** The weighted relevance score is deterministic and every decision is logged in `screening_audit_4402.csv`, but there was no second screener, no inter-rater kappa and no calibration of the score against a hand-screened reference set. The manuscript states this.
2. **The screening threshold is procedural.** A cut at `score >= 14` was chosen for tractability, not derived statistically. Cut at 12 would have yielded 403 records rather than 254. The sensitivity ladder is in `docs/METHODS_SPECIFICATION.md`.
3. **The ICC is estimator-dependent.** The manuscript reports 0.48 (variance-component estimator). A one-way ANOVA ICC(1) on the same data gives 0.15. A cluster bootstrap CI of 79.7–90.6% agrees with the reported 79.7–90.9%, so the substantive conclusion (the naive CI is too narrow) holds either way. See `docs/KNOWN_DISCREPANCIES.md`.
4. **Antibiotic class, mechanism and cidal/static labels are hard-coded**, not derived from an external pharmacology database. Reproduce them from `scripts/09_eskape_extract.py` or Table 1 will not reconcile.

---

## Publishing this repository (for the authors)

The repository ships with `publish.sh`, which validates the results, blocks any
copyrighted article text from being committed, initialises git and creates the
GitHub repository in one step.

```bash
bash publish.sh
```

Before publishing, complete the placeholders in `CITATION.cff` and `.zenodo.json`
and in the *Citation* section below (`[USERNAME]`, `[REPO-NAME]`, author names).

To mint a citable DOI, after pushing: Zenodo -> Settings -> GitHub -> enable this
repository, then cut a release on GitHub tagged `v1.0.0`.

## Licence and citation

- **Code:** MIT (see `LICENSE`).
- **Data files we generated** (`eskape_pas_dataset_clean.csv`, `crkp_pas_dataset_clean.csv`, `screening_audit_4402.csv`, `table1_corrected.json`, `within_study_paired.csv`): CC BY 4.0.
- **Third-party content is not redistributed.** Article abstracts and full texts remain the property of their publishers. `final_references.json` contains only bibliographic metadata (titles, DOIs, PMIDs), which is factual.

If you use this dataset or code, please cite the manuscript.
