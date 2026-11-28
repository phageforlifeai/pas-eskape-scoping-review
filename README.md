# Code and data for: Phage–antibiotic synergy in ESKAPE pathogens — a scoping review of antibiotic-partner selection

This repository contains the extraction dataset, derived data tables and analysis code that reproduce every number, table and statistic reported in the manuscript.

**It does not contain verbatim article text.** Abstracts and full texts retrieved during screening are third-party copyrighted material and are therefore not redistributed here; the screening audit files record each record's bibliographic identifiers, relevance score and scoring rationale, which is sufficient to reproduce the screening decisions and re-retrieve any record from Europe PMC.

**Version:** v1.1.0 — reanalysed after a post hoc sensitivity screen (see §0).

---

## 0. Read this first: two validated versions

The manuscript reports **166 studies**. An earlier version reported 157. The difference comes from a sensitivity screen, and both versions are retained so the arithmetic behind either can be verified.

| | Published version | Reanalysed version (current) |
|---|---|---|
| Studies included | 157 | **166** |
| Observations extracted | 419 | **442** |
| Determinate observations | 272 | **295** |
| Pooled synergy | 85.3% | **83.7%** |
| Antagonism | 5.5% | **5.4%** |
| Validated by | `validate_results.py` → **178 PASS** | `validate_reanalysis.py` → **82 PASS** |

**What happened.** The weighted relevance score penalised "off-target" topics, and the term list included `sewage` and `wastewater`. Because phage provenance sentences ("phage X was isolated from hospital wastewater") are routine in abstracts of on-topic papers, this suppressed relevant records. Removing the two terms and re-scoring all 4,402 records raised the prioritised set from 254 to 282; screening the extra 28 found 9 further eligible studies contributing 23 observations.

**No conclusion changed.** Full reasoning in `docs/KNOWN_DISCREPANCIES.md` §1 and §13.

---

## Quick start — verify both versions in two commands

```bash
pip install -r requirements.txt
python scripts/validate_results.py       # 157-study version
python scripts/validate_reanalysis.py    # 166-study version
```

Expected output ends with:

```
RESULT: 178 PASS, 0 FAIL, 0 WARN
RESULT: 82 PASS, 0 FAIL
```

Both scripts recompute everything from the CSVs and assert every reported number — Table 1 rows with Wilson CIs, the ICC and design effect, Fisher/McNemar results, the organism and assay breakdowns, and the full screening flow. Each prints a PASS/FAIL line per check and exits non-zero on any mismatch. Places where a naive recalculation deliberately differs are printed as `[INFO]` lines and explained in `docs/KNOWN_DISCREPANCIES.md`.

---

## Repository contents

```
data/
  eskape_pas_dataset_clean.csv        419 x 31   original extraction dataset (1 row = 1 study x antibiotic arm)
  priority4_combined_dataset.csv      295 x 10   analysis dataset for the 166-study version (published 272 + 23 new)
  priority4_extracted_arms.csv           23 x 15 the 23 arms added by the sensitivity screen, full extraction schema
  priority4_screening_log_28.csv         28 x 11 every screening decision for the 28 added records, with evidence
  crkp_pas_dataset_clean.csv           74 x 41   CRKP continuity-check dataset (35 studies)
  screening_audit_4402.csv           4,402 rows  screening decisions with score + rationale, no abstracts
  triage_rescore_4402.csv             4,402 rows  re-scored triage after removing the isolation-source penalty, no abstracts
  eligibility_exclusions_97.csv          97 rows  eligibility exclusions, published screen
  table1_corrected.json                            Table 1 values and statistics, 157-study version
  table1_reanalysed.json                           Table 1 values and statistics, 166-study version
  reanalysis_results.json                          every revised statistic, machine-readable
  gee_results.json                                 GEE fit, published vs reanalysed
  within_study_paired.csv                23 rows  head-to-head paired studies, published
  within_study_paired_26.csv/.json       26 rows  head-to-head paired studies, reanalysed
  final_references.json                  30 rows  verified references with PMID/DOI

scripts/
  06_eskape_search.py                  search of Europe PMC (6 organism queries, 2015-2026)
  07_eskape_triage.py                  weighted relevance scoring and threshold
  07b_triage_sensitivity.py            re-score with isolation-source terms removed (sensitivity screen)
  09_eskape_extract.py                 structured extraction (includes antibiotic label vocabulary)
  18_reconcile_and_reanalyse.py        Table 1 reconciliation, within-study design
  19_corrected_table1_and_design.py    Table 1, ICC, design effect, McNemar, GEE
  21_reanalysis_priority4.py           rebuilds every statistic on the 166-study dataset
  22_gee_reanalysis.py                 refits the GEE; asserts the published fit reproduces
  validate_results.py                  asserts the 157-study numbers   (175 checks)
  validate_reanalysis.py               asserts the 166-study numbers   (77 checks)

docs/
  METHODS_SPECIFICATION.md             full methods spec + reanalysis addendum
  KNOWN_DISCREPANCIES.md               13 documented estimand/implementation choices

figures/
  prisma_scr_flow.svg / .png           PRISMA-ScR flow diagram

environment_pip_freeze.txt             Python environment
```

---

## Full reproduction

```bash
# screening (requires network access to Europe PMC)
python scripts/06_eskape_search.py       # regenerates the 4,402-record pool
python scripts/07_eskape_triage.py       # regenerates the score and the >=14 selection
python scripts/07b_triage_sensitivity.py # regenerates the sensitivity re-score

# statistics (offline, uses data/ only)
python scripts/19_corrected_table1_and_design.py
python scripts/18_reconcile_and_reanalyse.py
python scripts/21_reanalysis_priority4.py
python scripts/22_gee_reanalysis.py
```

Re-running the search will return a different number of records unless run on the same date, because Europe PMC is continuously updated. The search date used for the manuscript is **20 November 2026** and is recorded in `docs/METHODS_SPECIFICATION.md`.

`07b_triage_sensitivity.py` needs `data/eskape_pool.json`, which embeds record abstracts and is therefore not distributed. It exits with a clear `[SKIP]` message if the file is absent; the committed `data/triage_rescore_4402.csv` holds the resulting scores and rationales without the abstracts.

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

`priority4_combined_dataset.csv` is the analysis set for the reanalysed manuscript: the same columns harmonised, plus the 23 added arms, one row per antibiotic arm.

**Analysis sets:** rows with non-null `Synergy_Outcome` — 272 of 419 (105 study–organism entries) in the original; 295 of 442 (113 entries) in the reanalysis.

---

## How the clustering and the design-based analysis work

| Subset | Definition | 157-version | 166-version |
|---|---|---|---|
| **Pooled** | all rows with a determinate `Synergy_Outcome` | 272 arms / 105 entries | 295 arms / 113 entries |
| **Design-based** | studies with `nunique(Antibiotic_Mechanism) >= 2` — selected on *design*, never on outcome | 209 arms / 50 studies | 226 arms / 55 studies |
| **Outcome-selected** | studies with `nunique(Antibiotic_Tested) >= 2` **and** `0 < synergy_count < arm_count`. **Exploratory only** | 17 studies | 21 studies |
| **Head-to-head** | studies measuring *both* cell-wall PBP and DNA-gyrase mechanisms | 23 studies | 26 studies |

The outcome-selected subset is retained purely to demonstrate the size of the selection effect — it inflates the fluoroquinolone–β-lactam gap from 9.3 to 37.1 percentage points. It is not an estimate of comparative class performance.

---

## Caveats we ask you to read before citing this work

1. **The ICC is an ICC(1), not η².** An earlier version reported "ICC = 0.48" computed as `SS_between/SS_total`, which is η², not an intraclass correlation, and is invalid in the design-effect formula. The correct one-way ANOVA **ICC(1) = 0.157**, giving design effect 1.25 and effective n 235 (vs 1.76 and 154 under the incorrect estimator). The substantive conclusion — the naive CI is too narrow — is unchanged, and the corrected interval now agrees closely with the cluster bootstrap. See `KNOWN_DISCREPANCIES.md` §2.

2. **Clustering is detectable but modest.** The null expectation of η² is `(G−1)/(n−1) = 0.381`, so the observed 0.477 is only ~9.6 pp above chance. Permutation P = 0.020. Do not read the ICC as evidence of strong within-study correlation. See §3.

3. **Screening was single-reviewer.** The weighted score is deterministic and every decision is logged, but there was no second screener, no inter-rater kappa and no prospective calibration. Post hoc, the score discriminated between included and excluded records within the prioritised set (AUC 0.747), and eligible-study yield fell monotonically from 100% at scores of 22–24 to 39% at the threshold of 14. See §1.

4. **The screening threshold is procedural.** A cut at `score >= 14` was chosen for tractability. Raising it to 16 would have excluded 42 of 157 included studies; lowering it to 12 would have added 149 records to screen. The ladder is in `docs/METHODS_SPECIFICATION.md`.

5. **One component of the screen was found to be miscalibrated and was corrected** — the environmental isolation-source penalty. This is why the manuscript reports 166 rather than 157 studies. See §1 and §13.

6. **The GEE is fitted on the design-based subset, not the paired subset.** Refit after the reanalysis: **OR 3.97 (0.59–26.59), P = 0.156** (previously 4.6 [0.45–47.3], P = 0.20). Fitting on the 26 paired studies instead gives a different estimate and a narrower interval, because that subset is outcome-selected. Both intervals are very wide — five or six discordant studies carry all the information. See §7.

7. **Antibiotic class, mechanism and cidal/static labels are hard-coded**, not derived from an external pharmacology database. Reproduce them from `scripts/09_eskape_extract.py` or Table 1 will not reconcile.

---

## Publishing this repository (for the authors)

The repository ships with `publish.sh`, which validates the results, blocks any copyrighted article text from being committed, initialises git and creates the GitHub repository in one step.

```bash
bash publish.sh
```

Before publishing, complete the placeholders in `CITATION.cff` and `.zenodo.json` and in the *Citation* section below (the repository URL and author names).

To mint a citable DOI, after pushing: Zenodo → Settings → GitHub → enable this repository, then cut a release on GitHub tagged `v1.1.0`.

## Licence and citation

- **Code:** MIT (see `LICENSE`).
- **Data files we generated** (`eskape_pas_dataset_clean.csv`, `priority4_combined_dataset.csv`, `priority4_extracted_arms.csv`, `priority4_screening_log_28.csv`, `crkp_pas_dataset_clean.csv`, `screening_audit_4402.csv`, `triage_rescore_4402.csv`, `table1_corrected.json`, `table1_reanalysed.json`, `reanalysis_results.json`, `gee_results.json`, `within_study_paired*.csv`): CC BY 4.0.
- **Third-party content is not redistributed.** Article abstracts and full texts remain the property of their publishers. `final_references.json` contains only bibliographic metadata (titles, DOIs, PMIDs), which is factual.

If you use this dataset or code, please cite the manuscript.
