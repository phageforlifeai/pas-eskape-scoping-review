# Code and data for: Phage–antibiotic synergy in ESKAPE pathogens — a scoping review of antibiotic-partner selection

This repository contains the extraction dataset, derived data tables and analysis code that reproduce every number, table and statistic reported in the manuscript.

**It does not contain verbatim article text.** Abstracts and full texts retrieved during screening are third-party copyrighted material and are therefore not redistributed here. The screening audit files record each record's bibliographic identifiers, relevance score and scoring rationale, which is sufficient to reproduce every screening decision and re-retrieve any record from Europe PMC.

---

## Quick start — verify every reported number in three commands

```bash
pip install -r requirements.txt
python scripts/validate_results.py
python scripts/validate_reanalysis.py
```

Expected output ends with:

```
RESULT: 178 PASS, 0 FAIL, 0 WARN
RESULT: 82 PASS, 0 FAIL
```

Both scripts recompute everything from the CSVs and assert every reported value — all nine Table 1 rows with Wilson confidence intervals, the ICC and design effect, the Fisher, McNemar and GEE results, the organism and assay breakdowns, and the full screening flow. Each prints a PASS/FAIL line per check and exits non-zero on any mismatch.

`pip install -r requirements.txt` must run first: the GEE checks depend on `statsmodels`, and if it is absent those checks are skipped with a printed notice rather than failing — you will see **77 PASS** instead of 82.

Two data files are worth naming up front, because the pipeline changed as the study set grew and the file names are historical:

| What you want | Read this | Note |
|---|---|---|
| The analysis dataset behind every reported statistic | `data/priority4_combined_dataset.csv` | 295 rows, one per determinate observation |
| The full extraction dataset | `data/eskape_pas_dataset_clean.csv` | 419 rows; 166 studies, 442 antibiotic arms |
| Screening scores for all 4,402 records | `data/triage_rescore_4402.csv` | no article text, re-retrievable from Europe PMC |

`scripts/21_reanalysis_priority4.py` and `scripts/22_gee_reanalysis.py` operate on `priority4_combined_dataset.csv`; `scripts/18_*` and `scripts/19_*` operate on `eskape_pas_dataset_clean.csv`. Some file and directory names retain earlier working titles (e.g. `eligibility_exclusions_97.csv`, which now holds 116 rows, and `priority4_*`, which is the label the sensitivity screen was given internally). They are kept so that the paths in `docs/` and in the scripts resolve, and they are listed here so nothing is ambiguous.

---

## Headline results reproduced here

| Quantity | Value |
|---|---|
| Records screened | 4,402 |
| Studies included | 166 |
| Antibiotic-level observations extracted | 442 |
| Observations with a determinate outcome | 295 (66.7%) |
| Pooled synergy proportion | 83.7% (247/295), 95% CI 79.1–87.5% |
| Antagonism | 5.4% (16/295), 95% CI 3.4–8.6% |
| ICC(1) / design effect / effective n | 0.157 / 1.25 / 235 |
| Design-based subset | 55 studies, 226 observations, 82.7% |
| Head-to-head paired studies | 26; b = 6, c = 1; exact McNemar P = 0.13 |
| GEE odds ratio (design-based subset) | 3.97 (95% CI 0.59–26.59), P = 0.156 |

---

## Repository contents

```
data/
  eskape_pas_dataset_clean.csv        419 x 31   full extraction dataset (1 row = 1 study x antibiotic arm;
                                                  442 extracted observations incl. indeterminate outcomes)
  priority4_combined_dataset.csv      295 x 10   analysis dataset: all determinate observations
  priority4_extracted_arms.csv           23 x 15 sensitivity-screen additions, full extraction schema
  priority4_screening_log_28.csv         28 x 11 every eligibility decision with the supporting evidence
  triage_rescore_4402.csv             4,402 rows screening scores and rationale, no article text
  screening_audit_4402.csv            4,402 rows screening decisions with score and rationale
  eligibility_exclusions_97.csv         116 rows eligibility exclusions with reason (filename kept as legacy)
  crkp_pas_dataset_clean.csv           74 x 41   CRKP continuity-check dataset (35 studies)
  table1_reanalysed.json                           Table 1 values and statistics
  reanalysis_results.json                          every reported statistic, machine-readable
  gee_results.json                                 GEE fit
  within_study_paired_26.csv/.json       26 rows the head-to-head paired studies
  final_references.json                  30 rows verified references with PMID/DOI

scripts/
  06_eskape_search.py                  search of Europe PMC (6 organism queries, 2015-2026)
  07_eskape_triage.py                  weighted relevance scoring and threshold
  07b_triage_sensitivity.py            re-score with environmental isolation-source terms removed
  09_eskape_extract.py                 structured extraction (includes the antibiotic label vocabulary)
  18_reconcile_and_reanalyse.py        Table 1 reconciliation, within-study design
  19_corrected_table1_and_design.py    Table 1, ICC, design effect, McNemar, GEE
  21_reanalysis_priority4.py           rebuilds every reported statistic from the analysis dataset
  22_gee_reanalysis.py                 refits the GEE and asserts the published fit reproduces
  validate_results.py                  asserts every reported value (178 checks)
  validate_reanalysis.py               asserts the analysis dataset end to end (82 checks)

docs/
  METHODS_SPECIFICATION.md             full methods specification
  KNOWN_DISCREPANCIES.md               estimand choices and implementation notes

figures/
  prisma_scr_flow.svg / .png           PRISMA-ScR flow diagram

environment_pip_freeze.txt             Python environment
```

---

## Full reproduction

```bash
# screening (requires network access to Europe PMC)
python scripts/06_eskape_search.py        # regenerates the record pool
python scripts/07_eskape_triage.py        # regenerates the score and the >=14 selection
python scripts/07b_triage_sensitivity.py  # regenerates the sensitivity re-score

# statistics (offline, uses data/ only) — run in this order
pip install -r requirements.txt                # statsmodels is required for the GEE step
python scripts/19_corrected_table1_and_design.py
python scripts/18_reconcile_and_reanalyse.py
python scripts/22_gee_reanalysis.py            # writes data/gee_results.json
python scripts/21_reanalysis_priority4.py      # rebuilds the master results JSON

# either order of the last two now gives the same numbers: 21 re-injects the GEE
# fits from gee_results.json, so running 21 after 22 does not drop them
python scripts/22_gee_reanalysis.py
```

Re-running the search will return a different number of records unless run on the same date, because Europe PMC is continuously updated. The search date used for the manuscript is **20 November 2026**, recorded in `docs/METHODS_SPECIFICATION.md`.

`07b_triage_sensitivity.py` requires `data/eskape_pool.json`, which embeds record abstracts and is therefore not distributed. It exits with a clear notice if the file is absent; the committed `data/triage_rescore_4402.csv` holds the resulting scores and rationales without the abstracts.

---

## Dataset structure

`priority4_combined_dataset.csv` — one row per antibiotic arm within one study; this is the analysis set for the manuscript.

| Column | Content |
|---|---|
| `study` | record identifier (PMID or DOI); the clustering unit throughout |
| `Organism` | one of the six ESKAPE pathogen groups |
| `Antibiotic_Name`, `Antibiotic_Class` | the agent tested |
| `Antibiotic_Mechanism` | mechanism class from the fixed vocabulary in `scripts/09_eskape_extract.py` |
| `Bactericidal_Static` | cidal or static |
| `Assay` | primary assay type |
| `Dosing_Sequence` | simultaneous, phage-first or antibiotic-first |
| `Interaction_Classification` | synergy, indifference, additive or antagonism |
| `Synergy_Outcome` | 1/0 used in all binomial analyses |

The fuller original extraction file `eskape_pas_dataset_clean.csv` carries 31 columns including phage identity, resistance genotype, concentration relative to MIC, quantitative outcome values and evidence snippets.

---

## How the analyses are defined

| Subset | Definition | n |
|---|---|---|
| **Pooled** | all rows with a determinate outcome | 295 arms / 113 study–organism entries |
| **Design-based** | studies with `nunique(Antibiotic_Mechanism) >= 2` — selected on *design*, never on outcome | 226 arms / 55 studies |
| **Outcome-selected** | studies testing ≥ 2 antibiotics with `0 < synergy_count < arm_count`. **Exploratory only** | 21 studies |
| **Head-to-head** | studies measuring *both* cell-wall PBP and DNA-gyrase mechanisms | 26 studies |

The outcome-selected subset is reported solely to demonstrate the size of the selection effect — it inflates the fluoroquinolone–β-lactam gap from 9.3 to 37.1 percentage points. It is not an estimate of comparative class performance.

The GEE is fitted on the **design-based subset** restricted to the two mechanism classes, with study as the clustering variable, a binomial/logit link, an exchangeable working correlation and robust sandwich variance. `scripts/22_gee_reanalysis.py` asserts that this specification reproduces the reported estimate before printing it.

---

## Caveats

1. **The ICC is a one-way ANOVA ICC(1), not a variance ratio.** The design-effect formula `1 + (m̄ − 1)·ICC` is only valid for a genuine intraclass correlation. A variance-partition statistic (η² = 0.48) would give a design effect of 1.76; the correct ICC(1) = 0.157 gives 1.25. Only the latter is reported. See `docs/KNOWN_DISCREPANCIES.md` §2.

2. **Clustering is detectable but modest.** Under independence the expected η² is (G−1)/(n−1) = 0.381, so the observed 0.477 is only ~9.6 percentage points above chance. A permutation test gives P = 0.020. The ICC should not be read as evidence of strong within-study correlation. See §3.

3. **Screening was single-reviewer.** The score is deterministic and every decision is logged, but there was no second screener and no inter-rater kappa. Post hoc, the score discriminated between included and excluded records within the prioritised set (AUC 0.747), and the yield of eligible studies declined monotonically from 100% at scores of 22–24 to 39% at the threshold of 14. See §1 and §4.

4. **The screening threshold is procedural.** A cut at `score >= 14` was chosen for tractability, not derived statistically. Raising it to 16 would have excluded 42 included studies; lowering it to 12 would have added 149 records to screen. The sensitivity ladder is in `docs/METHODS_SPECIFICATION.md`.

5. **Environmental isolation-source terms were excluded from the off-target penalty.** Phage provenance sentences are routine in abstracts, so penalising `sewage` and `wastewater` would suppress on-topic records. §1 records the effect of this decision and the records it changed.

6. **Antibiotic class, mechanism and cidal/static labels are hard-coded**, not derived from an external pharmacology database. Reproduce them from `scripts/09_eskape_extract.py` or Table 1 will not reconcile.

7. **Organism assignment uses regular-expression matching** over title and abstract rather than taxonomic lookup. Consequences and the one record affected are documented in §11.

---

## Licence and citation

- **Code:** MIT (see `LICENSE`).
- **Data files we generated:** CC BY 4.0.
- **Third-party content is not redistributed.** Article abstracts and full texts remain the property of their publishers. `final_references.json` contains only bibliographic metadata (titles, DOIs, PMIDs), which is factual.

If you use this dataset or code, please cite the manuscript.
