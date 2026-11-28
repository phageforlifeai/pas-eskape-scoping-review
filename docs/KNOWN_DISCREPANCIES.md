# Known discrepancies and estimand choices

This file records every place where a naive recalculation will **not** reproduce a reported number, and why. If you are checking our arithmetic, read this first — the difference is usually deliberate and documented here rather than an error.

Last updated: 27 November 2026 (reanalysis after the sensitivity screen).

---

## 1. The manuscript now reports 166 studies, not 157

The manuscript was revised after a **post hoc sensitivity screen**. A component of the weighted relevance score penalised "off-target" topics, and that term list included the environmental isolation-source words `sewage` and `wastewater`. Because phage provenance sentences ("phage X was isolated from hospital wastewater") are routine in abstracts of on-topic papers, this penalty suppressed otherwise relevant records.

Removing those two terms and re-scoring all 4,402 records raised the prioritised set from 254 to 282. Screening the 28 additional records found **9 further eligible studies contributing 23 antibiotic-level observations**.

| | Published version | Reanalysed version |
|---|---|---|
| Studies included | 157 | **166** |
| Observations extracted | 419 | **442** |
| Determinate observations | 272 | **295** |
| Pooled synergy | 85.3% | **83.7%** |
| Antagonism | 5.5% | **5.4%** |

**Both versions are retained and both are validated:**

- `scripts/validate_results.py` asserts the 157-study numbers → **175 PASS**
- `scripts/validate_reanalysis.py` asserts the 166-study numbers → **77 PASS**

The analysis dataset for the 166-study version is `data/priority4_combined_dataset.csv`. The original 419-row extraction file is untouched.

Full audit trail: `data/priority4_screening_log_28.csv` (all 28 decisions with evidence), `data/priority4_extracted_arms.csv` (the 23 arms), `data/triage_rescore_4402.csv` (the full re-score).

---

## 2. η² is not an ICC — this is the most important entry

An earlier version of this repository and the manuscript reported an "ICC of 0.48" computed as:

```
Σ nᵢ(pᵢ − p̄)² / [(n − 1) · p̄(1 − p̄)]
```

The numerator is the between-study sum of squares. The denominator is `(n−1) × p̄(1−p̄)`, and since `n·p̄(1−p̄)` **is** the total sum of squares, this expression evaluates to `SS_between / SS_total` — that is, **η²**, the proportion of total variance explained by study membership.

Verified on the data:

```
SS_between = 16.3224
SS_total   = 34.1176
SS_between / SS_total           = 0.4784   <- η²
SS_between / ((n-1) p̄(1-p̄))    = 0.4802   <- the 0.48 that was reported
```

**η² is not an intraclass correlation and is not valid in the design-effect formula.** The manuscript now reports a one-way ANOVA ICC(1):

```
ICC(1) = (MS_between − MS_within) / (MS_between + (n₀ − 1)·MS_within) = 0.155  (157-study version)
                                                                     = 0.157  (166-study version)
```

Consequences:

| Quantity | With η² (wrong) | With ICC(1) (correct) |
|---|---|---|
| Design effect | 1.76 | **1.25** |
| Effective n | 154 | **218** (157-study) / **235** (166-study) |
| Cluster-adjusted CI | 79.7–90.9 | **80.6–90.0** (157-study) / **79.0–88.4** (166-study) |

**The substantive conclusion is unchanged:** the cluster-adjusted interval remains wider than the unadjusted Wilson interval, so the claim that a naive CI overstates precision still holds. The corrected interval also agrees closely with the cluster bootstrap (79.8–90.6%), which the η²-based interval did not.

`validate_results.py` now computes ICC(1) and prints η² on an `[INFO]` line explicitly labelled descriptive-only.

---

## 3. Clustering is detectable but modest — do not read 0.48 (or 0.16) as "strong"

Under complete independence with these fixed cluster sizes, η² has an expected value of `(G−1)/(n−1)`, not zero:

```
(G−1)/(n−1) = 104/271 = 0.384
observed    =                  0.478
```

So the observed value is only ~9.5 percentage points above what chance alone produces. A permutation test (20,000 shuffles of outcome across studies, cluster sizes held fixed) gives **P = 0.032** (157-study) and **P = 0.020** (166-study).

Clustering is real but marginal, and 50 of 105 clusters (52 of 113 in the reanalysis) are singletons that carry no information about the ICC.

---

## 4. Cluster size range is 1–10, not 1–8

An earlier Supporting Information draft said "range 1–8". The actual distribution:

| Arms per study | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| Studies | 50 | 16 | 9 | 10 | 9 | 6 | 2 | 1 | 2 |

The maximum is 10.

---

## 5. Three different "gap" numbers, three different subsets

The fluoroquinolone–β-lactam difference appears as several values. They are not interchangeable:

| Value | Subset | Meaning |
|---|---|---|
| 9.3 pp | pooled | arm-level difference across all determinate observations |
| 12.0 pp | design-based, study-level mean | each study given one vote |
| 13.1 pp | design-based, arm-based | arm-level within the design-based subset |
| 37.1 pp | outcome-selected | **exploratory only** — selecting on a discordant outcome inflates the contrast |
| 16.7 pp | head-to-head paired | mean within-study difference of the 26 paired studies |

Citations of "the 40.7 pp contrast" in earlier drafts referred to the outcome-selected subset in the 157-study version; it is 37.1 pp after reanalysis.

---

## 6. Wash-out after adding observations

Several counts were derived from the 272-arm dataset and shift under the reanalysis. They are listed together here so nothing is silently stale:

| Item | Before | After |
|---|---|---|
| Studies testing ≥2 antibiotics | 55 | 57 |
| Discordant among those | 17 (30.9%) | 19 (33.3%) |
| Head-to-head paired studies | 23 | 26 |
| Paired b / c | 5 / 1 | 6 / 1 |
| McNemar exact P | 0.22 | 0.13 |
| Design-based studies / arms | 50 / 209 | 55 / 226 |
| Outcome-selected studies | 17 | 21 |
| Bactericidal vs static | 85.2% (243) vs 86.2% (29) | 83.5% (266) vs 86.2% (29) |
| β-lactam vs others | 84.8% (99) vs 85.5% (173), P=0.86 | 82.1% (112) vs 84.7% (183), P=0.63 |
| Envelope antagonism | 7/158 vs 8/114, P=0.42 | 8/174 vs 8/121, P=0.45 |
| Assay mix (checkerboard/biofilm/time-kill/in vivo) | 170/45/37/17 | 180/46/37/17 |

---

## 7. The GEE is fitted on the design-based subset, not the paired subset

This is the easiest place in the paper to mis-reproduce, so it is stated explicitly.

The GEE for the fluoroquinolone vs beta-lactam contrast is fitted on the **design-based subset** (studies testing >= 2 distinct mechanisms), restricted to the two mechanism classes of interest — **not** on the 26 paired studies.

```
published (157-study): 50 design-based studies, 105 arms (FQ 34, BL 71)
                       OR 4.62 (0.45-47.26), P = 0.197     [manuscript: 4.6 (0.45-47.3), P = 0.20]
reanalysed (166-study): 55 design-based studies, 117 arms (FQ 37, BL 80)
                       OR 3.97 (0.59-26.59), P = 0.156
```

Fitting instead on the 26 paired studies gives a **different** point estimate (~3.9) and a much narrower interval, because the paired subset is outcome-selected by construction. Switching estimands without saying so would break comparability with the previous version and would overstate precision.

Specification: Binomial family, logit link, exchangeable working correlation, robust (sandwich) variance, study as the clustering variable. Reproduce with `python scripts/22_gee_reanalysis.py`, which also asserts that its published-version fit matches the 4.6 (0.45-47.3) reported in the manuscript.

Both intervals are very wide. The estimate is imprecise because five or six discordant studies carry all the information; this is a limitation of the evidence base, not of the estimator.

## 8. The `no phage` gate is a substring match

The relevance score excludes any record with no `phag` substring in the title or abstract. That substring also matches `phagocytosis`, `phagocyte` and `macrophage`.

- 2,398 records scored −99 (54.5%).
- 2,384 genuinely contain no phage term — correct.
- **14** contain a phage substring but no real phage term. All 14 are true exclusions (the substring was `phagocyt*` / `macrophag*`), i.e. right decision, wrong reason.
- **55** records passed the gate on such an artefact alone, of which **5** reached the ≥14 threshold.

All were resolved at the eligibility stage, so the included set is unaffected — but the screen relies on the downstream eligibility assessment to correct for it.

---

## 9. Off-target penalties counted repeated terms twice

The weight table specifies "−3 per **distinct** term". In practice comma-joined duplicates were counted separately: `Biocontrol,biocontrol`, `Food,food`, `Sewage,sewage`, `Poultry,poultry`, `Plant,plant`, `Wastewater,wastewater`, `Biosensor,biosensor` — 24 records over-penalised by 3 points each.

This did not change the selected set (the threshold ladder at ≥12 and ≥10 shifted by 2 and 3 records respectively), but it is a documented divergence between the stated rule and the implementation. The reanalysis script fixes it with case-insensitive deduplication.

---

## 10. "Within-study" and "design-based" are not synonyms

Three subset definitions appear in the paper and are easy to confuse:

| Subset | Definition | 157-version | 166-version |
|---|---|---|---|
| **Pooled** | all rows with a determinate outcome | 272 arms / 105 entries | 295 arms / 113 entries |
| **Design-based** | `nunique(Antibiotic_Mechanism) >= 2` — selected on *design*, never outcome | 209 arms / 50 studies | 226 arms / 55 studies |
| **Outcome-selected** | `nunique(Antibiotic_Tested) >= 2` **and** `0 < synergy_count < arm_count` | 17 studies | 21 studies |
| **Head-to-head** | studies measuring **both** PBP and DNA-gyrase mechanisms | 23 studies | 26 studies |

The outcome-selected subset exists only to demonstrate the size of the selection effect and must never be used to estimate comparative class performance.

---

## 11. Organism assignment uses regex, not taxonomy

`_organisms` is populated by regular-expression matching over title and abstract, not by taxonomic lookup. Two consequences are documented:

- A record can match several organisms; the first match is used as the primary organism.
- **One sensitivity-screen record was excluded for exactly this reason.** Phage vB_KquU_φKuK6 is a *Klebsiella quasipneumoniae* / *K. variicola* phage that the regex labelled `K. pneumoniae` via the shared "Klebsiella" token. It was excluded on manual review because non-*pneumoniae* Klebsiella species are not ESKAPE organisms.

Any automated re-run of the screen will re-include that record. This is why the 28-record screening log is committed with explicit per-record decisions rather than left to the regex.

---

## 12. Assay counts are primary-assay, not exclusive

`Assay_Type` lists the primary assay first; 3 arms (157-version) / 15 arms (166-version) carry no assay label. The reported counts (170/45/37/17 → 180/46/37/17) are **primary-assay** counts and sum to fewer than the total number of arms. True non-exclusive counts are much larger (an arm can be both checkerboard and time-kill). Do not treat these as exclusive categories.

---

## 13. Why the reanalysis lowers the pooled estimate

This is not a discrepancy but a result worth recording, because it is the paper's central argument in miniature.

Adding 23 arms drawn from a **less-selected** pool moved:

- pooled synergy **down** 1.6 pp (85.3 → 83.7)
- β-lactam synergy **down** 2.7 pp (84.8 → 82.1)
- design-based synergy **down** only 2.0 pp (84.7 → 82.7)

The design-based estimate moved least. That is the pattern predicted if pooled proportions are inflated by selective reporting of positive combinations, and it independently corroborates the paper's methodological thesis.
