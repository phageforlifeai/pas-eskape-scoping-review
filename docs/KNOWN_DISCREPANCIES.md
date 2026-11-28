# Known discrepancies and estimand choices

This file records every place where a naive recalculation will **not** reproduce a reported number, and why. If you are checking our arithmetic, read this first — the difference is usually deliberate and documented here rather than an error.

Last updated: 28 November 2026.

---

## 1. How the final study set was reached

Screening used a weighted relevance score, described in full in `docs/METHODS_SPECIFICATION.md`. Two features of that score need to be stated explicitly because they explain the final counts.

**The off-target penalty excludes environmental isolation-source terms.** Phage provenance sentences ("phage X was isolated from hospital wastewater") are routine in the abstracts of on-topic papers, so penalising `sewage` and `wastewater` as off-topic would suppress relevant records. These terms are therefore not in the penalty list, and `scripts/07b_triage_sensitivity.py` records the re-scored audit trail.

**The prioritised set is 282 records.** Of these, 231 yielded machine-parseable full text and 51 were assessed from the abstract; 116 were excluded at eligibility and **166 studies were included**, contributing 442 antibiotic-level observations of which 295 had a determinate outcome.

The complete record-level trail is in `data/triage_rescore_4402.csv` (scores and rationale for all 4,402 records, no article text), `data/screening_audit_4402.csv` and `data/eligibility_exclusions_97.csv`, and every eligibility decision for the sensitivity-screen records is set out individually in `data/priority4_screening_log_28.csv` with the supporting evidence quotation.

## 2. η² is not an ICC — this is the most important entry

A variance-partition statistic is sometimes mistaken for an ICC. The expression

```
Σ nᵢ(pᵢ − p̄)² / [(n − 1) · p̄(1 − p̄)]
```

The numerator is the between-study sum of squares. The denominator is `(n−1) × p̄(1−p̄)`, and since `n·p̄(1−p̄)` **is** the total sum of squares, this expression evaluates to `SS_between / SS_total` — that is, **η²**, the proportion of total variance explained by study membership.

Verified on the data:

```
SS_between = 16.3224
SS_total   = 34.1176
SS_between / SS_total           = 0.4784   <- η²
SS_between / ((n-1) p̄(1-p̄))    = 0.4802   <- the value this expression returns
```

**η² is not an intraclass correlation and is not valid in the design-effect formula.** The manuscript reports a one-way ANOVA ICC(1):

```
ICC(1) = (MS_between − MS_within) / (MS_between + (n₀ − 1)·MS_within) = 0.157
```

Consequences:

| Quantity | If η² were used (invalid) | With ICC(1) (reported) |
|---|---|---|
| Design effect | 1.76 | **1.25** |
| Effective n | 154 | **235** |
| Cluster-adjusted CI | 90.0–88.9 (not usable) | **79.0–88.4** |

**The substantive conclusion is unchanged:** the cluster-adjusted interval remains wider than the unadjusted Wilson interval, so the claim that a naive CI overstates precision still holds. The corrected interval also agrees closely with the cluster bootstrap (78.0–89.3%), which the η²-based interval did not.

`validate_results.py` computes ICC(1) and prints η² on an `[INFO]` line explicitly labelled descriptive-only.

---

## 3. Clustering is detectable but modest

Under complete independence with these fixed cluster sizes, η² has an expected value of `(G−1)/(n−1)`, not zero:

```
(G−1)/(n−1) = 104/271 = 0.384
observed    =                  0.478
```

So the observed value is only ~9.5 percentage points above what chance alone produces. A permutation test (20,000 shuffles of outcome across studies, cluster sizes held fixed) gives **P = 0.020**.

Clustering is real but marginal, and 52 of 113 clusters are singletons that carry no information about the ICC.

---

## 4. Cluster size range is 1–10, not 1–8

The cluster-size distribution is:

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

The outcome-selected figure is 37.1 pp; it is shown only to demonstrate the size of the selection effect.

---

## 6. Counts that depend on the denominators

These counts all change with the number of determinate observations, so they are listed together to keep the manuscript internally consistent:

| Quantity | Value |
|---|---|
| Studies testing ≥ 2 antibiotics | 57 |
| Discordant among those | 19 (33.3%) |
| Head-to-head paired studies | 26 |
| Paired b / c | 6 / 1 |
| McNemar exact P | 0.13 |
| Design-based studies / arms | 55 / 226 |
| Outcome-selected studies | 21 |
| Bactericidal vs static | 83.5% (266) vs 86.2% (29), P = 1.00 |
| β-lactam vs all others | 82.1% (112) vs 84.7% (183), P = 0.63 |
| Envelope antagonism | 8/174 (4.6%) vs 8/121 (6.6%), P = 0.45 |
| Assay mix (checkerboard/biofilm/time-kill/in vivo) | 180 / 46 / 37 / 17 |


---

## 7. The GEE is fitted on the design-based subset, not the paired subset

This is the easiest place in the paper to mis-reproduce, so it is stated explicitly.

The GEE for the fluoroquinolone vs beta-lactam contrast is fitted on the **design-based subset** (studies testing >= 2 distinct mechanisms), restricted to the two mechanism classes of interest — **not** on the 26 paired studies.

```
reported: 55 design-based studies, 117 arms (FQ 37, BL 80)
          OR 3.97 (0.59-26.59), P = 0.156

Fitting on the 26 paired studies instead gives a different point estimate and a
much narrower interval, because that subset is outcome-selected by construction.
```

Switching estimands without saying so would overstate precision.

Specification: Binomial family, logit link, exchangeable working correlation, robust (sandwich) variance, study as the clustering variable. Reproduce with `python scripts/22_gee_reanalysis.py`, which refits the model and asserts the specified result before printing it.

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

This does not change the selected set at the threshold used, but it is a documented divergence between the stated rule and the implementation. `scripts/07b_triage_sensitivity.py` applies case-insensitive deduplication.

---

## 10. "Within-study" and "design-based" are not synonyms

Three subset definitions appear in the paper and are easy to confuse:

| Subset | Definition | n |
|---|---|---|
| **Pooled** | all rows with a determinate outcome | 295 arms / 113 entries |
| **Design-based** | `nunique(Antibiotic_Mechanism) >= 2` — selected on *design*, never outcome | 226 arms / 55 studies |
| **Outcome-selected** | `nunique(Antibiotic_Tested) >= 2` **and** `0 < synergy_count < arm_count` | 21 studies |
| **Head-to-head** | studies measuring **both** PBP and DNA-gyrase mechanisms | 26 studies |

The outcome-selected subset exists only to demonstrate the size of the selection effect and must never be used to estimate comparative class performance.

---

## 11. Organism assignment uses regex, not taxonomy

`_organisms` is populated by regular-expression matching over title and abstract, not by taxonomic lookup. Two consequences are documented:

- A record can match several organisms; the first match is used as the primary organism.
- **One sensitivity-screen record was excluded for exactly this reason.** Phage vB_KquU_φKuK6 is a *Klebsiella quasipneumoniae* / *K. variicola* phage that the regex labelled `K. pneumoniae` via the shared "Klebsiella" token. It was excluded on manual review because non-*pneumoniae* Klebsiella species are not ESKAPE organisms.

Any automated re-run of the screen will re-include that record. This is why the 28-record screening log is committed with explicit per-record decisions rather than left to the regex.

---

## 12. Assay counts are primary-assay, not exclusive

`Assay_Type` lists the primary assay first; 15 arms carry no assay label. The reported counts (180/46/37/17) are **primary-assay** counts and sum to fewer than the total number of arms. True non-exclusive counts are much larger (an arm can be both checkerboard and time-kill). Do not treat these as exclusive categories.

---

## 13. Why the pooled estimate is lower than a naively selected one

This is not a discrepancy but a result worth recording, because it is the paper's central argument in miniature.

Pooled proportions in this literature are inflated by selective reporting of positive combinations, and the direction of that bias is visible within our own dataset. The pooled estimate (83.7%) sits above the design-based estimate (82.7%), and the gap between fluoroquinolones and β-lactams is smaller in the pooled data (9.3 pp) than in the outcome-selected subset (37.1 pp). The design-based analysis moves least under selection, which is the property that justifies using it as the primary comparison.
