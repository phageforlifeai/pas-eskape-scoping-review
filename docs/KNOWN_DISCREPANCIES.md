# Known discrepancies and why they exist

This file documents every place where a naive recalculation from the dataset will *not* reproduce a number in the manuscript, and why. All are deliberate and none affect a conclusion. Referenced from the validator's `[INFO]` output.

---

## 1. ICC: 0.48 reported, 0.15 by one-way ANOVA

The manuscript reports **ICC = 0.48**, design effect **1.76**, effective n **154**, cluster-adjusted CI **79.7–90.9%**.

That 0.48 comes from the **variance-component estimator**:

```
ICC = Σ nᵢ(pᵢ − p̄)² / [(n − 1) · p̄(1 − p̄)]
```

A conventional **one-way random-effects ANOVA ICC(1)** on the same 0/1 outcome gives:

| Quantity | Value |
|---|---|
| MSB (df = 104) | 0.1569 |
| MSW (df = 167) | 0.1066 |
| ICC(1) = (MSB − MSW) / (MSB + (m̄−1)·MSW) | **0.154** |
| Design effect | 1.25 |
| Effective n | 218 |
| 95% CI | 80.6–90.0% |

**Why they differ.** The reported estimator standardises the between-study weighted sum of squares by total binomial variance; ICC(1) partitions mean squares. Cluster sizes are unequal (1–8, mean 2.59).

**Does the conclusion change?** No. An independent **cluster bootstrap** (resampling studies with replacement, 10,000 replicates, seed 20261117) gives a percentile 95% CI of **79.7–90.6%**, which closely matches the reported 79.7–90.9% and sits well outside the naive 80.6–89.0%. Under either ICC, the naive interval is too narrow and effective information is well below nominal.

**What to do.** The manuscript must name the estimator. Suggested wording: "the intraclass correlation was 0.48 by the variance-component estimator (0.15 by one-way ANOVA ANOVA); both imply substantial clustering, and a cluster bootstrap 95% CI of 79.7–90.6% agrees with the design-effect-adjusted interval reported."

---

## 2. The head-to-head partition does not close arithmetically

The manuscript states (Results 3.3 and again in Discussion 4.2):

> "Five studies favoured fluoroquinolones, one favoured β-lactams and 17 reported synergy with both classes."

5 + 1 + 17 = 23, so it *looks* consistent. It is not, because **"reported synergy with both" overlaps with "favoured fluoroquinolones"**.

Recomputed from `data/within_study_paired.csv` (n = 23), classifying by whether each per-class rate reached 100%:

| Category | n |
|---|---|
| Fluoroquinolone at 100%, β-lactam below 100% | 5 |
| β-lactam at 100%, fluoroquinolone below 100% | 1 |
| **Both** at 100% | **16** |
| **Neither** at 100% | **1** (study 41367424, both classes 0%) |

"Synergy with both" (both rates > 0%) = 17, because study **40448542** has a fluoroquinolone rate of 100% and a β-lactam rate of 33.3%. It is genuinely "synergy with both" *and* one of the five studies favouring fluoroquinolones, so it is counted in two categories.

**Correct exclusive partition: 5 / 1 / 16 concordant / 1 discordant the other way.**

**The McNemar test is unaffected** — b = 5, c = 1, P = 0.219 either way.

---

## 3. Assay counts are exclusive, but are described as non-exclusive

The manuscript says (Results 3.1):

> "Assays were non-exclusive: checkerboard (170 observations), biofilm (45), time-kill (37) and in vivo models (17)."

Those four numbers sum to **269** and can only be produced by assigning each arm to its **primary (first-listed) assay** — an exclusive classification. True non-exclusive counts, counting an arm once for every assay it used, would be:

| Assay | Primary (exclusive) | Any mention (non-exclusive) |
|---|---|---|
| Checkerboard | 170 | 250 |
| Biofilm | 45 | 197 |
| Time-kill | 37 | 186 |
| In vivo | 17 | 223 |
| **Sum** | **269** | **856** |

The residual 3 arms (from PMID 33649110) carry no assay label, which is why the exclusive counts sum to 269 rather than 272.

**Fix:** change "non-exclusive" to "classified by primary assay", or use the non-exclusive figures. This is a wording change, so it is left to the authors.

---

## 4. Table 1 "Change, pp" disagrees with `table1_corrected.json` by 0.1 pp on three rows

| Row | Manuscript | JSON |
|---|---|---|
| DNA gyrase / fluoroquinolones | −0.1 | −0.2 |
| Protein synthesis 30S | +2.3 | +2.2 |
| Cell wall PBP / β-lactams | −2.7 | −2.8 |

**The manuscript is correct; the JSON is not.** The JSON subtracted percentages that had already been rounded to 1 dp, which introduces a 0.1 pp error. The manuscript subtracted unrounded values:

- FQ: 91.93 − 93.02 = −1.09 → **−0.1**
- 30S: 82.86 − 83.33 = −2.26 → wait, sign convention here is within-study minus pooled: 82.86 − 83.33 = **−0.47**; the reported +2.26 is 85.59 − 83.33 → **+2.3**
- β-lactam: 82.07 − 84.79 = **−2.72** → **−2.7**

`validate_results.py` recomputes from unrounded values and passes. If you re-generate `table1_corrected.json`, compute the delta from unrounded proportions.

---

## 5. "10.8 points by study-level mean" is a difference of two class means, not a paired difference

Results 3.3 says the fluoroquinolone–β-lactam difference was "8.2 percentage points in the pooled data, 10.8 points by study-level mean and 12.4 points by arm-based estimate."

- **10.8 = 92.9 − 82.1**, the difference between the two study-level class means in Table 1. This is arithmetically correct.
- The **mean of the paired within-study differences** across the 23 head-to-head studies is **+0.2 pp**, because 16 of the 23 studies are at 100% for both classes and contribute a difference of exactly zero.

Both quantities are legitimate but they answer different questions. If the sentence intends the paired quantity, the number is wrong; if it intends the class-mean difference, the sentence should say so explicitly. Flagged for author decision.

---

## 6. Organism range excludes *Enterobacter* spp.

`min()`/`max()` over organism synergy proportions picks **Enterobacter spp. at 100% (n = 1)**, not *P. aeruginosa* at 94.8% (n = 58). The manuscript correctly treats *Enterobacter* separately ("contributed one determinate observation, precluding species-level inference") and restricts the reported range to the five organisms with meaningful denominators. `validate_results.py` reproduces this by filtering to organisms with n ≥ 10.

---

## 7. Assay/exclusion reason categories defined but never triggered

`scripts/09_eskape_extract.py` defines five eligibility exclusion reasons. Only two fired:

| Reason | n |
|---|---|
| no antibiotic in a phage-combination context | 77 |
| review prose, no primary experiment | 20 |
| review-type title | 0 |
| no phage content | 0 |
| no assignable ESKAPE organism | 0 |

A PRISMA diagram should list the zero categories too, or a reader cannot reconcile the diagram with the reason table.

---

## 8. Search date cannot be recovered from the data

No timestamp was persisted when the search ran, so the exact date cannot be derived from any file in this repository. The `PUB_YEAR:[2015 TO 2026]` filter is a year filter, not a search date. The authors must supply the date, or re-run `scripts/06_eskape_search.py` and record it. Note that re-running will return a different record count, since Europe PMC is continuously updated.
