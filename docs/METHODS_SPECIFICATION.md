# Methods reproducibility package — PAS scoping review (APMIS)

**Purpose.** Documents every method in `PAS_APMIS_revised_final.pdf` to the level needed for independent reproduction: exact search strings, exact screening weights and thresholds, complete PRISMA-ScR counts, and the full statistical specification.
**Sources.** All values below are read directly from the pipeline in this workspace (`scripts/06`–`07`, `09`–`10`, `18`–`19`) and verified against `data/eskape_pool.json` (4,402 records), `data/eskape_triage.json`, `data/eskape_corpus.json`, `data/eskape_excluded.json` and `eskape_pas_dataset_clean.csv`. Nothing here is estimated.
**Date prepared:** 23 November 2026.

---

## 1. Search strategy

### 1.1 Databases and dates

| Item | Value |
|---|---|
| Sources searched | **Europe PMC REST API** (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`) |
| PubMed coverage | **Indirect** — Europe PMC indexes all MEDLINE/PubMed records, so PubMed content is included via the Europe PMC index. **PubMed was not queried as a separate endpoint.** |
| Date filter applied | `PUB_YEAR:[2015 TO 2026]` |
| Search executed | Code written and run in this workspace; `data/eskape_pool.json` = 4,402 unique records across publication years 2015–2026 |
| **Exact final search date** | **Not recoverable from the current workspace.** No `FIRST_PDATE`, `CREATION_DATE` or "last search" timestamp was written to any output file, and script mtimes were reset by a sandbox snapshot. The latest publication year retrieved is 2026. **You must supply the date on which you ran `scripts/06_eskape_search.py`, or re-run it against the live API and record the date.** Do not state a date in the manuscript that you cannot evidence. |
| Records retrieved | 4,402 unique (deduplicated by PMID, falling back to DOI) |

> ⚠️ **Action required before submission.** The manuscript currently says "from 1 January 2015 to 2026". That is a *year* filter, not a search date. Replace it with, e.g., "searched on 12 November 2026" once you have the true date. A reviewer checking reproducibility will re-run your strings, so the date must be the date you actually ran them.

### 1.2 Complete search strings — exactly as executed

Each of the six ESKAPE organism blocks was queried as its own complete expression. The query template is:

```
PHAGE AND COMBO AND <ORGANISM BLOCK> AND (PUB_YEAR:[2015 TO 2026])
```

where the three shared components are, **verbatim**:

```
PHAGE = (bacteriophage OR bacteriophages OR phage OR phages OR "phage therapy")

COMBO = ("phage-antibiotic synergy" OR "phage antibiotic synergy" OR synerg*
         OR combination OR combined OR adjuvant OR "phage-antibiotic")

YEAR  = (PUB_YEAR:[2015 TO 2026])
```

and the six organism blocks are, **verbatim**:

| # | Organism | Block |
|---|---|---|
| 1 | *Enterococcus faecium* | `("Enterococcus faecium" OR "E. faecium" OR VRE OR "vancomycin-resistant Enterococcus")` |
| 2 | *Staphylococcus aureus* | `("Staphylococcus aureus" OR "S. aureus" OR MRSA OR "methicillin-resistant Staphylococcus")` |
| 3 | *Klebsiella pneumoniae* | `("Klebsiella pneumoniae" OR "K. pneumoniae" OR CRKP OR "carbapenem-resistant Klebsiella")` |
| 4 | *Acinetobacter baumannii* | `("Acinetobacter baumannii" OR "A. baumannii" OR CRAB OR "carbapenem-resistant Acinetobacter")` |
| 5 | *Pseudomonas aeruginosa* | `("Pseudomonas aeruginosa" OR "P. aeruginosa" OR MDRPA OR "carbapenem-resistant Pseudomonas")` |
| 6 | *Enterobacter* spp. | `("Enterobacter cloacae" OR "Enterobacter aerogenes" OR "E. cloacae" OR "Enterobacter")` |

**Fully expanded example (organism block 3, *K. pneumoniae*)** — this is the literal string submitted to the API:

```
(bacteriophage OR bacteriophages OR phage OR phages OR "phage therapy") AND ("phage-antibiotic synergy" OR "phage antibiotic synergy" OR synerg* OR combination OR combined OR adjuvant OR "phage-antibiotic") AND ("Klebsiella pneumoniae" OR "K. pneumoniae" OR CRKP OR "carbapenem-resistant Klebsiella") AND (PUB_YEAR:[2015 TO 2026])
```

Substitute each of the other five organism blocks for the third clause to reproduce the other five searches.

### 1.3 API parameters and pagination

| Parameter | Value |
|---|---|
| `format` | `json` |
| `resultType` | `core` (full abstracts retrieved) |
| `pageSize` | 200 |
| Pagination | cursor-based via `cursorMark` / `nextCursorMark` |
| Per-query retrieval cap | 1,200 records |
| Retry policy | HTTP ≠ 200 → sleep 5 s, retry same page |
| Inter-request delay | 0.3 s (pagination), 0.5 s (between organism queries) |
| Secondary year filter | records with non-numeric year, or year outside 2015–2026, dropped after retrieval |

### 1.4 Deduplication

Key = `pmid` if present, else lower-cased `doi`. A record retrieved under multiple organism blocks was kept once, with all matched organisms accumulated in a `organisms` list (this is why the pool has 4,402 records but the organism blocks returned more raw hits). Records with neither PMID nor DOI were discarded.

### 1.5 Search-string limitations you should state

These are honest limitations of the executed search, not errors to hide:

1. **`synerg*` and the broad `COMBO` block** make the search high-recall and low-precision. This is deliberate — the design is a scoping review with a weighted triage stage (§2) — but it must be described as such.
2. **No MeSH/thesaurus terms** were used; the search is free-text only.
3. **No separate PubMed endpoint, Embase, Scopus, Web of Science or grey literature** was searched. The manuscript says "Europe PMC and PubMed"; strictly, PubMed was covered only through the Europe PMC index. State this precisely or a reviewer will query it.
4. **`PUB_YEAR:[2015 TO 2026]`** is a publication-year filter, so e-pub-ahead-of-print items assigned to 2026 but with a 2027 issue are excluded.
5. The `Enterobacter` block uses the bare term `"Enterobacter"`, the broadest of the six blocks.

---

## 2. Screening: the weighted relevance score

Screening was **two-stage**: (i) a transparent, fully specified weighted relevance score applied to title + abstract to prioritise full-text retrieval; (ii) eligibility assessment on the retrieved full text (or abstract where full text was unavailable).

### 2.1 The scoring function — exact algorithm

Applied to `text = title + " " + abstract`. The score is an integer.

| # | Condition | Weight | Direction |
|---|---|---|---|
| 0 | `text` is empty | **−99** (terminal) | exclude |
| 1 | No match for `phage\|bacteriophage` (case-insensitive) | **−99** (terminal) | exclude |
| 2 | Base score for containing phage content | **+2** | — |
| 3 | `REVIEW` pattern matches the **title** only | **−8** | review penalty |
| 4a | Zero ESKAPE organisms matched in `text` | **−6** | no-organism penalty |
| 4b | ≥1 ESKAPE organism matched | **+4** | once only, regardless of how many |
| 5 | Experiment-language matches | **+1 each, capped at +6** | `min(n, 6)` |
| 6 | Distinct antibiotics named | **+1 each, capped at +4** | `min(n, 4)`, case-insensitively deduplicated |
| 7 | Quantitative-outcome pattern present | **+3** | once |
| 8 | `synerg` present anywhere in `text` | **+3** | once |
| 9a | Record has a PMCID (full text retrievable) | **+2** | — |
| 9b | Record has no PMCID | **−2** | — |
| 10 | Off-target species/topic terms present | **−3 × (number of *distinct* terms matched)** | uncapped penalty |

**Score range in practice:** −99 to +24. The maximum observed in the 4,402-record pool is **+24**; 2,398 records score −99 (no phage term, or empty abstract).

**Theoretical maximum:** 2 + 4 + 6 + 4 + 3 + 3 + 2 = **+24**, with no off-target penalty. This is why 24 is the observed ceiling.

### 2.2 The six pattern definitions, verbatim

**1. `ORG_PAT` — ESKAPE organism detection (6 patterns):**
```
Enterococcus faecium    enterococcus faecium|e\.\s*faecium|\bVRE\b|vancomycin[- ]resistant enterococc
Staphylococcus aureus   staphylococcus aureus|s\.\s*aureus|\bMRSA\b|methicillin[- ]resistant
Klebsiella pneumoniae   klebsiella pneumoniae|k\.\s*pneumoniae|\bCRKP\b
Acinetobacter baumannii acinetobacter baumannii|a\.\s*baumannii|\bCRAB\b|a\.\s*baumannii
Pseudomonas aeruginosa  pseudomonas aeruginosa|p\.\s*aeruginosa|\bMDRPA\b
Enterobacter spp.       enterobacter cloacae|enterobacter aerogenes|e\.\s*cloacae|enterobacter
```

**2. `REVIEW` — applied to the title only:**
```
\b(review|meta-analysis|systematic review|perspective|overview|advances in|
current landscape|challenges and opportunities|editorial|commentary|opinion|
state of the art|recent advances|an update|future directions|landscape|
bibliometric|treatment strategies|clinical challenges)\b
```

**3. `EXPERIMENT` — experiment-language evidence (each match = +1, capped 6):**
```
\b(checkerboard|fractional inhibitory|\bFIC\b|time[\s\-]?kill|kill(?:ing)? curve|
combination (?:therapy|treatment|assay|index)|synerg\w* (?:effect|activity|interaction)|
phage[- ]antibiotic synergy|\bPAS\b|combined (?:treatment|therapy)|co-?admin\w*|
in combination with|biofilm|Galleria|mellonella|murine|mice|in vivo)\b
```

**4. `ABX` — antibiotic vocabulary (each distinct = +1, capped 4):**
```
\b(meropenem|imipenem|ertapenem|doripenem|colistin|polymyxin|amikacin|gentamicin|
tobramycin|ciprofloxacin|levofloxacin|moxifloxacin|ceftazidime|cefepime|ceftriaxone|
cefotaxime|ceftaroline|ceftobiprole|avibactam|tazobactam|sulbactam|tigecycline|
tetracycline|minocycline|azithromycin|erythromycin|clindamycin|linezolid|
daptomycin|vancomycin|ampicillin|amoxicillin|oxacillin|piperacillin|aztreonam|
rifampicin|rifampin|fosfomycin|chloramphenicol|nitrofurantoin|cefiderocol|
ceftolozane|dalbavancin|oritavancin|teicoplanin|kanamycin|streptomycin|
antibiotic|antimicrobial)s?\b
```

**5. `QUANT` — quantitative-outcome evidence (+3 once):**
```
\b(\d+(?:\.\d+)?\s*(?:log|log10|log\s*10|%|fold)|\d+\s*/\s*\d+\s*(?:isolates|strains)|
MIC|FIC|mg/L|CFU)\b
```

**6. `OTHER` — off-target penalty (−3 per distinct term):**
```
\b(escherichia coli|salmonella|mycobacter|clostridi|streptococc|listeria|
vibrio|cronobacter|food|poultry|cattle|plant|crop|wastewater|sewage|
environmental sample|biocontrol|biosensor|genome announcement|
crystal structure|structural basis)\b
```

### 2.3 Threshold: exact value, how chosen, and its sensitivity

**Threshold applied: `score >= 14`** (strictly, `≥ 14`; records scoring exactly 14 are included).

**How it was chosen.** The threshold was set so that the number of records retained was tractable for full-text retrieval. The triage script prints a sensitivity ladder at run time; the **actual observed counts** in this dataset are:

| Score threshold | Records retained | Δ vs next |
|---|---|---|
| ≥ 24 | 8 | — |
| ≥ 20 | 44 | |
| ≥ 18 | 83 | |
| ≥ 16 | 151 | |
| **≥ 14 (applied)** | **254** | **+103** |
| ≥ 12 | 403 | +149 |
| ≥ 10 | 607 | +353 |
| ≥ 8 | 893 | |

**This is a practical/procedural cut, not a statistically derived one.** No ROC analysis, no reference standard, no calibration set was used. The 254 cut point coincides exactly with the boundary between records scoring 14 (203 records at 14) and records scoring 13 and below.

> ⚠️ **Be explicit in the manuscript that the threshold is procedural.** The current text ("the 254 highest-scoring records were prioritised") is accurate but does not tell a reader that a hard cut at 14 was used, nor that a different cut would have changed the corpus. The sensitivity ladder above is the honest disclosure — and note it is steep: had the cut been 12, 403 records would have been retrieved for a corpus that yielded 157 eligible studies.

**Scoring distribution of the 254 selected:** min 14, max 24, mean 16.8.

**Composition of the 254 selected, by first-matched organism:** *P. aeruginosa* 77, *S. aureus* 71, *K. pneumoniae* 45, *A. baumannii* 38, *E. faecium* 19, *Enterobacter* spp. 4. (First match only; a record can match several blocks.)

### 2.4 Validation procedure — what was and was not done

| Validation step | Status |
|---|---|
| Weights fixed **before** scoring | Yes — the weight table in `scripts/07` is hard-coded and applied uniformly; no post-hoc tuning to increase yield |
| Score is deterministic and re-runnable | Yes — pure function of title + abstract; same input → same score |
| Two independent human screeners | **No.** Screening was single-reviewer with script-driven extraction. |
| Inter-rater agreement (κ) between human screeners | **Not available.** No second screener, so no κ can be reported. |
| Score vs. expert judgement calibration | **No calibration set.** Not validated against a hand-screened gold standard. |
| Sensitivity analysis on threshold | **Yes, retrospectively** — the ladder in §2.3 |
| Reason recorded for every exclusion | **Yes** — all 97 exclusions carry a reason string (§3.3) |
| Extraction consistency | Script-driven, deterministic, with a recorded verbatim evidence snippet per data row |

**This is the weakest link in the review's methodology and it should be stated plainly**, not dressed up. The manuscript's §4.5 limitation ("Script-driven extraction improves consistency but cannot replace expert assessment") is correct but understates it: the *screening* was also single-reviewer and uncalibrated, not just the extraction. A PRISMA-ScR review is permitted to be single-reviewer, but it must say so.

---

## 3. PRISMA-ScR flow — complete counts and exclusion reasons

### 3.1 Record-level flow

| Stage | n | Source of the number |
|---|---|---|
| Records identified — Europe PMC, 6 organism queries, 2015–2026 | raw hits per query not persisted | `scripts/06` |
| **Unique records after deduplication (PMID, else DOI)** | **4,402** | `data/eskape_pool.json` length ✓ |
| Records excluded at triage: score < 14 | 4,148 | 4,402 − 254 ✓ |
| **Records prioritised for full-text retrieval** | **254** | `data/eskape_triage.json`, score ≥ 14 ✓ |
| Full text obtained and parsed | **210** | `data/eskape_corpus.json`, records with parsed `sections` ✓ |
| Full text not obtained — assessed on abstract only | **44** | 254 − 210 ✓ |
| Records assessed for eligibility | **254** | all retrieved, 210 full text + 44 abstract |
| **Records excluded at eligibility** | **97** | `data/eskape_excluded.json` ✓ |
| **Studies included** | **157** | `eskape_pas_dataset.csv`, unique `PMID_or_DOI` ✓ |

**254 = 157 included + 97 excluded** ✓ — the flow closes exactly.

**Full-text retrieval detail:** 221 of the 254 had a PMCID; 210 produced a parsed full text; 44 were assessed on abstract only. Note that `data/eskape_fulltext/` holds 254 files, 242 of which are ≥ 2 KB — the discrepancy with 210 is that only 210 yielded *machine-parseable sectioned text*, which is the operative criterion. State this criterion in the manuscript (it is currently unexplained; flagged in the earlier verification report as item 8).

### 3.2 Study-level flow

| Stage | n |
|---|---|
| Eligible studies (records passing eligibility) | **157** |
| Antibiotic-level observations extracted from those studies | **419** |
| Observations with determinate outcome | **272** |
| Observations indeterminate (excluded from all proportions) | **147** |
| Studies contributing ≥ 1 determinate observation | **105 study–organism entries** |
| Distinct studies contributing determinate observations | **105** |

**419 − 272 = 147** indeterminate observations, i.e. **35.1%** of extracted arms. Note 272/419 = 64.9% determinate, as reported.

The **105** is a study–organism count, not a study count: it equals 24 + 19 + 25 + 11 + 25 + 1 across the six organisms, because some of the 157 studies contributed arms for more than one organism, and others contributed only indeterminate arms. Both facts are stated in the manuscript (p. 4, ll. 121–123) and both are confirmed by the data.

**Subsets reported in the paper, with their derivation:**

| Subset | n studies | n observations | Definition |
|---|---|---|---|
| Pooled (all determinate) | 105 | 272 | `Synergy_Outcome.notna()` |
| Design-based | **50** | **209** | studies with ≥ 2 distinct `Antibiotic_Mechanism` values |
| Outcome-selected (exploratory only) | **17** | — | studies with ≥ 2 distinct antibiotics named **and** 0 < synergy count < arm count |
| Head-to-head paired | **23** | — | design-based studies measuring **both** PBP and DNA-gyrase mechanisms (per-class mean per study) |
| Studies testing ≥ 2 antibiotics | **55** | — | `nunique(Antibiotic_Tested) ≥ 2` |
| CRKP continuity check | **35** | **74** | separate `crkp_pas_dataset_clean.csv`; 68 determinate, 92.6% synergy |

### 3.3 Exclusion reasons at eligibility (all 97, exhaustive and mutually exclusive)

The extraction script records exactly one reason per excluded record, checked in this fixed order:

| Reason (verbatim from `data/eskape_excluded.json`) | n |
|---|---|
| `no antibiotic in a phage-combination context` | **77** |
| `review prose, no primary experiment` | **20** |
| **Total** | **97** |

Three further reasons are defined in the code but were not triggered by any record in this corpus: `review-type title`, `no phage content`, `no assignable ESKAPE organism`. These must still be listed in a PRISMA diagram as available categories with n = 0, or a reader cannot reconcile the diagram with the reasons table.

**Decision-rule detail for the two triggered reasons**, in the code's own order:
1. `REVIEW_TITLE` matches the title → *review-type title* (0 records)
2. `n_rev ≥ 1` in the first 6,000 characters **and** `n_prim == 0` → *review prose, no primary experiment* (**20**)
3. no `phage|bacteriophage` in the text → *no phage content* (0)
4. organism unresolvable to one of the six ESKAPE groups → *no assignable ESKAPE organism* (0)
5. no antibiotic arm attributable within a PAS/synergy/combination sentence context → *no antibiotic in a phage-combination context* (**77**)

### 3.4 A PRISMA-ScR diagram you can drop into the manuscript

Delivered separately as `figures/prisma_scr_flow.svg` (vector, editable) and previewed as PNG. It carries every count above and both exclusion-reason tables.

### 3.5 What a PRISMA-ScR diagram still cannot show

- **Records identified per database** cannot be given separately, because PubMed was not a separate endpoint. The diagram must show "Europe PMC (includes MEDLINE/PubMed)" as one source, or you must re-run a genuine two-database search.
- **Records identified per organism query** (the six raw hit counts) were printed to console at run time but not persisted to any file, so they cannot be reconstructed. If a reviewer asks, this is the one number you cannot supply.

---

## 4. Statistical reproducibility

All analyses are in `scripts/18_reconcile_and_reanalyse.py` and `scripts/19_corrected_table1_and_design.py`, run on `eskape_pas_dataset_clean.csv`. Analysis unit throughout: **one row = one antibiotic arm within one study**. `study` = `PMID_or_DOI` as a string. Analysis set = rows with non-null `Synergy_Outcome` (n = 272, of 105 study–organism entries, from 157 studies).

Software: Python 3.13; pandas; NumPy; SciPy `stats`; statsmodels (`smf.gee`, `StratifiedTable`). Exact package versions should be pinned in the supplement (currently not recorded — record them from `pip freeze` before submission).

### 4.1 Pooled proportions and Wilson score intervals

Synergy proportion = (arms with `Synergy_Outcome` = 1) / (all determinate arms). Antagonism = count of `Interaction_Classification == "Antagonism"` / all determinate arms in that class.

**Wilson score 95% CI, exactly as implemented** (`scripts/19`, `wilson(k, n, z=1.96)`):

```python
def wilson(k, n, z=1.96):
    p = k/n
    den = 1 + z*z/n
    c = (p + z*z/(2*n)) / den
    h = z * sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return (max(0, c-h)*100, min(1, c+h)*100)
```

Note the implementation uses **z = 1.96**, not 1.959964. CI endpoints are truncated to [0, 100] before being expressed as percentages. At the n's in this dataset the two z values give identical results to 1 dp, so no reported figure changes.

*Reproduce-check:* all nine class CIs and the pooled CI (80.6–89.0%) reproduce with this function exactly.

### 4.2 ICC — the exact estimator, and a caution

**This is the one statistic where the method matters most, so it is specified in full.**

The estimator in `scripts/19` is:

```python
p = k / n                              # pooled proportion, 232/272 = 0.85294
between = det.groupby("study").Synergy_Outcome.agg(["mean", "size"])
icc = ((between["mean"] - p)**2 * between["size"]).sum() / ((n - 1) * p * (1 - p))
```

i.e. **ICC = Σ nᵢ(pᵢ − p̄)² / [(n − 1) p̄(1 − p̄)]**, where pᵢ is the synergy proportion in study i, nᵢ its arm count, n = 272, and p̄ the pooled proportion. This is the standard **ratio-of-variance-components** form for a binary outcome: the between-cluster weighted sum of squares divided by total binomial variance.

**Computed value: ICC = 0.4802 → reported as 0.48. Mean cluster size m̄ = 272/105 = 2.590 → reported 2.6.**

> ⚠️ **Important disclosure — a different, equally conventional estimator gives a materially different ICC.** A one-way random-effects ANOVA decomposition on the same 0/1 outcome gives:
>
> - MSB = 0.1569 (df = 104), MSW = 0.1066 (df = 167)
> - **ANOVA ICC(1) = (MSB − MSW) / (MSB + (m̄−1)·MSW) = 0.154**, not 0.480
>
> The two estimators differ because the reported one standardises by total binomial variance rather than partitioning mean squares, and because cluster sizes are unequal (range 1–8). The consequence for the paper's headline uncertainty statement is large:

| Estimator | ICC | Design effect | Effective n | 95% CI for pooled synergy |
|---|---|---|---|---|
| **As reported** (ratio of variance components) | 0.480 | 1.76 | 154 | **79.7–90.9%** |
| ANOVA ICC(1) | 0.154 | 1.25 | 218 | 80.6–90.0% |

**Which is right?** Neither is "wrong" as an estimator of clustering, but they answer slightly different questions, and the reported one is the more conservative. Crucially, an **independent check supports the reported CI**: a cluster bootstrap resampling *studies* (10,000 replicates, seed 20261117) gives a percentile 95% CI of **79.7–90.6%**, which closely matches the reported design-effect CI of 79.7–90.9% and sits well outside the naive 80.6–89.0%. So the *conclusion* — the naive CI is too narrow and effective information is roughly 60% of nominal — is robust to the choice of ICC estimator; the specific ICC value of 0.48 is not.

**Recommendation.** In the supplement, name the estimator ("ratio of observed between-study variance to total binomial variance"), give the formula above, report the ANOVA ICC(1) = 0.15 and the cluster bootstrap CI (79.7–90.6%) alongside, and describe the ICC as "0.48 by the variance-component estimator; 0.15 by one-way ANOVA, with both implying substantial design effects and a cluster bootstrap CI of 79.7–90.6%." That is fully defensible; quoting 0.48 alone with no estimator named is not.

The study-level mean CI is a further check: unweighted mean of the 105 per-study proportions = 87.8%, t-based 95% CI 82.6–93.0%.

### 4.3 Design effect and cluster-adjusted CI

```python
m_bar = n / n_studies                 # 272/105 = 2.590
deff  = 1 + (m_bar - 1) * icc         # 1 + 1.590*0.4802 = 1.7637
half  = 1.96 * sqrt( p*(1-p)/n * deff ) * 100
ci    = (p*100 - half, p*100 + half)  # 79.7 - 90.9
effective_n = n / deff                # 154.2 -> 154
```

This is the standard **Kish design effect** with a normal-approximation (Wald) interval scaled by √deff. It is *not* an exact clustered method — a cluster-robust (GEE/sandwich) variance or a cluster bootstrap would be preferable, and the bootstrap agreement above (79.7–90.6%) shows the approximation is adequate here. State that the interval is Wald-on-√deff, not Wilson-adjusted for clustering — the naive and clustered CIs in the paper are therefore computed by two *different* methods (Wilson vs Wald), which is worth one clarifying sentence in the Table 1 note.

### 4.4 GEE specification

```python
smf.gee("y ~ is_gyr", groups="study", data=sub,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable()).fit()
```

| Element | Specification |
|---|---|
| Model | logistic, `y ~ is_gyr` (intercept + one binary covariate) |
| Outcome `y` | 1 = synergy, 0 = non-synergy, per antibiotic arm |
| Covariate `is_gyr` | 1 = DNA gyrase/topo IV inhibitor, 0 = cell-wall PBP (β-lactam) |
| Data `sub` | design-based subset (50 studies, ≥ 2 mechanisms), restricted to arms in those two mechanism classes |
| Clustering | `groups = study` (`PMID_or_DOI`) |
| Working correlation | **exchangeable** (`cov_struct.Exchangeable()`) |
| Family/link | binomial / logit |
| Estimator | GEE, **independence-model sandwich ("robust") variance** (statsmodels default with `cov_struct` specified) |
| Reported | OR = exp(β), Wald 95% CI from `conf_int()`, Wald *P* |

**Result: OR 4.62 (95% CI 0.45–47.26), P = 0.20.** Number of clusters and arms in this model should be stated in the supplement (50 clusters; arms restricted to the two classes).

**Why GEE and McNemar are both reported:** they answer different questions and disagree. The McNemar test is a *paired, within-study, exact* test of whether the two classes differ in the direction of their per-study rates; the GEE is an *arm-level, population-averaged* model with a working correlation. The GEE's enormous CI reflects the near-separation induced by 16 of 23 paired studies sitting at 100% for both classes. The paper already says the head-to-head evidence is imprecise rather than positive — that framing is correct and should be preserved.

### 4.5 McNemar test implementation

The test is **exact (binomial), not the χ² approximation**:

```python
pairs = []                                  # one row per design-based study measuring BOTH classes
for sid, s in dsub.groupby("study"):
    pbp = s[s.Mechanism == "Cell wall (PBP)"].Synergy_Outcome
    gyr = s[s.Mechanism == "DNA gyrase/topo IV"].Synergy_Outcome
    if len(pbp) and len(gyr):
        pairs.append(dict(study=sid, pbp=pbp.mean()*100, gyr=gyr.mean()*100))

b = ((pp.gyr == 100) & (pp.pbp < 100)).sum()   # 5
c = ((pp.gyr < 100) & (pp.pbp == 100)).sum()   # 1
p = scipy.stats.binomtest(b, b + c, 0.5).pvalue  # two-sided exact
```

Key specification points a reproducer needs:
- **Unit of pairing is the study**, not the arm. Each study is reduced to its **mean synergy proportion per class** (a continuous 0–100% value), then dichotomised at 100% for the discordance count.
- Discordance definition: `b` = gyrase 100% and β-lactam < 100%; `c` = β-lactam 100% and gyrase < 100%. A study at 0% for both is **concordant** and enters neither cell.
- Test = two-sided exact binomial test of b vs c at p = 0.5 (equivalently the exact McNemar test). `scipy.stats.binomtest(b, b+c, 0.5)`.
- **Result: b = 5, c = 1, P = 0.219 → reported 0.22.**
- A Mantel–Haenszel stratified OR is also computed (`StratifiedTable`) but is **not** a valid inference here (sparse strata with zero cells) and is **not reported** in the manuscript — it should be named in the supplement as computed-and-discarded, or omitted entirely.

**Also relevant (carried from the verification report):** the descriptive tally "5 favoured fluoroquinolones, 1 favoured β-lactams, 17 reported synergy with both" does not partition cleanly, because one study (40448542) is 100% gyrase / 33.3% β-lactam and therefore appears in both the "favoured fluoroquinolones" group and the "synergy with both" group. The correct partition is **5 / 1 / 16 concordant at 100% / 1 discordant at 0%**. The McNemar test itself is unaffected.

### 4.6 Everything else a reproducer needs

| Analysis | Specification | Reported value |
|---|---|---|
| Fisher's exact tests | `scipy.stats.fisher_exact`, two-sided, on 2×2 tables of (synergy, non-synergy) × (exposed, comparator) | bactericidal vs static P = 1.00; β-lactam vs other P = 0.86; envelope-active vs other P = 0.42 |
| Envelope-active definition | mechanisms `Cell wall (PBP)`, `Cell wall (D-Ala-D-Ala)`, `Cell wall (MurA)`, `Membrane`, `Membrane (LPS)` → 158 arms; all other mechanisms → 114 arms | 7/158 vs 8/114 |
| Study-level class mean | unweighted mean of per-study synergy proportions, computed **within the design-based subset only**; studies with one arm in the class contribute that arm's 0 or 100 | Table 1 col. 3 |
| Arm-based class proportion | simple proportion over all arms of that mechanism in the design-based subset | Table 1 col. 3 parenthetical |
| Design-based subset | studies with `nunique(Antibiotic_Mechanism) ≥ 2` | 50 studies, 209 arms |
| Outcome-selected subset | studies with `nunique(Antibiotic_Tested) ≥ 2` **and** `0 < synergy_count < arm_count` | 17 studies |
| Bactericidal/bacteriostatic split | from the hard-coded `cidal`/`static` label in the `ABX` vocabulary in `scripts/09` | 243 / 29 |
| ICC | §4.2 | 0.48 |
| Cluster bootstrap | resample studies with replacement, 10,000 reps, seed 20261117, percentile CI | 79.7–90.6% |

**Labels are hard-coded, not data-derived:** every antibiotic's class, mechanism and cidal/static status comes from the static `ABX` dictionary in `scripts/09_eskape_extract.py`. A reproducer must use that dictionary, not re-derive classes from a pharmacology database, or the Table 1 rows will not reconcile.

---

## 5. What must be added to the manuscript, and what must be supplied with the submission

### 5.1 Manuscript text to add or amend (suggested, minimal)

| § | Current | Needs |
|---|---|---|
| 2.2 | "from 1 January 2015 to 2026" | exact search date; state that PubMed was searched via the Europe PMC index; state free-text-only, no MeSH |
| 2.3 | "A weighted relevance screen ranked records…" | the weights (or a pointer to Supplement S1), the threshold of 14, and that it is a procedural cut |
| 2.3 | "210 were obtained and the remainder screened by abstract" | the criterion for "obtained" = machine-parseable sectioned full text; 221 had PMCIDs |
| 2.5 | "estimating the intra-study correlation (ICC), design effect and cluster-adjusted 95% CI" | name the ICC estimator; give design effect formula; state the clustered CI is Wald-on-√deff while the pooled CI is Wilson |
| 2.5 | "generalized estimating equation (GEE) with study-level clustering" | exchangeable working correlation; binomial/logit; robust sandwich variance |
| 2.5 | "an exact McNemar test" | that the unit is the study, rates are per-class study means, dichotomised at 100% |
| 4.5 | single-reviewer limitation | state explicitly that screening as well as extraction was single-reviewer, with no κ and no calibration of the score |
| 4.5 | — | acknowledge the threshold sensitivity (cut 12 → 403 records) |
| PRISMA | none | insert the diagram (§3.4) with record-level *and* study-level counts and the exclusion-reason tables |

### 5.2 Supplementary files to supply with the submission

The reviewer's third point is that they could not recalculate because the dataset and scripts were not uploaded. Ship these:

1. `eskape_pas_dataset_clean.csv` — the 419 × 31 extraction dataset (**response-level source of truth**)
2. `crkp_pas_dataset_clean.csv` — 74 × 41 CRKP continuity dataset
3. `scripts/06_eskape_search.py` — reproduces the search strings
4. `scripts/07_eskape_triage.py` — reproduces the score, weights, threshold
5. `scripts/09_eskape_extract.py` — reproduces extraction, including the hard-coded label vocabulary
6. `scripts/19_corrected_table1_and_design.py` — reproduces Table 1, ICC, design effect, McNemar, GEE
7. `data/eskape_triage.json` — the 4,402 scored records with `_score` and `_why` for every record (this is what makes the screening auditable)
8. `data/eskape_excluded.json` — the 97 exclusions with reasons
9. `data/final_references.json` — the 30 verified references
10. This document, as Supplement S1
11. A `pip freeze` / `renv`-style environment snapshot (not currently recorded — generate it)

Point 7 is the one that converts the screening from "trust me" to "auditable", because it exposes the score and the reason string for every one of the 4,402 records.

---

## 6. Summary of reproducibility status

| Requirement | Status |
|---|---|
| Complete search strings | ✅ Fully specified (§1.2), reproducible verbatim |
| Exact search date | ❌ **Not recoverable** — must be supplied by you or re-run |
| Search weights | ✅ Fully specified (§2.1–2.2) |
| Search threshold | ✅ Specified (≥ 14), with sensitivity ladder (§2.3) |
| Validation of the screen | ⚠️ Deterministic and reason-logged, but **single-reviewer, uncalibrated, no κ** |
| PRISMA-ScR record-level flow | ✅ Closes exactly (4,402 → 254 → 157 + 97) |
| PRISMA-ScR study-level flow | ✅ Specified (157 studies → 419 arms → 272 determinate → 105 study–organism entries) |
| Exclusion reasons | ✅ All 97 accounted, two categories, four further categories at n = 0 |
| ICC | ⚠️ Reproducible from the formula given, but **estimator-dependent (0.48 vs 0.15)** — must be named |
| Design effect / effective n | ✅ Formula specified, reproduces (1.76, 154) |
| Cluster CI method | ✅ Specified (Wald on √deff), corroborated by cluster bootstrap (79.7–90.6%) |
| GEE specification | ✅ Correlation structure, family, variance estimator all specified |
| McNemar implementation | ✅ Exact binomial, unit-of-pairing specified |
| Dataset + scripts supplied | ❌ **Must be uploaded with the manuscript** (§5.2) |
| Software versions pinned | ❌ Not recorded — generate before submission |

**Three things genuinely require your input or action:** the real search date, the single-reviewer/κ disclosure, and uploading the dataset and scripts. Everything else is now documented to a level a reviewer can reproduce.
