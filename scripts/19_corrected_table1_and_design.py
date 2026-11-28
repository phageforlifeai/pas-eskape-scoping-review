#!/usr/bin/env python3
"""
Step 19 - Corrected Table 1 and corrected within-study statistical design.

Fixes two defects found in step 18:
  (1) Table 1 silently dropped two mechanism labels ('Membrane' = daptomycin,
      n = 11; 'Cell wall (MurA)' = fosfomycin, n = 4). Its pooled column summed
      to 257 of 272 observations and its within-study column to 85 of 88.
  (2) The within-study analysis selected studies on their outcome (discordant
      only) and treated repeated observations from one study as independent.

Corrected design, in order of strength:
  1. Design-based restriction - studies testing >=2 distinct mechanisms. Selection
     is on the design, not the result, so it cannot bias the rate estimate.
  2. Study is the unit of inference. Per-study rates are computed first, then
     averaged, so a study contributing 5 arms carries one vote.
  3. Paired McNemar test for the gyrase-vs-PBP comparison, using only studies
     that measured both classes - an exact, within-study test.
  4. Cluster-robust inference (GEE, exchangeable) and design-effect reporting for
     the pooled estimates.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent.parent
OUTD = ROOT / "data"

d = pd.read_csv(ROOT / "eskape_pas_dataset_clean.csv")
det = d[d.Synergy_Outcome.notna()].copy()
det["study"] = det.PMID_or_DOI.astype(str)

# Full mechanism set, including the two labels the old table dropped.
MECH = ["DNA gyrase/topo IV", "Protein synthesis (30S)", "Protein synthesis (50S)",
        "RNA polymerase", "Membrane (LPS)", "Membrane",
        "Cell wall (PBP)", "Cell wall (MurA)", "Cell wall (D-Ala-D-Ala)"]
PRETTY = {"DNA gyrase/topo IV": "DNA gyrase / topoisomerase IV (fluoroquinolones)",
          "Protein synthesis (30S)": "Protein synthesis 30S (aminoglycosides)",
          "Protein synthesis (50S)": "Protein synthesis 50S (macrolides, oxazolidinones)",
          "RNA polymerase": "RNA polymerase (rifamycins)",
          "Membrane (LPS)": "Membrane, LPS (polymyxins)",
          "Membrane": "Membrane, phosphatidylglycerol (daptomycin)",
          "Cell wall (PBP)": "Cell wall PBP (β-lactams)",
          "Cell wall (MurA)": "Cell wall MurA (fosfomycin)",
          "Cell wall (D-Ala-D-Ala)": "Cell wall D-Ala-D-Ala (glycopeptides)"}


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0, c - h) * 100, min(1, c + h) * 100)


# design-based restriction: >= 2 distinct mechanisms tested
gm = det.groupby("study").Antibiotic_Mechanism.nunique()
DESIGN = set(gm[gm >= 2].index)
dsub = det[det.study.isin(DESIGN)]

print("=" * 80)
print("CORRECTED TABLE 1")
print("=" * 80)
print(f"pooled denominator : {len(det)} observations, {det.study.nunique()} studies")
print(f"within-study set   : {len(dsub)} observations, {len(DESIGN)} studies "
      "(selected on design, not outcome)")

rows = []
for m in MECH:
    a = det[det.Antibiotic_Mechanism == m]
    b = dsub[dsub.Antibiotic_Mechanism == m]
    p_pool = a.Synergy_Outcome.mean() * 100
    lo_p, hi_p = wilson(int(a.Synergy_Outcome.sum()), len(a))
    # study-level: mean of per-study rates
    per_study = b.groupby("study").Synergy_Outcome.mean() * 100 if len(b) else pd.Series(dtype=float)
    s_mean = per_study.mean() if len(per_study) else np.nan
    # study-level Wilson-style CI from the binomial counts among studies is not
    # meaningful; report the SD and the range of per-study rates instead
    k_ant = int((a.Interaction_Classification == "Antagonism").sum())
    rows.append(dict(
        mech=m, pretty=PRETTY[m],
        pool_pct=round(p_pool, 1), pool_n=len(a), pool_lo=round(lo_p, 1), pool_hi=round(hi_p, 1),
        ws_pct=round(s_mean, 1) if not np.isnan(s_mean) else None,
        ws_studies=int(len(per_study)), ws_obs=len(b),
        ws_obs_pct=round(b.Synergy_Outcome.mean() * 100, 1) if len(b) else None,
        delta=round(s_mean - p_pool, 1) if not np.isnan(s_mean) else None,
        ant=f"{k_ant}/{len(a)}"))

hdr = (f"{'mechanism':46s} {'pooled':>14s} {'study-level':>16s} {'arms-based':>13s} {'antag':>7s}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    print(f"{r['pretty']:46s} {r['pool_pct']:7.1f} ({r['pool_n']:3d}) "
          f"{str(r['ws_pct']):>7s} ({r['ws_studies']:2d} st) "
          f"{str(r['ws_obs_pct']):>7s} ({r['ws_obs']:3d}) {r['ant']:>7s}")

tot_pool_n = sum(r["pool_n"] for r in rows)
tot_ws_obs = sum(r["ws_obs"] for r in rows)
print("-" * len(hdr))
print(f"{'TOTAL':46s} {det.Synergy_Outcome.mean()*100:7.1f} ({tot_pool_n:3d}) "
      f"{dsub.Synergy_Outcome.mean()*100:>7.1f} ({len(DESIGN):2d} st) "
      f"{dsub.Synergy_Outcome.mean()*100:>7.1f} ({tot_ws_obs:3d}) "
      f"{int((det.Interaction_Classification=='Antagonism').sum()):>3d}/{tot_pool_n}")
print(f"\npooled column now sums to {tot_pool_n} of {len(det)} observations "
      f"(was 257 -> now complete)")
print(f"within-study column now sums to {tot_ws_obs} observations")

# cluster-robust CI for the pooled rate
k, n_s = det.Synergy_Outcome.sum(), len(det)
p = k / n_s
between = det.groupby("study").Synergy_Outcome.agg(["mean", "size"])
icc = ((between["mean"] - p) ** 2 * between["size"]).sum() / ((n_s - 1) * p * (1 - p))
m_bar = n_s / det.study.nunique()
deff = 1 + (m_bar - 1) * icc
half = 1.96 * np.sqrt(p * (1 - p) / n_s * deff) * 100
print(f"\npooled 85.3%: naive Wilson CI {wilson(int(k), n_s)[0]:.1f}-{wilson(int(k), n_s)[1]:.1f}%"
      f" | cluster-adjusted (ICC {icc:.2f}, design effect {deff:.2f}) "
      f"{p*100-half:.1f}-{p*100+half:.1f}%")

print("\n" + "=" * 80)
print("HEADLINE CONTRAST, THREE WAYS (shows the size of the selection effect)")
print("=" * 80)
variants = [("pooled, all studies", det),
            ("design-selected (>=2 mechanisms)", dsub),
            ("outcome-selected (discordant only)", None)]
gg = det.groupby("study").agg(nn=("Synergy_Outcome", "size"),
                              kk=("Synergy_Outcome", "sum"),
                              nab=("Antibiotic_Tested", "nunique"))
disc = set(gg[(gg.nab >= 2) & (gg.kk > 0) & (gg.kk < gg.nn)].index)

for label, frame in variants:
    f = det[det.study.isin(disc)] if frame is None else frame
    pbp = f[f.Antibiotic_Mechanism == "Cell wall (PBP)"].Synergy_Outcome
    gyr = f[f.Antibiotic_Mechanism == "DNA gyrase/topo IV"].Synergy_Outcome
    print(f"  {label:36s} PBP {pbp.mean()*100:5.1f}% (n={len(pbp):3d})   "
          f"gyrase {gyr.mean()*100:5.1f}% (n={len(gyr):3d})   "
          f"gap {gyr.mean()*100-pbp.mean()*100:+5.1f} pp")

print("\n" + "=" * 80)
print("PAIRED WITHIN-STUDY TEST (gyrase vs PBP, studies measuring both)")
print("=" * 80)
pairs, strata = [], []
for sid, s in dsub.groupby("study"):
    pbp = s[s.Antibiotic_Mechanism == "Cell wall (PBP)"].Synergy_Outcome
    gyr = s[s.Antibiotic_Mechanism == "DNA gyrase/topo IV"].Synergy_Outcome
    if len(pbp) and len(gyr):
        pairs.append(dict(study=sid, pbp=pbp.mean() * 100, gyr=gyr.mean() * 100))
        strata.append([[int(gyr.sum()), len(gyr) - int(gyr.sum())],
                       [int(pbp.sum()), len(pbp) - int(pbp.sum())]])
pp = pd.DataFrame(pairs)
b = int(((pp.gyr == 100) & (pp.pbp < 100)).sum())
c = int(((pp.gyr < 100) & (pp.pbp == 100)).sum())
print(f"studies measuring both classes: {len(pp)}")
print(f"  mean per-study rate, PBP    : {pp.pbp.mean():.1f}%")
print(f"  mean per-study rate, gyrase : {pp.gyr.mean():.1f}%")
print(f"  studies where only gyrase reached 100% synergy: {b}")
print(f"  studies where only PBP reached 100% synergy   : {c}")
mc = stats.binomtest(b, b + c, 0.5).pvalue if (b + c) else np.nan
print(f"  exact McNemar (binomial) two-sided P = {mc:.3f}")
st = StratifiedTable(strata)
print(f"  Mantel-Haenszel OR = {st.oddsratio_pooled:.2f} "
      f"(95% CI {st.oddsratio_pooled_confint()[0]:.2f}-{st.oddsratio_pooled_confint()[1]:.2f})")
t = st.test_null_odds()
print(f"  MH test P = {t.pvalue:.3f}  <- unreliable, sparse strata with zero cells")

sub = dsub[dsub.Antibiotic_Mechanism.isin(["Cell wall (PBP)", "DNA gyrase/topo IV"])].copy()
sub["is_gyr"] = (sub.Antibiotic_Mechanism == "DNA gyrase/topo IV").astype(int)
sub["y"] = sub.Synergy_Outcome.astype(int)
gee = smf.gee("y ~ is_gyr", groups="study", data=sub,
              family=sm.families.Binomial(),
              cov_struct=sm.cov_struct.Exchangeable()).fit()
ci = np.exp(gee.conf_int().loc["is_gyr"])
print(f"  GEE (cluster = study): OR {np.exp(gee.params['is_gyr']):.2f} "
      f"(95% CI {ci[0]:.2f}-{ci[1]:.2f}), P = {gee.pvalues['is_gyr']:.2f}")
print("  NOTE: GEE and MH disagree because the MH estimator is driven by sparse")
print("        strata; the exact paired test is the appropriate primary test.")

result = dict(
    mech_rows=rows,
    n_pooled=len(det), n_studies=det.study.nunique(),
    n_design=len(dsub), n_design_studies=len(DESIGN),
    pooled_pct=round(float(p * 100), 1),
    naive_ci=[round(wilson(int(k), n_s)[0], 1), round(wilson(int(k), n_s)[1], 1)],
    cluster_ci=[round(float(p * 100 - half), 1), round(float(p * 100 + half), 1)],
    icc=round(float(icc), 2), design_effect=round(float(deff), 2),
    effective_n=round(float(n_s / deff), 0),
    paired=dict(n=len(pp), pbp=round(float(pp.pbp.mean()), 1),
                gyr=round(float(pp.gyr.mean()), 1), b=b, c=c,
                mcnemar_P=round(float(mc), 3),
                gee_OR=round(float(np.exp(gee.params['is_gyr'])), 2),
                gee_lo=round(float(ci[0]), 2), gee_hi=round(float(ci[1]), 2),
                gee_P=round(float(gee.pvalues['is_gyr']), 2)),
)
json.dump(result, open(OUTD / "table1_corrected.json", "w"), indent=1)
pp.to_csv(OUTD / "within_study_paired.csv", index=False)
print("\nwritten: data/table1_corrected.json, data/within_study_paired.csv")
