#!/usr/bin/env python3
"""
Step 18 - Reconcile Table 1 against the extraction dataset, and rebuild the
within-study analysis with a defensible statistical design.

Part A  audit every cell of Table 1 against eskape_pas_dataset_clean.csv
Part B  diagnose the multiple-observation problem
Part C  compare candidate within-study designs
Part D  report the primary analysis with study-level clustering
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
OUTD = ROOT / "data"
OUTD.mkdir(exist_ok=True)

d = pd.read_csv(ROOT / "eskape_pas_dataset_clean.csv")
det = d[d.Synergy_Outcome.notna()].copy()
det["study"] = det.PMID_or_DOI.astype(str)

MECH = ["DNA gyrase/topo IV", "Protein synthesis (30S)", "Protein synthesis (50S)",
        "RNA polymerase", "Membrane (LPS)", "Cell wall (PBP)", "Cell wall (D-Ala-D-Ala)"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0, c - h) * 100, min(1, c + h) * 100)


print("=" * 78)
print("PART A  AUDIT OF TABLE 1 AGAINST THE EXTRACTION DATASET")
print("=" * 78)
print(f"records with a determinate outcome: {len(det)}")
print(f"distinct studies                   : {det.study.nunique()}")
print(f"outcome distribution               : "
      + ", ".join(f"{k}={v}" for k, v in
                  det.Interaction_Classification.value_counts().items()))
print(f"synergy rate, all determinate      : {det.Synergy_Outcome.mean()*100:.1f}%")

# reproduce the manuscript's discordant subset exactly as described in 2.5
g = det.groupby("study").agg(n=("Synergy_Outcome", "size"),
                              k=("Synergy_Outcome", "sum"),
                              n_ab=("Antibiotic_Tested", "nunique"))
multi_ab = g[g.n_ab >= 2]
disc = multi_ab[(multi_ab.k > 0) & (multi_ab.k < multi_ab.n)]
ws = det[det.study.isin(disc.index)]

print(f"\nstudies testing >=2 distinct antibiotics : {len(multi_ab)}")
print(f"  of which discordant                    : {len(disc)} "
      f"({len(disc)/len(multi_ab)*100:.1f}%)")
print(f"  observations in discordant subset      : {len(ws)}")
print(f"  synergy in discordant subset           : {ws.Synergy_Outcome.mean()*100:.1f}%")

TABLE1 = {
    "DNA gyrase/topo IV":        (93.0, 43, 88.9, 18,  "1/43"),
    "RNA polymerase":            (87.5,  8, 66.7,  3,  "1/8"),
    "Membrane (LPS)":            (86.7, 30, 55.6,  9,  "1/30"),
    "Cell wall (PBP)":           (84.8, 99, 48.1, 27,  "5/99"),
    "Protein synthesis (30S)":   (83.3, 48, 68.8, 16,  "4/48"),
    "Protein synthesis (50S)":   (80.0, 15, 66.7,  9,  "2/15"),
    "Cell wall (D-Ala-D-Ala)":   (71.4, 14, 33.3,  3,  "1/14"),
}

print("\nmechanism                     pooled (manuscript -> recomputed)      "
      "within-study (manuscript -> recomputed)")
audit_rows = []
for m in MECH:
    a = det[det.Antibiotic_Mechanism == m]
    b = ws[ws.Antibiotic_Mechanism == m]
    m_p, m_n, m_w, m_wn, m_ant = TABLE1[m]
    p_pool = a.Synergy_Outcome.mean() * 100
    p_ws = b.Synergy_Outcome.mean() * 100 if len(b) else np.nan
    k_ant = int((a.Interaction_Classification == "Antagonism").sum())
    ok_pool = abs(p_pool - m_p) < 0.15 and len(a) == m_n
    ok_ws = (len(b) == m_wn) and (np.isnan(p_ws) or abs(p_ws - m_w) < 0.15)
    ok_ant = m_ant == f"{k_ant}/{len(a)}"
    flag = "OK" if (ok_pool and ok_ws and ok_ant) else "MISMATCH"
    print(f"{m:28s} {m_p:5.1f} (n={m_n:3d}) -> {p_pool:5.1f} (n={len(a):3d})   "
          f"{m_w:5.1f} (n={m_wn:3d}) -> {p_ws:5.1f} (n={len(b):3d})   {flag}")
    audit_rows.append(dict(mech=m, man_pooled=m_p, man_n=m_n,
                           rec_pooled=round(p_pool, 1), rec_n=len(a),
                           man_ws=m_w, man_wn=m_wn,
                           rec_ws=round(p_ws, 1) if len(b) else None, rec_wn=len(b),
                           antagonism=f"{k_ant}/{len(a)}", status=flag))

# headline totals quoted in the manuscript
print("\n--- headline numbers quoted in the text ---")
checks = {
    "272 determinate calls": len(det) == 272,
    "232 synergistic": int((det.Interaction_Classification == "Synergy").sum()) == 232,
    "21 indifferent": int((det.Interaction_Classification == "Indifferent").sum()) == 21,
    "4 additive": int((det.Interaction_Classification == "Additive").sum()) == 4,
    "15 antagonistic": int((det.Interaction_Classification == "Antagonism").sum()) == 15,
    "85.3% overall synergy": abs(det.Synergy_Outcome.mean() * 100 - 85.3) < 0.1,
    "88 observations in discordant subset": len(ws) == 88,
    "62.5% synergy in discordant subset": abs(ws.Synergy_Outcome.mean() * 100 - 62.5) < 0.1,
    "55 studies testing >=2 antibiotics": len(multi_ab) == 55,
    "17 discordant studies": len(disc) == 17,
}
for label, ok in checks.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {label}")

print("\n" + "=" * 78)
print("PART B  THE MULTIPLE-OBSERVATION PROBLEM")
print("=" * 78)
print(f"observations : {len(det)}")
print(f"studies      : {det.study.nunique()}")
print(f"mean observations per study: {len(det)/det.study.nunique():.1f}")
per = det.groupby("study").size()
print(f"max observations from one study: {per.max()}")
print(f"studies contributing >=10 observations: {int((per >= 10).sum())}")
print(f"share of all observations from the top 5 studies: "
      f"{per.nlargest(5).sum()/len(det)*100:.1f}%")

# design effect for a clustered binary outcome
k, n_s = det.Synergy_Outcome.sum(), len(det)
p = k / n_s
between = det.groupby("study").Synergy_Outcome.agg(["mean", "size"])
icc_num = ((between["mean"] - p) ** 2 * between["size"]).sum()
icc_den = (n_s - 1) * p * (1 - p)
icc = icc_num / icc_den
m_bar = n_s / det.study.nunique()
deff = 1 + (m_bar - 1) * icc
print(f"\nintra-study correlation coefficient (ICC) approx : {icc:.3f}")
print(f"mean cluster size                                : {m_bar:.2f}")
print(f"design effect                                    : {deff:.2f}")
print(f"effective sample size for the pooled 85.3%       : {n_s/deff:.1f} "
      f"(nominal {n_s})")

ci_naive = wilson(k, n_s)
z = 1.96
half = z * np.sqrt(p * (1 - p) / n_s * deff) * 100
print(f"naive Wilson 95% CI          : {ci_naive[0]:.1f}-{ci_naive[1]:.1f}%")
print(f"cluster-adjusted 95% CI      : {p*100-half:.1f}-{p*100+half:.1f}%")

# how concentrated is each mechanism in single studies?
print("\nconcentration of evidence by mechanism (largest single-study share):")
for m in MECH:
    s = det[det.Antibiotic_Mechanism == m]
    biggest = s.groupby("study").size().max()
    print(f"  {m:28s} n={len(s):3d}  largest study contributes {biggest:2d} "
          f"({biggest/len(s)*100:4.1f}%)")

json.dump(audit_rows, open(OUTD / "table1_audit.json", "w"), indent=1)

print("\n" + "=" * 78)
print("PART C  CANDIDATE WITHIN-STUDY DESIGNS")
print("=" * 78)

gmech = det.groupby("study").Antibiotic_Mechanism.nunique()
multi_mech = gmech[gmech >= 2]

designs = {}

# design 0: the manuscript's outcome-selected discordant subset (for comparison)
designs["0. Outcome-selected (manuscript)"] = disc.index

# design 1: design-selected - any study testing >=2 distinct mechanisms
designs["1. Design-selected (>=2 mechanisms)"] = multi_mech.index

# design 2: design-selected AND internally discordant (conservative)
both = [s for s in multi_mech.index if s in disc.index]
designs["2. >=2 mechanisms AND discordant"] = pd.Index(both)

# design 3: >=2 mechanisms, analysed as study-level difference (paired)
print("design                                   studies  obs   synergy%")
for name, idx in designs.items():
    sub = det[det.study.isin(idx)]
    print(f"{name:40s} {len(idx):5d} {len(sub):5d}   {sub.Synergy_Outcome.mean()*100:5.1f}")

print("\n--- per-mechanism rates under each design ---")
print(f"{'mechanism':28s}" + "".join(f"{k.split('.')[0]:>10s}" for k in designs))
rows_c = []
for m in MECH:
    line = f"{m:28s}"
    rec = {"mech": m}
    for name, idx in designs.items():
        sub = det[(det.study.isin(idx)) & (det.Antibiotic_Mechanism == m)]
        v = sub.Synergy_Outcome.mean() * 100 if len(sub) else np.nan
        rec[name] = (round(v, 1) if not np.isnan(v) else None, len(sub))
        line += f"{v:6.1f}({len(sub):2d})" if not np.isnan(v) else f"{'-':>10s}"
    print(line)
    rows_c.append(rec)

print("\n" + "=" * 78)
print("PART D  PRIMARY ANALYSIS: STUDY AS THE UNIT OF INFERENCE")
print("=" * 78)

# For each study that tested both PBP and gyrase, compute the within-study
# difference. Study is the experimental unit, so n = number of studies.
pairs = []
for sid, s in det[det.study.isin(multi_mech.index)].groupby("study"):
    pbp = s[s.Antibiotic_Mechanism == "Cell wall (PBP)"].Synergy_Outcome
    gyr = s[s.Antibiotic_Mechanism == "DNA gyrase/topo IV"].Synergy_Outcome
    if len(pbp) and len(gyr):
        pairs.append(dict(study=sid, n_pbp=len(pbp), n_gyr=len(gyr),
                          pbp=pbp.mean() * 100, gyr=gyr.mean() * 100,
                          diff=gyr.mean() * 100 - pbp.mean() * 100))
pp = pd.DataFrame(pairs)
print(f"studies testing both PBP and gyrase: {len(pp)}")
print(f"  mean PBP synergy  : {pp.pbp.mean():.1f}%")
print(f"  mean gyrase synergy: {pp.gyr.mean():.1f}%")
print(f"  mean paired difference (gyrase - PBP): {pp['diff'].mean():+.1f} pp "
      f"(SD {pp['diff'].std():.1f})")

w_stat, w_p = stats.wilcoxon(pp.pbp, pp.gyr)
t_stat, t_p = stats.ttest_rel(pp.gyr, pp.pbp)
print(f"\nWilcoxon signed-rank: W={w_stat:.1f}, P={w_p:.3f}")
print(f"paired t-test       : t={t_stat:.2f}, df={len(pp)-1}, P={t_p:.3f}")

# bootstrap CI on the paired mean difference, resampling studies
rng = np.random.default_rng(42)
boots = [rng.choice(pp["diff"].values, len(pp), replace=True).mean()
         for _ in range(10000)]
lo, hi = np.percentile(boots, [2.5, 97.5])
print(f"bootstrap 95% CI of the difference (studies resampled): "
      f"{lo:+.1f} to {hi:+.1f} pp")
print(f"studies favouring gyrase: {int((pp['diff'] > 0).sum())}, "
      f"favouring PBP: {int((pp['diff'] < 0).sum())}, "
      f"tied: {int((pp['diff'] == 0).sum())}")

# how many studies actually discriminate?
print(f"\nstudies where the two classes truly differ (|diff|>0): "
      f"{int((pp['diff'] != 0).sum())} of {len(pp)}")

pp.to_csv(OUTD / "within_study_paired.csv", index=False)

print("\n--- full per-mechanism within-study comparison, study as unit ---")
print(f"{'mechanism':28s} {'studies':>8s} {'study-level mean':>18s} {'pooled-obs':>12s}")
for m in MECH:
    sub = det[(det.study.isin(multi_mech.index)) & (det.Antibiotic_Mechanism == m)]
    if not len(sub):
        continue
    per_study = sub.groupby("study").Synergy_Outcome.mean() * 100
    print(f"{m:28s} {per_study.size:8d} {per_study.mean():17.1f}% "
          f"{sub.Synergy_Outcome.mean()*100:11.1f}%")

res = dict(
    n_studies=len(pp),
    pbp_mean=round(float(pp.pbp.mean()), 1),
    gyr_mean=round(float(pp.gyr.mean()), 1),
    diff_mean=round(float(pp["diff"].mean()), 1),
    diff_sd=round(float(pp["diff"].std()), 1),
    wilcoxon_P=round(float(w_p), 3),
    ttest_P=round(float(t_p), 3),
    boot_lo=round(float(lo), 1),
    boot_hi=round(float(hi), 1),
    n_favour_gyr=int((pp["diff"] > 0).sum()),
    n_favour_pbp=int((pp["diff"] < 0).sum()),
    icc=round(float(icc), 3),
    deff=round(float(deff), 2),
)
json.dump(res, open(OUTD / "within_study_primary.json", "w"), indent=1)
print("\nwritten: data/table1_audit.json, data/within_study_paired.csv, "
      "data/within_study_primary.json")
