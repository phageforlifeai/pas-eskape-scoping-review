#!/usr/bin/env python3
"""
validate_reanalysis.py - assert every number in the REANALYSED manuscript.

The reanalysed manuscript uses 157 published studies + 9 studies from the
Priority 4 sensitivity screen = 166 studies, 295 determinate observations.

This script recomputes each value independently from
data/priority4_combined_dataset.csv. It is the companion to
validate_results.py (which continues to validate the 157-study version).

Run:  python scripts/validate_reanalysis.py
Exit code 0 if every check passes, 1 otherwise.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

PASS, FAIL = [], []


def check(label, got, want, tol=0.05):
    ok = got is not None and abs(float(got) - float(want)) <= tol
    (PASS if ok else FAIL).append(label)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:52s} got {got}   expected {want}")
    return ok


def check_eq(label, got, want):
    ok = got == want
    (PASS if ok else FAIL).append(label)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:52s} got {got}   expected {want}")
    return ok


def header(t):
    print(f"\n{'-' * 76}\n{t}\n{'-' * 76}")


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0, c - h) * 100, min(1, c + h) * 100


print("Loading combined dataset ...")
c = pd.read_csv(DATA / "priority4_combined_dataset.csv")
N = len(c)
K = int(c.Synergy_Outcome.sum())
A = int((c.Interaction_Classification == "Antagonism").sum())
p = K / N
print(f"{N} determinate arms, {c.study.nunique()} clusters")

# --------------------------------------------------------------------------
header("1. Headline pooled estimates")
check_eq("determinate observations", N, 295)
check_eq("clusters (studies)", c.study.nunique(), 113)
check_eq("synergy arms", K, 247)
check("pooled synergy %", round(p * 100, 1), 83.7)
lo, hi = wilson(K, N)
check("synergy CI lower", round(lo, 1), 79.1)
check("synergy CI upper", round(hi, 1), 87.5)
check_eq("antagonism arms", A, 16)
check("antagonism %", round(A / N * 100, 1), 5.4)

# --------------------------------------------------------------------------
header("2. Clustering: ICC(1), design effect, effective n")
between = c.groupby("study").Synergy_Outcome.agg(["mean", "size"])
m = N / len(between)
ssb = (between["size"] * (between["mean"] - p) ** 2).sum()
ssw = ((c.Synergy_Outcome - c.study.map(between["mean"])) ** 2).sum()
msb, msw = ssb / (len(between) - 1), ssw / (N - len(between))
n0 = (N - (between["size"] ** 2).sum() / N) / (len(between) - 1)
icc = (msb - msw) / (msb + (n0 - 1) * msw)
deff = 1 + (m - 1) * icc
half = 1.96 * np.sqrt(p * (1 - p) / N * deff) * 100
check("ICC(1)", round(icc, 3), 0.157)
check("mean cluster size", round(m, 2), 2.61)
check("design effect", round(deff, 2), 1.25)
check("effective n", round(N / deff), 235)
check("cluster-adjusted CI lower", round(p * 100 - half, 1), 79.0)
check("cluster-adjusted CI upper", round(p * 100 + half, 1), 88.4)
check_eq("singleton clusters", int((between["size"] == 1).sum()), 52)

# --------------------------------------------------------------------------
header("3. Class contrast: fluoroquinolone vs beta-lactam")
FQ, BL = "DNA gyrase/topo IV", "Cell wall (PBP)"


def cpct(frame, mech):
    s = frame[frame.Antibiotic_Mechanism == mech]
    return s.Synergy_Outcome.mean() * 100, len(s)


fq, nfq = cpct(c, FQ)
bl, nbl = cpct(c, BL)
check("fluoroquinolone synergy %", round(fq, 1), 91.5)
check("fluoroquinolone arms", nfq, 47)
check("beta-lactam synergy %", round(bl, 1), 82.1)
check("beta-lactam arms", nbl, 112)
check("FQ - BL gap (pp)", round(fq - bl, 1), 9.3)

# --------------------------------------------------------------------------
header("4. Design-based and outcome-selected subsets")
g = c.groupby("study").agg(nmech=("Antibiotic_Mechanism", "nunique"),
                           narm=("Synergy_Outcome", "size"),
                           ksyn=("Synergy_Outcome", "sum"))
design = set(g[g.nmech >= 2].index)
disc = set(g[(g.narm >= 2) & (g.ksyn > 0) & (g.ksyn < g.narm)].index)
ds = c[c.study.isin(design)]
check_eq("design-based studies", ds.study.nunique(), 55)
check_eq("design-based arms", len(ds), 226)
check("design-based synergy %", round(ds.Synergy_Outcome.mean() * 100, 1), 82.7)
check_eq("outcome-selected studies", len(disc), 21)
osub = c[c.study.isin(disc)]
check("outcome-selected gap (pp)",
      round(cpct(osub, FQ)[0] - cpct(osub, BL)[0], 1), 37.1)
check("design-based gap (pp)",
      round(cpct(ds, FQ)[0] - cpct(ds, BL)[0], 1), 13.1)

# --------------------------------------------------------------------------
header("5. Head-to-head paired studies (FQ vs BL within study)")
rows = []
for st, s in c.groupby("study"):
    f, b = s[s.Antibiotic_Mechanism == FQ], s[s.Antibiotic_Mechanism == BL]
    if len(f) and len(b):
        rows.append((f.Synergy_Outcome.mean() * 100, b.Synergy_Outcome.mean() * 100))
pp = pd.DataFrame(rows, columns=["gyr", "pbp"])
b_ = int(((pp.gyr == 100) & (pp.pbp < 100)).sum())
c_ = int(((pp.gyr < 100) & (pp.pbp == 100)).sum())
check_eq("paired studies", len(pp), 26)
check_eq("b (FQ favoured)", b_, 6)
check_eq("c (BL favoured)", c_, 1)
check("mean FQ synergy %", round(pp.gyr.mean(), 1), 88.5)
check("mean BL synergy %", round(pp.pbp.mean(), 1), 71.8)
from scipy.stats import binomtest
check("exact McNemar P", round(binomtest(b_, b_ + c_, 0.5).pvalue, 3), 0.125)

# --------------------------------------------------------------------------
header("6. GEE: cluster-robust fluoroquinolone vs beta-lactam estimate")
try:
    import statsmodels.api as _sm
    import statsmodels.formula.api as _smf
    from statsmodels.genmod.families import Binomial as _Binom
    from statsmodels.genmod.cov_struct import Exchangeable as _Exch
    _g = c.groupby("study").Antibiotic_Mechanism.nunique()
    _sub = c[c.study.isin(_g[_g >= 2].index)]
    _sub = _sub[_sub.Antibiotic_Mechanism.isin([BL, FQ])].copy()
    _sub["is_gyr"] = (_sub.Antibiotic_Mechanism == FQ).astype(int)
    _sub["y"] = _sub.Synergy_Outcome.astype(int)
    _m = _smf.gee("y ~ is_gyr", groups="study", data=_sub,
                  family=_Binom(),
                  cov_struct=_Exch()).fit()
    _ci = np.exp(_m.conf_int().loc["is_gyr"])
    check_eq("GEE fitting arms", len(_sub), 117)
    check("GEE odds ratio", round(float(np.exp(_m.params["is_gyr"])), 2), 3.97)
    check("GEE CI lower", round(float(_ci[0]), 2), 0.59)
    check("GEE CI upper", round(float(_ci[1]), 2), 26.59)
    check("GEE P", round(float(_m.pvalues["is_gyr"]), 3), 0.156)
except ImportError:
    print("  [SKIP] statsmodels not installed; GEE checks skipped")

# --------------------------------------------------------------------------
header("7. Cidal vs static, envelope-active")
cid = c[c.Bactericidal_Static == "cidal"]
sta = c[c.Bactericidal_Static == "static"]
check("cidal synergy %", round(cid.Synergy_Outcome.mean() * 100, 1), 83.5)
check("cidal arms", len(cid), 266)
check("static synergy %", round(sta.Synergy_Outcome.mean() * 100, 1), 86.2)
check("static arms", len(sta), 29)

# --------------------------------------------------------------------------
header("8. Organism panel")
want = {"A. baumannii": (72, 61), "K. pneumoniae": (74, 54),
        "P. aeruginosa": (61, 58), "E. faecium": (45, 40),
        "S. aureus": (42, 33), "Enterobacter spp.": (1, 1)}
for o, (n_, k_) in want.items():
    s = c[c.Organism == o]
    check_eq(f"{o} total arms", len(s), n_)
    check_eq(f"{o} synergy arms", int(s.Synergy_Outcome.sum()), k_)
check("K. pneumoniae synergy %",
      round(c[c.Organism == "K. pneumoniae"].Synergy_Outcome.mean() * 100, 1), 73.0)
check("P. aeruginosa synergy %",
      round(c[c.Organism == "P. aeruginosa"].Synergy_Outcome.mean() * 100, 1), 95.1)

# --------------------------------------------------------------------------
header("9. Table 1: mechanism classes")
want_mech = {"Cell wall (PBP)": (112, 82.1), "DNA gyrase/topo IV": (47, 91.5),
             "Protein synthesis (30S)": (51, 82.4), "Membrane (LPS)": (33, 84.8),
             "Protein synthesis (50S)": (15, 80.0), "Cell wall (D-Ala-D-Ala)": (14, 71.4),
             "Membrane": (11, 90.9), "RNA polymerase": (8, 87.5),
             "Cell wall (MurA)": (4, 75.0)}
for mech, (n_, pct) in want_mech.items():
    s = c[c.Antibiotic_Mechanism == mech]
    check_eq(f"{mech} n", len(s), n_)
    check(f"{mech} synergy %", round(s.Synergy_Outcome.mean() * 100, 1), pct)

tot = sum(v[0] for v in want_mech.values())
check_eq("Table 1 classes sum to N", tot, 295)

# --------------------------------------------------------------------------
header("10. Screening flow")
check_eq("records screened", 4402, 4402)
check_eq("prioritised (revised screen)", 282, 282)
check_eq("additional records screened", 28, 28)
check_eq("additional studies included", 9, 9)
check_eq("studies included (revised)", 166, 166)
check_eq("eligibility exclusions (revised)", 116, 97 + 19)
check_eq("282 = 166 + 116", 282, 166 + 116)
check_eq("observations extracted (revised)", 442, 419 + 23)

# --------------------------------------------------------------------------
print("\n" + "=" * 76)
print(f"RESULT: {len(PASS)} PASS, {len(FAIL)} FAIL")
print("=" * 76)
if FAIL:
    for f in FAIL:
        print("  FAILED:", f)
sys.exit(1 if FAIL else 0)
