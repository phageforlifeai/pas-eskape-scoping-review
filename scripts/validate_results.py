#!/usr/bin/env python3
"""
validate_results.py - independently recompute and assert every headline number
reported in the manuscript.

Run:  python scripts/validate_results.py
Exit code 0 if every check passes, 1 otherwise.

This script is deliberately self-contained: it depends only on
data/eskape_pas_dataset_clean.csv, data/crkp_pas_dataset_clean.csv and
data/within_study_paired.csv, and re-derives everything from scratch rather
than reading table1_corrected.json.
"""
import sys
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

PASS, FAIL, WARN = [], [], []


def check(label, got, want, tol=0.05, warn_only=False):
    """Assert got == want within tol. Records PASS / FAIL / WARN."""
    ok = (got is not None) and (abs(float(got) - float(want)) <= tol)
    tag = "PASS" if ok else ("WARN" if warn_only else "FAIL")
    (PASS if ok else (WARN if warn_only else FAIL)).append(label)
    print(f"  [{tag}] {label:58s} got {got}   expected {want}")
    return ok


def check_eq(label, got, want, warn_only=False):
    ok = got == want
    tag = "PASS" if ok else ("WARN" if warn_only else "FAIL")
    (PASS if ok else (WARN if warn_only else FAIL)).append(label)
    print(f"  [{tag}] {label:58s} got {got}   expected {want}")
    return ok


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0, c - h) * 100, min(1, c + h) * 100


def header(t):
    print(f"\n{'-' * 78}\n{t}\n{'-' * 78}")


# --------------------------------------------------------------------------
print("Loading data ...")
d = pd.read_csv(DATA / "eskape_pas_dataset_clean.csv")
det = d[d.Synergy_Outcome.notna()].copy()
det["study"] = det.PMID_or_DOI.astype(str)
print(f"  dataset rows            : {len(d)}")
print(f"  determinate observations: {len(det)}")
print(f"  distinct studies        : {det.study.nunique()}")

# --------------------------------------------------------------------------
header("1. Evidence base  (Methods / Results 3.1)")
check_eq("eligible studies", d.PMID_or_DOI.nunique(), 157)
check_eq("antibiotic-level observations", len(d), 419)
check_eq("determinate outcomes", len(det), 272)
check("determinate proportion (%)", round(len(det) / len(d) * 100, 1), 64.9)
check_eq("study-organism entries with determinate outcome", det.study.nunique(), 105)

# --------------------------------------------------------------------------
header("2. Outcome distribution  (Results 3.2)")
for cls, want_n, want_pct in [("Synergy", 232, 85.3), ("Indifferent", 21, 7.7),
                              ("Additive", 4, 1.5), ("Antagonism", 15, 5.5)]:
    n = int((det.Interaction_Classification == cls).sum())
    check_eq(f"{cls} count", n, want_n)
    check(f"{cls} percent", round(n / len(det) * 100, 1), want_pct)
check_eq("outcome counts sum to determinate total",
         int(det.Interaction_Classification.isin(
             ["Synergy", "Indifferent", "Additive", "Antagonism"]).sum()), 272)

# --------------------------------------------------------------------------
header("3. Organism breakdown  (Results 3.1)")
ORGS = [("A. baumannii", 65, 24), ("K. pneumoniae", 61, 19), ("P. aeruginosa", 58, 25),
        ("E. faecium", 45, 11), ("S. aureus", 42, 25), ("Enterobacter spp.", 1, 1)]
tot_obs = tot_st = 0
for org, want_n, want_st in ORGS:
    s = det[det.Organism == org]
    check_eq(f"{org}: observations", len(s), want_n)
    check_eq(f"{org}: studies", s.study.nunique(), want_st)
    tot_obs += len(s)
    tot_st += s.study.nunique()
check_eq("organism observations sum", tot_obs, 272)
check_eq("organism studies sum", tot_st, 105)

# --------------------------------------------------------------------------
header("4. Pooled synergy by mechanism class + Wilson 95% CI  (Table 1, col 2)")
# row: mechanism, pooled %, n, CI lo, CI hi, within-study %, n studies, arms %, delta pp, antagonism
TABLE1 = [
    ("DNA gyrase/topo IV", 93.0, 43, 81.4, 97.6, 92.9, 28, 94.1, -0.1, (1, 43)),
    ("Membrane", 90.9, 11, 62.3, 98.4, 88.9, 9, 88.9, -2.0, (0, 11)),
    ("RNA polymerase", 87.5, 8, 52.9, 97.8, 87.5, 8, 87.5, 0.0, (1, 8)),
    ("Membrane (LPS)", 86.7, 30, 70.3, 94.7, 83.3, 21, 81.8, -3.4, (1, 30)),
    ("Cell wall (PBP)", 84.8, 99, 76.5, 90.6, 82.1, 39, 81.7, -2.7, (5, 99)),
    ("Protein synthesis (30S)", 83.3, 48, 70.4, 91.3, 85.6, 30, 87.5, 2.3, (4, 48)),
    ("Protein synthesis (50S)", 80.0, 15, 54.8, 93.0, 83.3, 11, 78.6, 3.3, (2, 15)),
    ("Cell wall (MurA)", 75.0, 4, 30.1, 95.4, 75.0, 4, 75.0, 0.0, (0, 4)),
    ("Cell wall (D-Ala-D-Ala)", 71.4, 14, 45.4, 88.3, 71.4, 7, 71.4, 0.0, (1, 14)),
]

# design-based subset: studies testing >= 2 distinct mechanisms (selected on design)
gm = det.groupby("study").Antibiotic_Mechanism.nunique()
DESIGN = set(gm[gm >= 2].index)
dsub = det[det.study.isin(DESIGN)]
check_eq("design-based subset: studies", len(DESIGN), 50)
check_eq("design-based subset: observations", len(dsub), 209)

sum_pool = sum_arms = sum_ant = 0
for mech, p_pct, p_n, lo, hi, ws, ws_st, arms_pct, delta, (k_ant, n_ant) in TABLE1:
    a = det[det.Antibiotic_Mechanism == mech]
    b = dsub[dsub.Antibiotic_Mechanism == mech]
    k = int(a.Synergy_Outcome.sum())
    w_lo, w_hi = wilson(k, len(a))
    per = b.groupby("study").Synergy_Outcome.mean() * 100
    short = mech[:26]
    check_eq(f"{short}: pooled n", len(a), p_n)
    check(f"{short}: pooled %", round(k / len(a) * 100, 1), p_pct)
    check(f"{short}: CI lower", round(w_lo, 1), lo)
    check(f"{short}: CI upper", round(w_hi, 1), hi)
    if len(per):
        check(f"{short}: within-study %", round(per.mean(), 1), ws)
        check_eq(f"{short}: within-study studies", len(per), ws_st)
        check(f"{short}: arms %", round(b.Synergy_Outcome.mean() * 100, 1), arms_pct)
        check(f"{short}: delta pp", round(per.mean() - k / len(a) * 100, 1), delta, tol=0.15)
    n_ant_got = int((a.Interaction_Classification == "Antagonism").sum())
    check_eq(f"{short}: antagonism", (n_ant_got, len(a)), (k_ant, n_ant))
    sum_pool += len(a)
    sum_arms += len(b)
    sum_ant += n_ant_got

header("5. Table 1 totals and column sums")
check_eq("pooled column sums to", sum_pool, 272)
check_eq("within-study column sums to", sum_arms, 209)
check_eq("antagonism column sums to", sum_ant, 15)
check("All classes: pooled %", round(det.Synergy_Outcome.mean() * 100, 1), 85.3)
check("All classes: design-based arms %", round(dsub.Synergy_Outcome.mean() * 100, 1), 84.7)
check("All classes: delta pp", round(dsub.Synergy_Outcome.mean() * 100
                                    - det.Synergy_Outcome.mean() * 100, 1), -0.6)

# --------------------------------------------------------------------------
header("6. Clustering: ICC(1), design effect, effective n, cluster CI  (Results 3.2)")
k, n = det.Synergy_Outcome.sum(), len(det)
p = k / n
between = det.groupby("study").Synergy_Outcome.agg(["mean", "size"])
m_bar = n / det.study.nunique()
naive_lo, naive_hi = wilson(int(k), n)

# --- one-way ANOVA ICC(1): the estimator now used in the manuscript --------
ssb = (between["size"] * (between["mean"] - p) ** 2).sum()
ssw = ((det.Synergy_Outcome - det.study.map(between["mean"])) ** 2).sum()
K = len(between)
msb, msw = ssb / (K - 1), ssw / (n - K)
n0 = (n - (between["size"] ** 2).sum() / n) / (K - 1)
icc = (msb - msw) / (msb + (n0 - 1) * msw)
deff = 1 + (m_bar - 1) * icc
half = 1.96 * np.sqrt(p * (1 - p) / n * deff) * 100

check("ICC(1) one-way ANOVA", round(icc, 2), 0.16)
check("mean cluster size", round(m_bar, 2), 2.59)
check("design effect", round(deff, 2), 1.25)
check("effective n", round(n / deff), 218)
check("naive Wilson CI lower", round(naive_lo, 1), 80.6)
check("naive Wilson CI upper", round(naive_hi, 1), 89.0)
check("cluster-adjusted CI lower", round(p * 100 - half, 1), 80.6)
check("cluster-adjusted CI upper", round(p * 100 + half, 1), 90.0)

print("\n  -- eta-squared is NOT an ICC and is no longer used in the design effect --")
eta2 = ssb / ((det.Synergy_Outcome - p) ** 2).sum()
print(f"  [INFO] eta^2 (variance explained by study)  = {eta2:.3f}  "
      f"- descriptive only; would give deff {1 + (m_bar - 1) * eta2:.2f} (not used)")
null_eta2 = (K - 1) / (n - 1)
print(f"  [INFO] expected eta^2 with no clustering    = {null_eta2:.3f}  "
      f"(so {eta2:.3f} is only {(eta2 - null_eta2) * 100:.1f} pp above chance)")
check("cluster sizes min", int(between["size"].min()), 1)
check("cluster sizes max", int(between["size"].max()), 10)
check_eq("singleton clusters", int((between["size"] == 1).sum()), 50)

# --------------------------------------------------------------------------
header("7. Multifactorial and discordance counts  (Results 3.2)")
gg = det.groupby("study").agg(nab=("Antibiotic_Tested", "nunique"),
                              nn=("Synergy_Outcome", "size"),
                              kk=("Synergy_Outcome", "sum"))
multi = gg[gg.nab >= 2]
disc = multi[(multi.kk > 0) & (multi.kk < multi.nn)]
check_eq("studies testing >=2 antibiotics", len(multi), 55)
check_eq("discordant among them", len(disc), 17)
check("discordant percent", round(len(disc) / len(multi) * 100, 1), 30.9)

# --------------------------------------------------------------------------
header("8. Outcome-selection effect: the 8.2 / 12.4 / 40.7 pp contrast  (Results 3.3)")
PBP, GYR = "Cell wall (PBP)", "DNA gyrase/topo IV"


def gap(frame):
    b = frame[frame.Antibiotic_Mechanism == PBP].Synergy_Outcome
    g = frame[frame.Antibiotic_Mechanism == GYR].Synergy_Outcome
    return g.mean() * 100 - b.mean() * 100


DISC = set(disc.index)
check("pooled gap (pp)", round(gap(det), 1), 8.2, tol=0.1)
check("arm-based design gap (pp)", round(gap(dsub), 1), 12.4, tol=0.1)
check("outcome-selected gap (pp)", round(gap(det[det.study.isin(DISC)]), 1), 40.7, tol=0.15)
sel = det[det.study.isin(DISC)]
check("outcome-selected beta-lactam %", round(
    sel[sel.Antibiotic_Mechanism == PBP].Synergy_Outcome.mean() * 100, 1), 48.1, tol=0.2)
check("'five-fold' ratio 40.7/8.2", round(40.7 / 8.2, 1), 5.0, tol=0.1)
check("study-level class-mean gap (92.9-82.1)", round(92.9 - 82.1, 1), 10.8, tol=0.1)
print("  [INFO] note: 10.8 = difference between the two study-level class means "
      "(92.9 - 82.1),")
print("  [INFO]       NOT the mean of paired within-study differences (which is +0.2 pp)")

# --------------------------------------------------------------------------
header("9. Head-to-head paired comparison  (Results 3.3)")
pp = pd.read_csv(DATA / "within_study_paired.csv")
pp["study"] = pp.study.astype(str)
check_eq("studies measuring both classes", len(pp), 23)
check("mean per-study synergy, fluoroquinolone", round(pp.gyr.mean(), 1), 91.3)
check("mean per-study synergy, beta-lactam", round(pp.pbp.mean(), 1), 75.4)
check("paired difference (pp)", round(pp.gyr.mean() - pp.pbp.mean(), 1), 15.9)
b = int(((pp.gyr == 100) & (pp.pbp < 100)).sum())
c = int(((pp.gyr < 100) & (pp.pbp == 100)).sum())
check_eq("discordant pairs (b, c)", (b, c), (5, 1))
check("exact McNemar P", round(stats.binomtest(b, b + c, 0.5).pvalue, 2), 0.22)

# the partition problem flagged in the verification report
both100 = int(((pp.gyr == 100) & (pp.pbp == 100)).sum())
neither = int(((pp.gyr < 100) & (pp.pbp < 100)).sum())
both_pos = int(((pp.gyr > 0) & (pp.pbp > 0)).sum())
print(f"\n  -- discordance partition --")
print(f"  [INFO] strictly exclusive: {b} favour FQ + {c} favour PBP + "
      f"{both100} both at 100% + {neither} neither = {b + c + both100 + neither}")
print(f"  [INFO] 'reported synergy with both' (both >0%) = {both_pos}, which overlaps")
print(f"  [INFO] with the '{b} favour FQ' group by one study -> see KNOWN_DISCREPANCIES")
check_eq("exclusive partition sums to n", b + c + both100 + neither, 23)

# GEE
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    sub = dsub[dsub.Antibiotic_Mechanism.isin([PBP, GYR])].copy()
    sub["is_gyr"] = (sub.Antibiotic_Mechanism == GYR).astype(int)
    sub["y"] = sub.Synergy_Outcome.astype(int)
    gee = smf.gee("y ~ is_gyr", groups="study", data=sub,
                  family=sm.families.Binomial(),
                  cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = np.exp(gee.conf_int().loc["is_gyr"])
    check("GEE OR", round(float(np.exp(gee.params["is_gyr"])), 1), 4.6, tol=0.1)
    check("GEE CI lower", round(float(ci[0]), 2), 0.45, tol=0.02)
    check("GEE CI upper", round(float(ci[1]), 1), 47.3, tol=0.3)
except ImportError:
    print("  [SKIP] statsmodels not installed - GEE checks skipped")

# --------------------------------------------------------------------------
header("10. Bactericidal vs bacteriostatic  (Results 3.4)")
CIDAL = {"Fluoroquinolone", "Aminoglycoside", "Polymyxin", "Carbapenem", "Cephalosporin",
         "Penicillin", "Glycopeptide", "Lipopeptide", "Cephalosporin/BLI", "Penicillin/BLI",
         "Rifamycin", "Siderophore cephalosporin", "Fosfomycin"}
cid = det[det.Antibiotic_Class.isin(CIDAL)]
stat_ = det[~det.Antibiotic_Class.isin(CIDAL)]
check_eq("bactericidal n", len(cid), 243)
check("bactericidal synergy %", round(cid.Synergy_Outcome.mean() * 100, 1), 85.2)
check_eq("bacteriostatic n", len(stat_), 29)
check("bacteriostatic synergy %", round(stat_.Synergy_Outcome.mean() * 100, 1), 86.2)
tab = [[int(cid.Synergy_Outcome.sum()), len(cid) - int(cid.Synergy_Outcome.sum())],
       [int(stat_.Synergy_Outcome.sum()), len(stat_) - int(stat_.Synergy_Outcome.sum())]]
check("Fisher exact P", round(stats.fisher_exact(tab)[1], 2), 1.00, tol=0.01)

# --------------------------------------------------------------------------
header("11. Beta-lactam hypothesis specificity tests  (Results 3.7)")
bl = det[det.Antibiotic_Mechanism == PBP]
ot = det[det.Antibiotic_Mechanism != PBP]
kb, ko = int(bl.Synergy_Outcome.sum()), int(ot.Synergy_Outcome.sum())
check_eq("beta-lactam synergy", (kb, len(bl)), (84, 99))
check_eq("all other agents synergy", (ko, len(ot)), (148, 173))
check("Fisher exact P (beta-lactam vs other)",
      round(stats.fisher_exact([[kb, len(bl) - kb], [ko, len(ot) - ko]])[1], 2), 0.86, tol=0.01)

ENVELOPE = {PBP, "Cell wall (D-Ala-D-Ala)", "Cell wall (MurA)", "Membrane", "Membrane (LPS)"}
env = det[det.Antibiotic_Mechanism.isin(ENVELOPE)]
oth = det[~det.Antibiotic_Mechanism.isin(ENVELOPE)]
ea = int((env.Interaction_Classification == "Antagonism").sum())
oa = int((oth.Interaction_Classification == "Antagonism").sum())
check_eq("envelope-active antagonism", (ea, len(env)), (7, 158))
check_eq("other antagonism", (oa, len(oth)), (8, 114))
check("Fisher exact P (envelope vs other)",
      round(stats.fisher_exact([[ea, len(env) - ea], [oa, len(oth) - oa]])[1], 2), 0.42, tol=0.01)

# --------------------------------------------------------------------------
header("12. Ranges and assay breakdown  (Results 3.1, 3.5, Discussion)")
# exclude organisms with a single observation (Enterobacter spp., n=1): the
# manuscript reports these separately as precluding species-level inference
org_counts = det.groupby("Organism").Synergy_Outcome.size()
syn_by_org = (det.groupby("Organism").Synergy_Outcome.mean() * 100)[org_counts >= 10]
check("lowest organism synergy (S. aureus)", round(syn_by_org.min(), 1), 78.6, tol=0.1)
check_eq("  ...which organism", syn_by_org.idxmin(), "S. aureus")
check("highest organism synergy (P. aeruginosa)", round(syn_by_org.max(), 1), 94.8, tol=0.1)
check_eq("  ...which organism", syn_by_org.idxmax(), "P. aeruginosa")
check_eq("antagonistic observations", int((det.Interaction_Classification == "Antagonism").sum()), 15)
check_eq("30S antagonism", (int(((det.Antibiotic_Mechanism == "Protein synthesis (30S)")
                                 & (det.Interaction_Classification == "Antagonism")).sum()),
                            int((det.Antibiotic_Mechanism == "Protein synthesis (30S)").sum())), (4, 48))
print("\n  -- assay counts (primary-assay classification) --")
first = det.Assay_Type.astype(str).str.split(";").str[0].str.strip()
assay_got = {a: int((first == a).sum()) for a in
             ["Checkerboard", "Biofilm", "Time-Kill", "In Vivo"]}
for a, want in [("Checkerboard", 170), ("Biofilm", 45), ("Time-Kill", 37), ("In Vivo", 17)]:
    check_eq(f"primary assay: {a}", assay_got[a], want)
check_eq("primary-assay counts sum", sum(assay_got.values()), 269)
print("  [INFO] these four counts are EXCLUSIVE (primary assay per arm) and sum to 269,")
print("  [INFO] not 272: 3 arms (PMID 33649110) carry no assay label. The manuscript")
print("  [INFO] calls them 'non-exclusive' - see KNOWN_DISCREPANCIES.")
nonexcl = {a: int(det.Assay_Types_All.astype(str).str.contains(a, regex=False).sum())
           for a in ["Checkerboard", "Biofilm", "Time-Kill", "In Vivo"]}
print(f"  [INFO] true non-exclusive counts would be {nonexcl} (sum {sum(nonexcl.values())})")

# --------------------------------------------------------------------------
header("13. CRKP continuity check  (Results 3.2)")
crkp = pd.read_csv(DATA / "crkp_pas_dataset_clean.csv")
cdet = crkp[crkp.Synergy_Outcome.notna()]
check_eq("CRKP observations", len(crkp), 74)
check_eq("CRKP studies", crkp.PMID_or_DOI.nunique(), 35)
check("CRKP pooled synergy %", round(cdet.Synergy_Outcome.mean() * 100, 1), 92.6)

# --------------------------------------------------------------------------
header("14. Screening and eligibility flow  (PRISMA-ScR)")
sa = pd.read_csv(DATA / "screening_audit_4402.csv")
check_eq("unique records screened", len(sa), 4402)
check_eq("records prioritised (score >= 14)", int((sa.prioritised_full_text == "Y").sum()), 254)
print(f"  [INFO] 254 = 157 included studies + 97 eligibility exclusions")
ex = pd.read_csv(DATA / "eligibility_exclusions_97.csv")
check_eq("eligibility exclusions", len(ex), 97)
check_eq("exclusions + included = prioritised", len(ex) + 157, 254)
reasons = ex.exclusion_reason.value_counts().to_dict()
check_eq("  reason: no antibiotic in a phage-combination context",
         reasons.get("no antibiotic in a phage-combination context"), 77)
check_eq("  reason: review prose, no primary experiment",
         reasons.get("review prose, no primary experiment"), 20)

# --------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"RESULT: {len(PASS)} PASS, {len(FAIL)} FAIL, {len(WARN)} WARN")
if FAIL:
    print("\nFAILED CHECKS:")
    for f in FAIL:
        print(f"  - {f}")
print("=" * 78)
sys.exit(1 if FAIL else 0)
