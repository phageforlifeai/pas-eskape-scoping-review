#!/usr/bin/env python3
"""
Step 21 - REANALYSIS after the Priority 4 sensitivity screen.

Rebuilds every reported statistic on the combined dataset:
  published 272 determinate observations + 23 arms from 9 newly eligible studies.

Writes scripts/reanalysis_results.json, which the manuscript/SI update step reads.
Nothing published is overwritten. The original dataset is left untouched.

Usage:  python scripts/21_reanalysis_priority4.py
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = {}

# --------------------------------------------------------------------------
# 1. Build the combined arm-level dataset
# --------------------------------------------------------------------------
raw = pd.read_csv(DATA / "eskape_pas_dataset_clean.csv")
base = raw[raw.Synergy_Outcome.notna()].copy()
base["study"] = "OLD:" + base.PMID_or_DOI.astype(str)
base = base.rename(columns={"Antibiotic_Tested": "Antibiotic_Name", "Assay_Type": "Assay"})

new = pd.read_csv(DATA / "priority4_extracted_arms.csv")
ORG = {"Acinetobacter baumannii": "A. baumannii",
       "Klebsiella pneumoniae": "K. pneumoniae",
       "Pseudomonas aeruginosa": "P. aeruginosa"}
new["Organism"] = new.organism.map(lambda x: ORG.get(x, x))
new["Synergy_Outcome"] = (new.interaction_classification == "Synergy").astype(int)
new["Interaction_Classification"] = new.interaction_classification
new["Antibiotic_Mechanism"] = new.antibiotic_mechanism
new["Antibiotic_Class"] = new.antibiotic
new["Bactericidal_Static"] = new.cidal_static.str.lower()
new["Assay"] = new.assay
new["Dosing_Sequence"] = new.dosing_sequence
new["Antibiotic_Name"] = new.antibiotic
new["study"] = "NEW:" + new.screen_id

COLS = ["study", "Organism", "Antibiotic_Class", "Antibiotic_Mechanism",
        "Bactericidal_Static", "Assay", "Dosing_Sequence", "Antibiotic_Name",
        "Synergy_Outcome", "Interaction_Classification"]
comb = pd.concat([base[COLS], new[COLS]], ignore_index=True)
comb.to_csv(DATA / "priority4_combined_dataset.csv", index=False)

OUT["n_published"] = int(len(base))
OUT["n_new"] = int(len(new))
OUT["n_combined"] = int(len(comb))
OUT["studies_published"] = int(base.study.nunique())
OUT["studies_new"] = int(new.screen_id.nunique())
OUT["studies_combined"] = int(comb.study.nunique())


# --------------------------------------------------------------------------
# 2. Estimators
# --------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h) * 100, min(1.0, c + h) * 100)


def icc1(df):
    """One-way ANOVA ICC(1) - the estimator the manuscript now uses."""
    g = df.groupby("study").Synergy_Outcome.agg(["sum", "size"])
    g["p"] = g["sum"] / g["size"]
    N, K = len(df), len(g)
    pbar = df.Synergy_Outcome.mean()
    ssb = float(np.sum(g["size"] * (g["p"] - pbar) ** 2))
    ssw = float(np.sum((df.Synergy_Outcome - df.study.map(g["p"])) ** 2))
    msb, msw = ssb / (K - 1), ssw / (N - K)
    n0 = (N - np.sum(g["size"] ** 2) / N) / (K - 1)
    icc = (msb - msw) / (msb + (n0 - 1) * msw)
    return dict(icc=icc, N=N, K=K, mean_cluster=N / K,
                deff=1 + (N / K - 1) * icc, eff_n=(N / (1 + (N / K - 1) * icc)),
                sizes=g["size"].tolist())


def meta2(df, label):
    N = len(df)
    k = int(df.Synergy_Outcome.sum())
    a = int((df.Interaction_Classification == "Antagonism").sum())
    lo, hi = wilson(k, N)
    alo, ahi = wilson(a, N)
    ic = icc1(df)
    d = dict(label=label, n_determinate=N, n_studies=ic["K"],
             synergy_n=k, synergy_pct=round(k / N * 100, 2),
             synergy_ci=[round(lo, 1), round(hi, 1)],
             antagonism_n=a, antagonism_pct=round(a / N * 100, 2),
             antagonism_ci=[round(alo, 1), round(ahi, 1)],
             icc1=round(ic["icc"], 3), mean_cluster=round(ic["mean_cluster"], 2),
             deff=round(ic["deff"], 2), eff_n=round(ic["eff_n"]))
    # cluster-adjusted interval on the ICC(1) deff
    p = k / N
    se = np.sqrt(p * (1 - p) / N) * np.sqrt(ic["deff"])
    d["cluster_ci"] = [round(max(0, p - 1.96 * se) * 100, 1),
                       round(min(1, p + 1.96 * se) * 100, 1)]
    # cluster bootstrap
    rng = np.random.default_rng(42)
    cl = list(zip(ic["sizes"], df.groupby("study").Synergy_Outcome.sum().tolist()))
    bs = []
    for _ in range(10000):
        idx = rng.integers(0, len(cl), len(cl))
        kk = sum(cl[i][1] for i in idx)
        nn = sum(cl[i][0] for i in idx)
        bs.append(kk / nn)
    blo, bhi = np.percentile(bs, [2.5, 97.5])
    d["bootstrap_ci"] = [round(blo * 100, 1), round(bhi * 100, 1)]
    # permutation test for clustering
    y = df.Synergy_Outcome.values.astype(float)
    gi = pd.factorize(df.study)[0]
    G = gi.max() + 1
    sizes = np.bincount(gi, minlength=G).astype(float)

    def eta2(yv):
        pb = yv.mean()
        m = np.bincount(gi, weights=yv, minlength=G) / sizes
        return np.sum(sizes * (m - pb) ** 2) / np.sum((yv - pb) ** 2)

    obs = eta2(y)
    perm = np.array([eta2(rng.permutation(y)) for _ in range(20000)])
    d["eta2"] = round(float(obs), 4)
    d["eta2_null_mean"] = round(float(perm.mean()), 4)
    d["perm_p"] = round(float((perm >= obs).mean()), 4)
    d["size_range"] = [int(min(ic["sizes"])), int(max(ic["sizes"]))]
    d["singletons"] = int(sum(1 for s in ic["sizes"] if s == 1))
    return d


OUT["pooled_before"] = meta2(base, "published")
OUT["pooled_after"] = meta2(comb, "combined")

# --------------------------------------------------------------------------
# 3. By antibiotic class and mechanism
# --------------------------------------------------------------------------
def by(df, col, min_n=1):
    rows = []
    for v in df[col].dropna().unique():
        s = df[df[col] == v]
        k, n = int(s.Synergy_Outcome.sum()), len(s)
        if n < min_n:
            continue
        lo, hi = wilson(k, n)
        rows.append(dict(category=v, synergy_n=k, total_n=n,
                         synergy_pct=round(k / n * 100, 2),
                         ci=[round(lo, 1), round(hi, 1)],
                         antagonism_n=int((s.Interaction_Classification == "Antagonism").sum())))
    return sorted(rows, key=lambda r: -r["synergy_pct"])


OUT["by_class_after"] = by(comb, "Antibiotic_Class")
OUT["by_mechanism_after"] = by(comb, "Antibiotic_Mechanism")
OUT["by_class_before"] = by(base, "Antibiotic_Class")
OUT["by_mechanism_before"] = by(base, "Antibiotic_Mechanism")

# --------------------------------------------------------------------------
# 4. By organism
# --------------------------------------------------------------------------
rows = []
for o in comb.Organism.dropna().unique():
    s = comb[comb.Organism == o]
    k, n = int(s.Synergy_Outcome.sum()), len(s)
    lo, hi = wilson(k, n)
    rows.append(dict(organism=o, synergy_n=k, total_n=n, studies=int(s.study.nunique()),
                     synergy_pct=round(k / n * 100, 2), ci=[round(lo, 1), round(hi, 1)],
                     antagonism_n=int((s.Interaction_Classification == "Antagonism").sum())))
OUT["by_organism_after"] = sorted(rows, key=lambda r: -r["total_n"])

# --------------------------------------------------------------------------
# 5. Design-based and outcome-selected subsets
# --------------------------------------------------------------------------
def subsets(df):
    g = df.groupby("study").agg(nmech=("Antibiotic_Mechanism", "nunique"),
                                narm=("Synergy_Outcome", "size"),
                                ksyn=("Synergy_Outcome", "sum"))
    design = set(g[g.nmech >= 2].index)
    disc = set(g[(g.narm >= 2) & (g.ksyn > 0) & (g.ksyn < g.narm)].index)
    out = {}
    for name, ids in [("design_based", design), ("outcome_selected", disc)]:
        sub = df[df.study.isin(ids)]
        k, n = int(sub.Synergy_Outcome.sum()), len(sub)
        per_study = sub.groupby("study").Synergy_Outcome.mean()
        out[name] = dict(studies=int(sub.study.nunique()), arms=n,
                         synergy_arms_pct=round(k / n * 100, 2),
                         synergy_studylevel_pct=round(per_study.mean() * 100, 2))
    return out


OUT["subsets_before"] = subsets(base)
OUT["subsets_after"] = subsets(comb)

# --------------------------------------------------------------------------
# 6. Class contrasts (FQ vs beta-lactam) under each subset
# --------------------------------------------------------------------------
FQ, BL = "DNA gyrase/topo IV", "Cell wall (PBP)"
def contrast(df, label):
    gm = df.groupby("study").Antibiotic_Mechanism.nunique()
    design = set(gm[gm >= 2].index)
    g2 = df.groupby("study").agg(narm=("Synergy_Outcome", "size"), ksyn=("Synergy_Outcome", "sum"))
    disc = set(g2[(g2.narm >= 2) & (g2.ksyn > 0) & (g2.ksyn < g2.narm)].index)
    r = {}
    for name, ids in [("pooled", None), ("outcome_selected", disc), ("design_based", design)]:
        sub = df if ids is None else df[df.study.isin(ids)]
        f = sub[sub.Antibiotic_Mechanism == FQ]
        b = sub[sub.Antibiotic_Mechanism == BL]
        fp = f.Synergy_Outcome.mean() * 100 if len(f) else float("nan")
        bp = b.Synergy_Outcome.mean() * 100 if len(b) else float("nan")
        r[name] = dict(fq_pct=round(fp, 2), fq_n=len(f),
                       bl_pct=round(bp, 2), bl_n=len(b),
                       gap_pp=round(fp - bp, 2))
    return r

OUT["contrast_before"] = contrast(base, "before")
OUT["contrast_after"] = contrast(comb, "after")

# --------------------------------------------------------------------------
# 7. Head-to-head paired analysis (FQ vs BL within the same study)
# --------------------------------------------------------------------------
def paired(df):
    rows = []
    for st, s in df.groupby("study"):
        f = s[s.Antibiotic_Mechanism == FQ]
        b = s[s.Antibiotic_Mechanism == BL]
        if len(f) and len(b):
            rows.append(dict(study=st, gyr=f.Synergy_Outcome.mean() * 100,
                             pbp=b.Synergy_Outcome.mean() * 100))
    if len(rows) < 3:
        return None
    pdf = pd.DataFrame(rows)
    b = int(((pdf.gyr == 100) & (pdf.pbp < 100)).sum())
    c = int(((pdf.gyr < 100) & (pdf.pbp == 100)).sum())
    both = int(((pdf.gyr == 100) & (pdf.pbp == 100)).sum())
    neither = int(((pdf.gyr < 100) & (pdf.pbp < 100)).sum())
    # exact McNemar
    from scipy.stats import binomtest
    p = binomtest(b, b + c, 0.5).pvalue if (b + c) > 0 else float("nan")
    return dict(n_studies=len(pdf), b=b, c=c, both100=both, neither100=neither,
                mcnemar_exact_p=round(float(p), 4),
                mean_fq=round(pdf.gyr.mean(), 1), mean_bl=round(pdf.pbp.mean(), 1),
                diff_pp=round(pdf.gyr.mean() - pdf.pbp.mean(), 1),
                data=pdf.to_dict("records"))

OUT["paired_before"] = paired(base)
OUT["paired_after"] = paired(comb)

# --------------------------------------------------------------------------
# 8. Cidal vs static, envelope-active, assay mix
# --------------------------------------------------------------------------
def fisher(df, mask, label):
    a = df[mask]
    b = df[~mask]
    ka, na = int(a.Synergy_Outcome.sum()), len(a)
    kb, nb = int(b.Synergy_Outcome.sum()), len(b)
    if na == 0 or nb == 0:
        return None
    tab = [[ka, na - ka], [kb, nb - kb]]
    _, p = stats.fisher_exact(tab)
    return dict(label=label, group_n=na, group_pct=round(ka / na * 100, 2),
                rest_n=nb, rest_pct=round(kb / nb * 100, 2), fisher_p=round(float(p), 3))

for tag, df in [("before", base), ("after", comb)]:
    OUT[f"cidal_static_{tag}"] = fisher(df, df.Bactericidal_Static == "cidal", "cidal vs static")
    env = df.Antibiotic_Mechanism.isin(
        ["Cell wall (PBP)", "Cell wall (D-Ala-D-Ala)", "Cell wall (MurA)", "Membrane", "Membrane (LPS)"])
    OUT[f"envelope_{tag}"] = fisher(df, env, "envelope-active vs other")

def assay_mix(df):
    first = df.Assay.astype(str).str.split(";").str[0].str.strip()
    rows = []
    for a in ["Checkerboard", "Biofilm", "Time-Kill", "In Vivo"]:
        n = int((first == a).sum())
        k = int(df[first == a].Synergy_Outcome.sum())
        rows.append(dict(assay=a, n=n, synergy_n=k,
                         synergy_pct=round(k / n * 100, 2) if n else None))
    return rows

OUT["assay_before"] = assay_mix(base)
OUT["assay_after"] = assay_mix(comb)

# --------------------------------------------------------------------------
# 9. PRISMA-ScR flow (updated)
# --------------------------------------------------------------------------
OUT["flow_after"] = dict(
    records_identified=4402,
    excluded_at_triage_score_lt_14=None,  # depends on adopted screen
    prioritised_published=254,
    prioritised_revised=282,
    additional_screened=28,
    additional_included=9,
    fulltext_parsed_published=210,
    fulltext_parsed_revised=210 + 21,
    abstract_only_published=44,
    abstract_only_revised=44 + 7,
    excluded_at_eligibility_published=97,
    excluded_at_eligibility_revised=97 + 19,
    studies_included_published=157,
    studies_included_revised=157 + 9,
    observations_published=419,
    observations_revised=419 + 23,
    determinate_published=272,
    determinate_revised=272 + 23,
)

# --------------------------------------------------------------------------
# 10. Sensitivity ladder (unchanged values, restated for the SI)
# --------------------------------------------------------------------------
OUT["threshold_ladder"] = {">=20": 44, ">=18": 83, ">=16": 151,
                           ">=14": 254, ">=12": 403, ">=10": 607}
OUT["threshold_ladder_revised"] = {">=20": 48, ">=18": 93, ">=16": 167,
                                   ">=14": 282, ">=12": 440, ">=10": 660}

(DATA / "reanalysis_results.json").write_text(json.dumps(OUT, indent=1))

# --------------------------------------------------------------------------
# console summary
# --------------------------------------------------------------------------
def show(tag, d):
    print(f"\n[{tag}]  n={d['n_determinate']}  studies={d['n_studies']}")
    print(f"   synergy {d['synergy_n']}/{d['n_determinate']} = {d['synergy_pct']}% "
          f"[{d['synergy_ci'][0]}-{d['synergy_ci'][1]}]")
    print(f"   antagonism {d['antagonism_n']} = {d['antagonism_pct']}% "
          f"[{d['antagonism_ci'][0]}-{d['antagonism_ci'][1]}]")
    print(f"   ICC(1)={d['icc1']}  mean cluster={d['mean_cluster']}  deff={d['deff']}  eff n={d['eff_n']}")
    print(f"   cluster CI {d['cluster_ci']}   bootstrap {d['bootstrap_ci']}")
    print(f"   eta2={d['eta2']} (null {d['eta2_null_mean']}) perm P={d['perm_p']}")
    print(f"   cluster sizes {d['size_range']}, singletons {d['singletons']}")

show("BEFORE", OUT["pooled_before"])
show("AFTER ", OUT["pooled_after"])
print("\n[CONTRASTS]")
for k in ["pooled", "outcome_selected", "design_based"]:
    b, a = OUT["contrast_before"][k], OUT["contrast_after"][k]
    print(f"  {k:18s} before gap {b['gap_pp']:5.1f} pp -> after {a['gap_pp']:5.1f} pp")
print("\n[PAIRED]")
for t in ["paired_before", "paired_after"]:
    p = OUT[t]
    print(f"  {t}: n={p['n_studies']} b={p['b']} c={p['c']} P={p['mcnemar_exact_p']} "
          f"FQ {p['mean_fq']}% BL {p['mean_bl']}%")
print("\nFLOW:", json.dumps(OUT["flow_after"], indent=1))
print(f"\nwrote {DATA / 'reanalysis_results.json'}")
