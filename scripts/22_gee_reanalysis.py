#!/usr/bin/env python3
"""
Step 22 - GEE for the fluoroquinolone vs beta-lactam head-to-head contrast.

This refits the cluster-robust generalised estimating equation on the
reanalysed dataset. It exists because the McNemar test was re-run after the
sensitivity screen but the GEE had not been, so the previously reported odds
ratio could not simply be carried forward.

Specification (unchanged from the manuscript):
  response  : y  = Synergy_Outcome (1 = synergy, 0 = not)
  predictor : is_gyr = 1 if mechanism is DNA gyrase/topo IV, 0 if Cell wall (PBP)
  cluster   : study  (PMID or DOI)
  family    : Binomial
  link      : logit
  working correlation : exchangeable
  variance  : robust (sandwich)

IMPORTANT - the fitting sample. This is fitted on the DESIGN-BASED subset
(studies testing >= 2 distinct antibiotic mechanisms), restricted to the two
mechanism classes of interest, NOT on the 26 paired studies. That is how the
published estimate was produced. Fitting on the paired subset instead gives a
different point estimate (~3.9 vs ~4.0) and a much narrower interval; switching
estimands silently would break comparability with the previous version.
See docs/KNOWN_DISCREPANCIES.md section 7.

Usage:  python scripts/22_gee_reanalysis.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

FQ, BL = "DNA gyrase/topo IV", "Cell wall (PBP)"


def fit(df, label):
    g = df.groupby("study").Antibiotic_Mechanism.nunique()
    dsub = df[df.study.isin(g[g >= 2].index)]
    sub = dsub[dsub.Antibiotic_Mechanism.isin([BL, FQ])].copy()
    sub["is_gyr"] = (sub.Antibiotic_Mechanism == FQ).astype(int)
    sub["y"] = sub.Synergy_Outcome.astype(int)

    m = smf.gee("y ~ is_gyr", groups="study", data=sub,
                family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = np.exp(m.conf_int().loc["is_gyr"])
    out = dict(label=label,
               design_based_studies=int(dsub.study.nunique()),
               gee_arms=int(len(sub)),
               OR=round(float(np.exp(m.params["is_gyr"])), 2),
               ci=[round(float(ci[0]), 2), round(float(ci[1]), 2)],
               p=round(float(m.pvalues["is_gyr"]), 3),
               fq_arms=int((sub.is_gyr == 1).sum()),
               bl_arms=int((sub.is_gyr == 0).sum()))
    print(f"{label}: design-based studies={out['design_based_studies']}, "
          f"GEE arms={out['gee_arms']} (FQ {out['fq_arms']}, BL {out['bl_arms']})")
    print(f"  OR = {out['OR']} (95% CI {out['ci'][0]}-{out['ci'][1]}), P = {out['p']}")
    return out


def main():
    raw = pd.read_csv(DATA / "eskape_pas_dataset_clean.csv")
    base = raw[raw.Synergy_Outcome.notna()].copy()
    base["study"] = "OLD:" + base.PMID_or_DOI.astype(str)

    comb = pd.read_csv(DATA / "priority4_combined_dataset.csv")

    print("Refitting GEE (exchangeable, logit, robust sandwich)\n")
    pub = fit(base, "published (157-study)")
    new = fit(comb, "reanalysed (166-study)")

    print("\nPublished manuscript reports: OR 4.6 (95% CI 0.45-47.3), P = 0.20")
    print(f"Replication above            : OR {pub['OR']} "
          f"(95% CI {pub['ci'][0]}-{pub['ci'][1]}), P = {pub['p']}")
    ok = (abs(pub["OR"] - 4.6) < 0.1 and abs(pub["ci"][0] - 0.45) < 0.02
          and abs(pub["ci"][1] - 47.3) < 0.2)
    print(f"Replication {'MATCHES' if ok else 'DOES NOT MATCH'} the published estimate")

    res = dict(published=pub, reanalysed=new)
    (DATA / "gee_results.json").write_text(json.dumps(res, indent=1))

    # merge into the master reanalysis file so one file carries everything
    rp = DATA / "reanalysis_results.json"
    if rp.exists():
        master = json.loads(rp.read_text())
        master["gee_before"] = pub
        master["gee_after"] = new
        rp.write_text(json.dumps(master, indent=1))
        print(f"\nmerged into {rp.name}")

    print(f"\nwrote {DATA / 'gee_results.json'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
