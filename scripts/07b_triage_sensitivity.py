#!/usr/bin/env python3
"""
Step 7b - ESKAPE triage, SENSITIVITY RE-RUN (Priority 4).

Purpose
-------
The published triage (scripts/07_eskape_triage.py) applied -3 per off-target
term, where the off-target list included the environmental isolation-source
terms 'wastewater' and 'sewage'. Phage provenance sentences ("phage X was
isolated from hospital wastewater") are routine in abstracts of on-topic
papers, so this penalty fired on a methods sentence rather than on the topic.

This script:
  1. reproduces the published scoring EXACTLY (asserted, 4402/4402);
  2. re-scores with 'wastewater' and 'sewage' removed from the off-target list;
  3. re-scores with the duplicate-term bug also fixed (case-insensitive dedupe),
     which the published rule ("per distinct term") already implied;
  4. reports the threshold ladder, the newly-crossing records, and whether any
     of them are already in the included set.

It writes no production file. Output paths are under data/ with a _sens suffix
so the published artefacts are never silently overwritten.

Usage:  python scripts/07b_triage_sensitivity.py
"""
import json
import re
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"

# --- constants copied verbatim from 07_eskape_triage.py -------------------
ORG_PAT = {
    "Enterococcus faecium": r"enterococcus faecium|e\.\s*faecium|\bVRE\b|vancomycin[- ]resistant enterococc",
    "Staphylococcus aureus": r"staphylococcus aureus|s\.\s*aureus|\bMRSA\b|methicillin[- ]resistant",
    "Klebsiella pneumoniae": r"klebsiella pneumoniae|k\.\s*pneumoniae|\bCRKP\b",
    "Acinetobacter baumannii": r"acinetobacter baumannii|a\.\s*baumannii|\bCRAB\b|a\.\s*baumannii",
    "Pseudomonas aeruginosa": r"pseudomonas aeruginosa|p\.\s*aeruginosa|\bMDRPA\b",
    "Enterobacter spp.": r"enterobacter cloacae|enterobacter aerogenes|e\.\s*cloacae|enterobacter",
}

REVIEW = re.compile(
    r"\b(review|meta-analysis|systematic review|perspective|overview|advances in|"
    r"current landscape|challenges and opportunities|editorial|commentary|opinion|"
    r"state of the art|recent advances|an update|future directions|landscape|"
    r"bibliometric|treatment strategies|clinical challenges)\b", re.I)

EXPERIMENT = re.compile(
    r"\b(checkerboard|fractional inhibitory|\bFIC\b|time[\s\-]?kill|kill(?:ing)? curve|"
    r"combination (?:therapy|treatment|assay|index)|synerg\w* (?:effect|activity|interaction)|"
    r"phage[- ]antibiotic synergy|\bPAS\b|combined (?:treatment|therapy)|co-?admin\w*|"
    r"in combination with|biofilm|Galleria|mellonella|murine|mice|in vivo)\b", re.I)

ABX = re.compile(
    r"\b(meropenem|imipenem|ertapenem|doripenem|colistin|polymyxin|amikacin|gentamicin|"
    r"tobramycin|ciprofloxacin|levofloxacin|moxifloxacin|ceftazidime|cefepime|ceftriaxone|"
    r"cefotaxime|ceftaroline|ceftobiprole|avibactam|tazobactam|sulbactam|tigecycline|"
    r"tetracycline|minocycline|azithromycin|erythromycin|clindamycin|linezolid|"
    r"daptomycin|vancomycin|ampicillin|amoxicillin|oxacillin|piperacillin|aztreonam|"
    r"rifampicin|rifampin|fosfomycin|chloramphenicol|nitrofurantoin|cefiderocol|"
    r"ceftolozane|dalbavancin|oritavancin|teicoplanin|kanamycin|streptomycin|"
    r"antibiotic|antimicrobial)s?\b", re.I)

# PUBLISHED off-target list (includes environmental isolation-source terms)
OTHER_PUBLISHED = re.compile(
    r"\b(escherichia coli|salmonella|mycobacter|clostridi|streptococc|listeria|"
    r"vibrio|cronobacter|food|poultry|cattle|plant|crop|wastewater|sewage|"
    r"environmental sample|biocontrol|biosensor|genome announcement|"
    r"crystal structure|structural basis)\b", re.I)

# REVISED off-target list: 'wastewater' and 'sewage' removed
OTHER_REVISED = re.compile(
    r"\b(escherichia coli|salmonella|mycobacter|clostridi|streptococc|listeria|"
    r"vibrio|cronobacter|food|poultry|cattle|plant|crop|"
    r"environmental sample|biocontrol|biosensor|genome announcement|"
    r"crystal structure|structural basis)\b", re.I)

QUANT = re.compile(
    r"\b(\d+(?:\.\d+)?\s*(?:log|log10|log\s*10|%|fold)|\d+\s*/\s*\d+\s*(?:isolates|strains)|"
    r"MIC|FIC|mg/L|CFU)\b", re.I)


def score(r, other_pat, dedupe_ci=False):
    """Scoring rule. dedupe_ci=True applies 'per distinct term' case-insensitively."""
    text = f"{r['title']} {r['abstract']}"
    if not text.strip():
        return -99, []
    s, why = 0, []
    if not re.search(r"phage|bacteriophage", text, re.I):
        return -99, ["no phage"]
    s += 2
    if REVIEW.search(r["title"]):
        s -= 8
        why.append("review title")
    n_org = len({k for k, p in ORG_PAT.items() if re.search(p, text, re.I)})
    if n_org == 0:
        s -= 6
        why.append("no ESKAPE organism")
    else:
        s += 4
        why.append(f"{n_org} organism(s)")
    n_exp = len(EXPERIMENT.findall(text))
    s += min(n_exp, 6)
    n_abx = len(set(m.lower() for m in ABX.findall(text)))
    s += min(n_abx, 4)
    if QUANT.search(text):
        s += 3
        why.append("quantitative")
    if re.search(r"synerg", text, re.I):
        s += 3
        why.append("synergy term")
    s += 2 if r.get("pmcid") else -2
    o = other_pat.findall(text)
    if o:
        uniq = {m.lower() for m in o} if dedupe_ci else set(o)
        s -= 3 * len(uniq)
        why.append("off-target:" + ",".join(sorted(uniq)[:2]))
    return s, why


def main():
    pool_path = DATA / "eskape_pool.json"
    if not pool_path.exists():
        print("[SKIP] data/eskape_pool.json is not present in this repository.")
        print("       It is a screening-stage artefact regenerable with")
        print("       scripts/06_eskape_search.py (requires network access to Europe PMC).")
        print("       The re-scored output is committed as data/triage_rescore_4402.csv,")
        print("       and the 28 newly-crossing records as")
        print("       data/priority4_screening_log_28.csv.")
        return
    pool = json.loads(pool_path.read_text())

    # ---- 0. reproduce the published result exactly ----
    published = {r["pmid"] + "|" + (r["doi"] or ""): r["_score"]
                 for r in json.loads((DATA / "eskape_triage.json").read_text())}
    mism = 0
    for r in pool:
        s, _ = score(r, OTHER_PUBLISHED, dedupe_ci=False)
        if s != published.get(r["pmid"] + "|" + (r["doi"] or "")):
            mism += 1
    print(f"[0] reproduction of published scoring: "
          f"{len(pool) - mism}/{len(pool)} exact, {mism} mismatched")
    assert mism == 0, "baseline does not reproduce; aborting before any change"

    # ---- scenarios ----
    scen = {
        "A published": (OTHER_PUBLISHED, False),
        "B published + distinct-term fix": (OTHER_PUBLISHED, True),
        "C isolation terms removed": (OTHER_REVISED, False),
        "D C + distinct-term fix": (OTHER_REVISED, True),
    }
    print(f"\n{'scenario':34s}{'>=18':>6}{'>=16':>6}{'>=14':>6}{'>=12':>6}{'>=10':>6}")
    ladders = {}
    for name, (pat, dd) in scen.items():
        s = [score(r, pat, dd)[0] for r in pool]
        ladders[name] = s
        print(f"{name:34s}{sum(x>=18 for x in s):6d}{sum(x>=16 for x in s):6d}"
              f"{sum(x>=14 for x in s):6d}{sum(x>=12 for x in s):6d}{sum(x>=10 for x in s):6d}")

    # ---- which records newly cross, under the RECOMMENDED scenario D ----
    base, fixed = ladders["A published"], ladders["D C + distinct-term fix"]
    newly = [(fixed[i], base[i], r) for i, r in enumerate(pool)
             if base[i] < 14 <= fixed[i]]
    newly.sort(key=lambda x: -x[0])

    included = {s.lower() for s in
                __import__("pandas").read_csv(DATA.parent / "eskape_pas_dataset_clean.csv")
                ["PMID_or_DOI"].astype(str)}
    ex = __import__("pandas").read_csv(DATA.parent / "github_repo/data/eligibility_exclusions_97.csv")
    excluded = {str(k).lower().replace("/", "_") for k in ex.iloc[:, 0]}

    def already(rec):
        keys = {str(rec.get("pmid", "")).lower(), str(rec.get("doi", "")).lower()}
        if keys & included:
            return "INCLUDED"
        if {k.replace("/", "_") for k in keys} & excluded:
            return "EXCLUDED"
        return "NEW-SCREEN"

    print(f"\nRecords newly crossing the >=14 threshold: {len(newly)}")
    from collections import Counter
    stat = Counter(already(r) for _, _, r in newly)
    print(f"  status: {dict(stat)}")

    out = []
    for new_s, old_s, r in newly:
        out.append({
            "score_published": old_s,
            "score_revised": new_s,
            "delta": new_s - old_s,
            "screening_status": already(r),
            "pmid": r.get("pmid") or "",
            "doi": r.get("doi") or "",
            "year": r.get("year"),
            "journal": r.get("journal") or "",
            "title": r.get("title") or "",
        })
    (DATA / "triage_priority4_newly_crossing.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA / 'triage_priority4_newly_crossing.json'}")

    # ---- revised full triage, written under a sensitivity name ----
    scored = []
    for r in pool:
        s, why = score(r, OTHER_REVISED, dedupe_ci=True)
        rr = dict(r)
        rr["_score"] = s
        rr["_why"] = why
        scored.append(rr)
    scored.sort(key=lambda x: -x["_score"])
    # Write the redistributable form only: identifiers, scores and rationale.
    # Abstracts are third-party copyright and are deliberately not reproduced.
    import csv as _csv
    rows = [{"pmid": r.get("pmid", ""), "doi": r.get("doi", ""),
             "year": r.get("year", ""), "journal": r.get("journal", ""),
             "pmcid": r.get("pmcid", ""), "score": r["_score"],
             "why": "; ".join(r.get("_why") or []),
             "organisms": "; ".join(r.get("_organisms") or [])} for r in scored]
    out_csv = DATA / "triage_rescore_4402.csv"
    with open(out_csv, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv} ({len(rows)} rows, no article text) - published file untouched")

    n_a = sum(1 for r in pool if score(r, OTHER_PUBLISHED, False)[0] >= 14)
    n_d = sum(1 for r in pool if score(r, OTHER_REVISED, True)[0] >= 14)
    print(f"\nHEADLINE: prioritised set {n_a} -> {n_d} "
          f"(+{len(newly)} to screen). Published 254 is unchanged.")


if __name__ == "__main__":
    main()
