#!/usr/bin/env python3
"""
Step 7 - ESKAPE triage: keep records that plausibly report an original
phage-antibiotic combination experiment with a measurable interaction outcome.
Assigns organism + study design signals. Writes data/eskape_triage.json.
"""
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
pool = json.loads((DATA / "eskape_pool.json").read_text())

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

# evidence that an actual combination experiment was done
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

OTHER = re.compile(
    r"\b(escherichia coli|salmonella|mycobacter|clostridi|streptococc|listeria|"
    r"vibrio|cronobacter|food|poultry|cattle|plant|crop|wastewater|sewage|"
    r"environmental sample|biocontrol|biosensor|genome announcement|"
    r"crystal structure|structural basis)\b", re.I)

QUANT = re.compile(
    r"\b(\d+(?:\.\d+)?\s*(?:log|log10|log\s*10|%|fold)|\d+\s*/\s*\d+\s*(?:isolates|strains)|"
    r"MIC|FIC|mg/L|CFU)\b", re.I)


def score(r):
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
    if n_exp == 0:
        why.append("no experiment language")
    n_abx = len(set(m.lower() for m in ABX.findall(text)))
    s += min(n_abx, 4)
    if n_abx == 0:
        why.append("no antibiotic named")
    if QUANT.search(text):
        s += 3
        why.append("quantitative")
    if re.search(r"synerg", text, re.I):
        s += 3
        why.append("synergy term")
    if r.get("pmcid"):
        s += 2
    else:
        s -= 2
    o = OTHER.findall(text)
    if o:
        s -= 3 * len(set(o))
        why.append("off-target:" + ",".join(sorted(set(o))[:2]))
    return s, why


scored = []
for r in pool:
    s, why = score(r)
    r["_score"] = s
    r["_why"] = why
    # assign primary organism by most specific match
    txt = f"{r['title']} {r['abstract']}".lower()
    hits = [k for k, p in ORG_PAT.items() if re.search(p, r["title"] + " " + r["abstract"], re.I)]
    r["_organisms"] = hits or r["organisms"]
    scored.append(r)

scored.sort(key=lambda x: -x["_score"])
for thr in (18, 16, 14, 12, 10):
    print(f"score >= {thr}: {sum(1 for r in scored if r['_score'] >= thr)}")

(DATA / "eskape_triage.json").write_text(json.dumps(scored, indent=1))
sel = [r for r in scored if r["_score"] >= 14]
print(f"\nselected for full-text fetch: {len(sel)}")
print("\nBy organism:")
for k, v in Counter(o for r in sel for o in r["_organisms"][:1]).most_common():
    print(f"  {k:26s} {v}")
print("\nBy year:", dict(sorted(Counter(r["year"] for r in sel).items())))
