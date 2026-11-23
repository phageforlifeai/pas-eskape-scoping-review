#!/usr/bin/env python3
"""
Step 9 - ESKAPE-wide structured extraction of phage-antibiotic interaction data.

Reuses the validated CRKP extraction logic (clause-level drug attribution,
negation handling, per-arm quantification) and adds:
  * organism assignment (ESKAPE member)
  * phage family / genus
  * drug mechanism class (for the "which antibiotic" question)
  * resistance genotype (carbapenemase / MRSA / VRE / ESBL)
  * a mechanistically-informed feature set for predicting synergy

Output: eskape_pas_dataset.csv (one row per study x antibiotic arm)
"""
import csv
import json
import re
from collections import Counter
from pathlib import Path

DATA = Path("/home/user/data")
OUT = Path("/home/user")
corpus = json.loads((DATA / "eskape_corpus.json").read_text())

# ------------------------------- vocabularies -------------------------------
ORG_PAT = {
    "E. faecium": r"enterococcus faecium|e\.\s*faecium|\bVRE\b|vancomycin[- ]resistant enterococc",
    "S. aureus": r"staphylococcus aureus|s\.\s*aureus|\bMRSA\b",
    "K. pneumoniae": r"klebsiella pneumoniae|k\.\s*pneumoniae|\bCRKP\b",
    "A. baumannii": r"acinetobacter baumannii|a\.\s*baumannii|\bCRAB\b",
    "P. aeruginosa": r"pseudomonas aeruginosa|p\.\s*aeruginosa",
    "Enterobacter spp.": r"enterobacter cloacae|enterobacter aerogenes|e\.\s*cloacae",
}

# antibiotic -> (canonical, class, mechanism, bacteriostatic?)
ABX = {
    "meropenem": ("Meropenem", "Carbapenem", "Cell wall (PBP)", "cidal"),
    "imipenem": ("Imipenem", "Carbapenem", "Cell wall (PBP)", "cidal"),
    "ertapenem": ("Ertapenem", "Carbapenem", "Cell wall (PBP)", "cidal"),
    "doripenem": ("Doripenem", "Carbapenem", "Cell wall (PBP)", "cidal"),
    "biapenem": ("Biapenem", "Carbapenem", "Cell wall (PBP)", "cidal"),
    "ceftazidime/avibactam": ("Ceftazidime-avibactam", "Cephalosporin/BLI", "Cell wall (PBP)", "cidal"),
    "ceftazidime-avibactam": ("Ceftazidime-avibactam", "Cephalosporin/BLI", "Cell wall (PBP)", "cidal"),
    "ceftazidime avibactam": ("Ceftazidime-avibactam", "Cephalosporin/BLI", "Cell wall (PBP)", "cidal"),
    "ceftolozane/tazobactam": ("Ceftolozane-tazobactam", "Cephalosporin/BLI", "Cell wall (PBP)", "cidal"),
    "ceftolozane-tazobactam": ("Ceftolozane-tazobactam", "Cephalosporin/BLI", "Cell wall (PBP)", "cidal"),
    "aztreonam/avibactam": ("Aztreonam-avibactam", "Monobactam/BLI", "Cell wall (PBP)", "cidal"),
    "ampicillin/sulbactam": ("Ampicillin-sulbactam", "Penicillin/BLI", "Cell wall (PBP)", "cidal"),
    "piperacillin/tazobactam": ("Piperacillin-tazobactam", "Penicillin/BLI", "Cell wall (PBP)", "cidal"),
    "ceftazidime": ("Ceftazidime", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "cefepime": ("Cefepime", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "ceftriaxone": ("Ceftriaxone", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "cefotaxime": ("Cefotaxime", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "ceftaroline": ("Ceftaroline", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "ceftobiprole": ("Ceftobiprole", "Cephalosporin", "Cell wall (PBP)", "cidal"),
    "cefiderocol": ("Cefiderocol", "Siderophore cephalosporin", "Cell wall (PBP)", "cidal"),
    "colistin": ("Colistin", "Polymyxin", "Membrane (LPS)", "cidal"),
    "polymyxin b": ("Polymyxin B", "Polymyxin", "Membrane (LPS)", "cidal"),
    "polymyxin": ("Polymyxin B", "Polymyxin", "Membrane (LPS)", "cidal"),
    "amikacin": ("Amikacin", "Aminoglycoside", "Protein synthesis (30S)", "cidal"),
    "gentamicin": ("Gentamicin", "Aminoglycoside", "Protein synthesis (30S)", "cidal"),
    "tobramycin": ("Tobramycin", "Aminoglycoside", "Protein synthesis (30S)", "cidal"),
    "kanamycin": ("Kanamycin", "Aminoglycoside", "Protein synthesis (30S)", "cidal"),
    "streptomycin": ("Streptomycin", "Aminoglycoside", "Protein synthesis (30S)", "cidal"),
    "ciprofloxacin": ("Ciprofloxacin", "Fluoroquinolone", "DNA gyrase/topo IV", "cidal"),
    "levofloxacin": ("Levofloxacin", "Fluoroquinolone", "DNA gyrase/topo IV", "cidal"),
    "moxifloxacin": ("Moxifloxacin", "Fluoroquinolone", "DNA gyrase/topo IV", "cidal"),
    "tigecycline": ("Tigecycline", "Tetracycline", "Protein synthesis (30S)", "static"),
    "tetracycline": ("Tetracycline", "Tetracycline", "Protein synthesis (30S)", "static"),
    "minocycline": ("Minocycline", "Tetracycline", "Protein synthesis (30S)", "static"),
    "doxycycline": ("Doxycycline", "Tetracycline", "Protein synthesis (30S)", "static"),
    "azithromycin": ("Azithromycin", "Macrolide", "Protein synthesis (50S)", "static"),
    "erythromycin": ("Erythromycin", "Macrolide", "Protein synthesis (50S)", "static"),
    "clindamycin": ("Clindamycin", "Lincosamide", "Protein synthesis (50S)", "static"),
    "linezolid": ("Linezolid", "Oxazolidinone", "Protein synthesis (50S)", "static"),
    "chloramphenicol": ("Chloramphenicol", "Amphenicol", "Protein synthesis (50S)", "static"),
    "rifampicin": ("Rifampicin", "Rifamycin", "RNA polymerase", "cidal"),
    "rifampin": ("Rifampicin", "Rifamycin", "RNA polymerase", "cidal"),
    "fosfomycin": ("Fosfomycin", "Fosfomycin", "Cell wall (MurA)", "cidal"),
    "vancomycin": ("Vancomycin", "Glycopeptide", "Cell wall (D-Ala-D-Ala)", "cidal"),
    "teicoplanin": ("Teicoplanin", "Glycopeptide", "Cell wall (D-Ala-D-Ala)", "cidal"),
    "daptomycin": ("Daptomycin", "Lipopeptide", "Membrane", "cidal"),
    "oxacillin": ("Oxacillin", "Penicillin", "Cell wall (PBP)", "cidal"),
    "ampicillin": ("Ampicillin", "Penicillin", "Cell wall (PBP)", "cidal"),
    "amoxicillin": ("Amoxicillin", "Penicillin", "Cell wall (PBP)", "cidal"),
    "piperacillin": ("Piperacillin", "Penicillin", "Cell wall (PBP)", "cidal"),
    "nitrofurantoin": ("Nitrofurantoin", "Nitrofuran", "Multi-target", "cidal"),
}
ABX_SKIP = {"polymyxin"}

FAMILY = {
    "myoviridae": "Myoviridae (legacy)", "siphoviridae": "Siphoviridae (legacy)",
    "podoviridae": "Podoviridae (legacy)", "ackermannviridae": "Ackermannviridae",
    "straboviridae": "Straboviridae", "drexlerviridae": "Drexlerviridae",
    "autographiviridae": "Autographiviridae", "demerecviridae": "Demerecviridae",
    "casjensviridae": "Casjensviridae", "chaseviridae": "Chaseviridae",
    "mesyanzhinovviridae": "Mesyanzhinovviridae", "zobellviridae": "Zobellviridae",
    "kayfunaviridae": "Kayfunaviridae", "schitoviridae": "Schitoviridae",
    "webervirus": "Webervirus", "jedunavirus": "Jedunavirus", "przondovirus": "Przondovirus",
    "slopekvirus": "Slopekvirus", "tequatrovirus": "Tequatrovirus",
    "sugarlandvirus": "Sugarlandvirus", "drulisvirus": "Drulisvirus",
    "kp34virus": "Kp34virus", "marthavirus": "Marthavirus", "menlowvirus": "Menlowvirus",
    "miltonvirus": "Miltonvirus", "taipeivirus": "Taipeivirus", "gaprivervirus": "Gaprivervirus",
    "dietzellavirus": "Dietzellavirus", "fionnavirus": "Fionnavirus", "kuttervirus": "Kuttervirus",
    "viunavirus": "Viunavirus", "seongnamvirus": "Seongnamvirus", "kuravirus": "Kuravirus",
    "nonagvirus": "Nonagvirus", "jiaoyazivirus": "Jiaoyazivirus", "biseptimavirus": "Biseptimavirus",
    "melkshamvirus": "Melkshamvirus", "phikmvvirus": "Phikmvvirus", "lit1virus": "Lit1virus",
    "pikminvirus": "Pikminvirus", "chewzvirus": "Chewzvirus", "baoshanvirus": "Baoshanvirus",
}
MORPHO = [(r"\bjumbo ?phage", "Jumbo phage"), (r"\bsiphovir\w*", "Siphoviridae (morphotype)"),
          (r"\bmyovir\w*", "Myoviridae (morphotype)"), (r"\bpodovir\w*", "Podoviridae (morphotype)")]

REVIEW_TITLE = re.compile(
    r"\b(review|meta-analysis|systematic review|perspective|overview|advances in|"
    r"current landscape|challenges and opportunities|editorial|commentary|opinion|"
    r"state of the art|recent advances|an update|future directions|landscape|"
    r"bibliometric|treatment strategies|clinical challenges)\b", re.I)
REVIEW_BODY = re.compile(r"\b(this review|in this review|we summaris\w*|we summariz\w*|"
                         r"in the present review|literature review)\b", re.I)
PRIMARY = re.compile(r"\b(?:we|our (?:study|group|results|work))\s+(?:\w+\s+){0,2}"
                     r"(?:isolated|tested|performed|conducted|evaluated|assessed|measured|found|"
                     r"observed|constructed|cultured|screened|report|describe|characteri[sz]ed|developed|established)", re.I)

BACKGROUND = re.compile(
    r"\(\s*[A-Z][A-Za-z.\-]+(?:\s+et al\.)?,?\s*\d{4}\s*\)|et al\.|for instance|"
    r"studies (?:suggest|have shown|show)|has been (?:reported|shown|demonstrated)|"
    r"is (?:known|thought) to|mechanism of action|conventionally|traditionally|"
    r"such as |e\.g\.|reviewed in", re.I)
EXP_CONTEXT = re.compile(
    r"synerg|combination|combined|\bPAS\b|co-?administ|together with|with phage|"
    r"phage[- ]antibiotic|antagonis|additive|FIC|checkerboard|time[\s\-]kill|"
    r"biofilm|fold|log|reduction|inhibit|kill|eradica", re.I)


def sentences(t):
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", re.sub(r"\s+", " ", t)) if s.strip()]


def uniq(seq):
    out, seen = [], set()
    for s in seq:
        if s.lower() not in seen:
            seen.add(s.lower())
            out.append(s)
    return out


def organism_of(title, text):
    for o, p in ORG_PAT.items():
        if re.search(p, title, re.I):
            return o
    found = [o for o, p in ORG_PAT.items() if re.search(p, text[:6000], re.I)]
    return found[0] if found else "Unassigned"


def family_of(text):
    tl = text.lower()
    hits = uniq([v for k, v in FAMILY.items() if re.search(r"\b" + k, tl)])
    if not hits:
        hits = uniq([v for p, v in MORPHO if re.search(p, tl)])
    return "; ".join(hits[:2]) or "Not reported"


PHAGE_TOKEN = re.compile(r"(vB[_\-][A-Za-z0-9]+(?:[_\-][A-Za-z0-9\-]+)?|\bphi[A-Z0-9]\w{0,10}"
                         r"|φ[A-Z0-9]\w{0,10}|\bpKp\d\w{0,6}|\bKp[A-Za-z0-9]{1,3}\d{1,3}\b"
                         r"|\bSH\s?Kp\w+|\bKPP?\d{1,3}\b|\bP\d{2,3}[A-Za-z]?\b"
                         r"|\b(vB|phi|PA|Ab|S|EF|Ent)\w{0,2}\d{1,3}\b)")
PHAGE_STOP = re.compile(r"^(kpc|ndm|oxa|mic|fic|pas|mrsa|vre|crab|crkp|mdr|atcc|pbs|pcr|sd\d?|ci\d?)$", re.I)


def phage_name(text):
    c = Counter()
    for m in PHAGE_TOKEN.finditer(text):
        t = m.group(1).strip().rstrip(".,;:)")
        if len(t) < 3 or PHAGE_STOP.match(t):
            continue
        c[t] += 1
    if not c:
        return ""
    def rank(it):
        t, n = it
        b = 5 if t.lower().startswith("vb_") else (3 if t.lower().startswith(("phi", "φ")) else 0)
        return (b + min(n, 5), n)
    return "; ".join(t for t, _ in sorted(c.items(), key=rank, reverse=True)[:3])


GENO = [
    (r"kpc[- ]?2", "blaKPC-2"), (r"kpc[- ]?3", "blaKPC-3"), (r"\bkpc\b", "blaKPC"),
    (r"ndm[- ]?1", "blaNDM-1"), (r"ndm[- ]?5", "blaNDM-5"), (r"\bndm\b", "blaNDM"),
    (r"oxa[- ]?48", "blaOXA-48"), (r"oxa[- ]?23", "blaOXA-23"), (r"oxa[- ]?24", "blaOXA-24"),
    (r"oxa[- ]?58", "blaOXA-58"), (r"\bvim\b", "blaVIM"), (r"imp[- ]?\d+", "blaIMP"),
    (r"\bmeca\b|\bmecA\b", "mecA (MRSA)"), (r"\bvanA\b", "vanA (VRE)"),
    (r"\bvanB\b", "vanB (VRE)"), (r"ctx-m[- ]?\d+", "blaCTX-M"), (r"\bampc\b", "AmpC"),
    (r"\bges[- ]?\d+", "blaGES"), (r"\bsme\b", "blaSME"),
]


def genotype(text):
    tl = text.lower()
    found = uniq([n for p, n in GENO if re.search(p, tl)])
    return "/".join(found[:3])


ASSAY = [("Checkerboard", r"checkerboard|fractional inhibitory|\bFIC\b|FIC\s*index|FICi"),
         ("Time-Kill", r"time[\s\-]?kill|killing curve|kill(?:ing)? curve|growth curve|lysis curve"),
         ("Biofilm", r"biofilm|antibiofilm|sessile|\bMBEC\b|crystal violet|biomass"),
         ("In Vivo", r"in vivo|murine|mice\b|mouse\b|Galleria|mellonella|rat\b|rabbit|zebrafish|"
                     r"larvae|neutropenic|survival rate|wound")]


def assays(text):
    return [n for n, p in ASSAY if re.search(p, text, re.I)]


def sequence(text):
    tl = text.lower()
    if re.search(r"pre-?phage\s+(?:treatment|exposure|therapy)\s*(?:followed by|then|prior to)|"
                 r"phage\s+pre-?treatment\s*(?:followed by|then)|"
                 r"phage (?:was |were )?(?:added|administered|given)\s*(?:\d+\s*(?:h|hours?))?\s*(?:before|prior to)|"
                 r"phage[- ]first", tl):
        return "Phage-first"
    if re.search(r"antibiotic\s*(?:was |were )?(?:added|administered|given)\s*(?:\d+\s*(?:h|hours?))?\s*(?:before|prior to)|"
                 r"antibiotic[- ]first|antibiotic\s+pre-?treatment", tl):
        return "Antibiotic-first"
    if re.search(r"simultane\w*|concurrent\w*|at the same time|co-?administ\w*|co-?incubat\w*|"
                 r"in combination with|together with|added together", tl):
        return "Simultaneous"
    return "Not reported"


CONC = [r"(?:\d+\s*/\s*\d+|\d*\.\d+|\d+)\s*[×x]\s*(?:the\s*)?MIC",
        r"sub[\s\-]?(?:MIC|inhibitory|lethal)\w*(?:\s*concentrations?)?",
        r"\d+(?:\.\d+)?\s*(?:mg/L|mg/kg|μg/mL|µg/mL|ug/mL|μg/ml)",
        r"MIC(?:s)?\b"]


def concentration(text):
    for p in CONC:
        m = [x.group(0).strip() for x in re.finditer(p, text, re.I)]
        m = [x for x in m if not re.fullmatch(r"MICs?", x, re.I)]
        if m:
            return "; ".join(uniq(m)[:3])
    return ""


NUM = r"\d+(?:\.\d+)?"
OUTCOME_PATS = [
    (rf"FIC(?:\s*index|\s*i)?\s*(?:of|was|=|:|values? of)?\s*({NUM})", "FIC"),
    (rf"({NUM})\s*[-–]?\s*(?:log\s*10|log10|log)[\s\-]*(?:fold\s*)?(?:reduction|decrease|killing|inactivation|CFU)", "log10"),
    (rf"(?:reduc\w+|decreas\w+|lower\w*)\s*(?:of|by|in)?\s*(?:>|≥|~)?\s*({NUM})\s*[-–]?\s*(?:log\s*10|log10|log)\b", "log10"),
    (rf"({NUM})\s*%\s*(?:reduction|decrease|inhibition|inhibit\w*|killing|elimination|removal|eradication)", "percent"),
    (rf"(?:reduc\w+|inhibit\w*|decreas\w*)\s*(?:of|by|in)?\s*(?:>|≥|~)?\s*({NUM})\s*%", "percent"),
]


def outcome(text):
    out = []
    for p, lab in OUTCOME_PATS:
        for m in re.finditer(p, text, re.I):
            out.append(re.sub(r"\s+", " ", m.group(0)).strip())
    if not out:
        m = re.search(r"(complete (?:eradication|clearance|inhibition|killing)|bacterial eradication|"
                      r"multi-?log\w*\s*(?:CFU\s*)?(?:reduction|decrease)|complete lysis)", text, re.I)
        if m:
            return "Qualitative: " + m.group(0)
        return ""
    return " | ".join(uniq(out)[:3])


def classify(text, fics):
    if fics:
        f = min(fics)
        return ("Synergy" if f <= 0.5 else "Additive" if f <= 1 else
                "Indifferent" if f <= 4 else "Antagonism"), f"FIC {f}"
    tl = text.lower()
    if re.search(r"antagonis", tl):
        return "Antagonism", "text: antagonism"
    neg = re.search(r"(?:no|not|without|did not|failed to|absence of)\s+(\w+\s+){0,3}synerg|"
                    r"non-?synerg|no synergy|indifferent", tl)
    if neg:
        return "Indifferent", f"text negation: {neg.group(0)[:40]}"
    if re.search(r"synerg|enhanc\w+ (?:killing|activity|efficacy|bactericidal)", tl):
        return "Synergy", "text: synergy"
    if re.search(r"\badditive\b", tl):
        return "Additive", "text: additive"
    return "Not classified", ""


def pas_context(sents):
    return [s for s in sents if EXP_CONTEXT.search(s)
            and not (BACKGROUND.search(s) and not re.search(r"we |our |was (?:observed|found)|FIC", s, re.I))]


def rank_sents(sents, key):
    scored = []
    for s in sents:
        sl = s.lower()
        if key not in sl:
            continue
        sc = 0
        if re.search(r"synerg", sl):
            sc += 4
        if re.search(r"\bPAS\b|phage[- ]antibiotic synergy", s, re.I):
            sc += 3
        if re.search(r"combination|combined|co-?administ|together", sl):
            sc += 2
        if re.search(r"\d", s) and re.search(r"log|fold|%|FIC|mg/L|MIC", s, re.I):
            sc += 3
        if re.search(r"reduc|decreas|inhibit|kill|eradica|surviv", sl):
            sc += 2
        if BACKGROUND.search(s):
            sc -= 5
        if len(s) < 40:
            sc -= 3
        scored.append((sc, s))
    scored.sort(key=lambda x: -x[0])
    return [s for sc, s in scored if sc > 0]


# ------------------------------- main loop ---------------------------------
rows, excluded = [], []
for rec in corpus:
    text = rec["abstract"] + " " + " ".join(rec["sections"].values())
    sents = sentences(text)
    title = rec["title"]

    n_rev = len(REVIEW_BODY.findall(text[:6000]))
    n_prim = len(PRIMARY.findall(text[:6000]))
    if REVIEW_TITLE.search(title):
        excluded.append([rec["key"], "review-type title"])
        continue
    if n_rev >= 1 and n_prim == 0:
        excluded.append([rec["key"], "review prose, no primary experiment"])
        continue
    if not re.search(r"phage|bacteriophage", text, re.I):
        excluded.append([rec["key"], "no phage content"])
        continue

    org = organism_of(title, text)
    if org == "Unassigned":
        excluded.append([rec["key"], "no assignable ESKAPE organism"])
        continue

    SYN_TOK = re.compile(r"synerg|\bPAS\b|antagonis|additiv|phage[- ]antibiotic", re.I)
    COMBO_TOK = re.compile(r"combin|co-?administ|together with|plus|with phage|in the presence of", re.I)
    OTHER_SPECIES = re.compile(r"\b(?:S\.\s*aureus|MRSA|Staphylococcus|E\.\s*coli|Escherichia|"
                               r"Acinetobacter|Pseudomonas|Enterococcus|Salmonella|Candida)\b", re.I)

    pas_sents = pas_context(sents)
    arm_hits = {}
    for sent in pas_sents:
        has_syn = bool(SYN_TOK.search(sent))
        has_combo = bool(COMBO_TOK.search(sent))
        if not (has_syn or has_combo):
            continue
        # allow the study's own organism; reject cross-species only
        if OTHER_SPECIES.search(sent) and not re.search(
                r"klebsiella|staphylococc|aureus|acinetobacter|baumannii|pseudomonas|aeruginosa|"
                r"enterococc|enterobacter", sent, re.I):
            continue
        if re.search(r"minimum inhibitory concentrations?\s*\(MICs?\)|MICs? (?:was|were) determined|"
                     r"were determined (?:using|by) (?:broth )?microdilution|selection marker|"
                     r"conferred resistance|resistance marker", sent, re.I):
            continue
        for key, (canon, cls, mech, act) in ABX.items():
            if key in ABX_SKIP:
                continue
            if re.search(r"(?<![A-Za-z])" + re.escape(key) + r"(?![A-Za-z])", sent, re.I):
                e = arm_hits.setdefault(canon, {"key": key, "score": 0, "sents": [],
                                                "cls": cls, "mech": mech, "act": act})
                e["score"] += 3 if has_syn else 1
                e["sents"].append(sent)

    names = list(arm_hits)
    for g, prod in [("Ceftazidime", "Ceftazidime-avibactam"), ("Ampicillin", "Amoxicillin"),
                    ("Piperacillin", "Piperacillin-tazobactam")]:
        if g in arm_hits and prod in arm_hits:
            del arm_hits[g]

    found = sorted(arm_hits.items(), key=lambda kv: -kv[1]["score"])
    if not found:
        excluded.append([rec["key"], "no antibiotic in a phage-combination context"])
        continue

    study_assays = assays(text)
    study_outcome = outcome(" ".join(pas_sents) if pas_sents else text)
    seq_all = sequence(text)

    base = {
        "PMID_or_DOI": rec["pmid"] or rec["doi"],
        "DOI": rec["doi"], "PMCID": rec.get("pmcid", ""),
        "Publication_Year": rec["year"], "Organism": org,
        "Study_Title": title, "Journal": rec["journal"],
        "Open_Access": "Yes" if rec.get("is_open_access") else "No",
        "Fulltext_Available": "Yes" if rec["section_chars"] > 500 else "No",
        "Phage_Name": phage_name(text),
        "Phage_Family_or_Genus": family_of(text),
        "Resistance_Genotype": genotype(text),
        "Assay_Types_All": "; ".join(study_assays),
        "Dosing_Sequence_All": seq_all,
    }

    for canon, v in found:
        key = v["key"]
        ranked = rank_sents(v["sents"], key) or rank_sents(sents, key)
        ev = ranked[:3]
        ctx = " ".join(ev) if ev else (" ".join(pas_sents) or text[:2000])
        fics = []
        for m in re.finditer(r"FIC(?:\s*index|\s*i)?\s*(?:of|was|=|:)?\s*(\d*\.?\d+)", ctx, re.I):
            try:
                fics.append(float(m.group(1)))
            except ValueError:
                pass
        out = outcome(ctx) or study_outcome
        inter, why = classify(ctx, fics)
        seq = sequence(ctx) if sequence(ctx) != "Not reported" else seq_all
        asy = assays(ctx) or study_assays

        rows.append({
            **base,
            "Antibiotic_Tested": canon,
            "Antibiotic_Class": v["cls"],
            "Antibiotic_Mechanism": v["mech"],
            "Bactericidal_Static": v["act"],
            "Concentration_or_MIC_Fraction": concentration(ctx) or concentration(text),
            "Assay_Type": "; ".join(asy) if asy else "Not reported",
            "Dosing_Sequence": seq,
            "Quantitative_Outcome": out,
            "Interaction_Classification": inter,
            "Classification_Basis": why,
            "Evidence_Snippet": (ev[0][:600] if ev else ""),
        })

FIELDS = ["PMID_or_DOI", "DOI", "PMCID", "Publication_Year", "Organism", "Study_Title",
          "Journal", "Open_Access", "Fulltext_Available", "Phage_Name",
          "Phage_Family_or_Genus", "Resistance_Genotype", "Antibiotic_Tested",
          "Antibiotic_Class", "Antibiotic_Mechanism", "Bactericidal_Static",
          "Concentration_or_MIC_Fraction", "Assay_Type", "Dosing_Sequence",
          "Quantitative_Outcome", "Interaction_Classification", "Classification_Basis",
          "Assay_Types_All", "Dosing_Sequence_All", "Evidence_Snippet"]

with open(OUT / "eskape_pas_dataset.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in FIELDS})

(DATA / "eskape_excluded.json").write_text(json.dumps(excluded, indent=0))
print(f"rows={len(rows)} studies={len({r['PMID_or_DOI'] for r in rows})} excluded={len(excluded)}")
print("interaction:", dict(Counter(r["Interaction_Classification"] for r in rows)))
print("organisms:", dict(Counter(r["Organism"] for r in rows)))
