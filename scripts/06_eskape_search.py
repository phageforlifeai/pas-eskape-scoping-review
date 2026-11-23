#!/usr/bin/env python3
"""
Step 6 - ESKAPE-wide PAS search (Europe PMC + PubMed).

Searches phage-antibiotic synergy literature for all six ESKAPE pathogens plus
the Klebsiella data already collected. Writes data/eskape_pool.json.
"""
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests
# (PubMed via Europe PMC; Biopython not required for this sweep)

DATA = Path("/home/user/data")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

# ESKAPE organism blocks
ORG = {
    "Enterococcus faecium": '("Enterococcus faecium" OR "E. faecium" OR VRE OR "vancomycin-resistant Enterococcus")',
    "Staphylococcus aureus": '("Staphylococcus aureus" OR "S. aureus" OR MRSA OR "methicillin-resistant Staphylococcus")',
    "Klebsiella pneumoniae": '("Klebsiella pneumoniae" OR "K. pneumoniae" OR CRKP OR "carbapenem-resistant Klebsiella")',
    "Acinetobacter baumannii": '("Acinetobacter baumannii" OR "A. baumannii" OR CRAB OR "carbapenem-resistant Acinetobacter")',
    "Pseudomonas aeruginosa": '("Pseudomonas aeruginosa" OR "P. aeruginosa" OR MDRPA OR "carbapenem-resistant Pseudomonas")',
    "Enterobacter spp.": '("Enterobacter cloacae" OR "Enterobacter aerogenes" OR "E. cloacae" OR "Enterobacter")',
}
PHAGE = '(bacteriophage OR bacteriophages OR phage OR phages OR "phage therapy")'
COMBO = '("phage-antibiotic synergy" OR "phage antibiotic synergy" OR synerg* OR combination OR combined OR adjuvant OR "phage-antibiotic")'
YEAR = '(PUB_YEAR:[2015 TO 2026])'   # widen the window for the ESKAPE sweep


def epmc(q, page=200, cap=1200):
    hits, cur = [], "*"
    while True:
        r = requests.get(f"{EPMC}/search", params={
            "query": q, "format": "json", "resultType": "core",
            "pageSize": page, "cursorMark": cur}, timeout=90)
        if r.status_code != 200:
            time.sleep(5)
            continue
        js = r.json()
        b = js.get("resultList", {}).get("result", [])
        hits.extend(b)
        nxt = js.get("nextCursorMark")
        if not b or not nxt or nxt == cur or len(hits) >= cap:
            break
        cur = nxt
        time.sleep(0.3)
    return hits


pool = {}
for org, block in ORG.items():
    q = f"{PHAGE} AND {COMBO} AND {block} AND {YEAR}"
    res = epmc(q)
    for e in res:
        pmid = e.get("pmid") or ""
        doi = (e.get("doi") or "").lower()
        key = pmid or doi
        if not key:
            continue
        yr = str(e.get("pubYear") or "")
        if not (yr.isdigit() and 2015 <= int(yr) <= 2026):
            continue
        if key in pool:
            if org not in pool[key]["organisms"]:
                pool[key]["organisms"].append(org)
            continue
        pool[key] = {
            "pmid": pmid, "doi": e.get("doi") or "", "title": e.get("title") or "",
            "year": yr,
            "journal": (e.get("journalInfo", {}).get("journal", {}) or {}).get("title", ""),
            "abstract": e.get("abstractText") or "",
            "is_open_access": e.get("isOpenAccess") == "Y",
            "pmcid": e.get("pmcid") or "",
            "cited_by": e.get("citedByCount", 0),
            "organisms": [org],
        }
    print(f"{org:26s} hits={len(res):5d}  pool={len(pool)}")
    time.sleep(0.5)

out = list(pool.values())
(DATA / "eskape_pool.json").write_text(json.dumps(out, indent=1))
print(f"\nESKAPE pool: {len(out)} unique records")
print(f"  with PMCID (full text possible): {sum(1 for r in out if r['pmcid'])}")
print(f"  open access: {sum(1 for r in out if r['is_open_access'])}")
