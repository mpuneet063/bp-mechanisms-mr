#!/usr/bin/env python3
"""
Fetch every summary statistic the pipeline needs from the OpenGWAS API.

Run this on your own machine (it needs internet access to api.opengwas.io):

    pip install requests pandas openpyxl
    export OPENGWAS_JWT="eyJ..."          # token from https://api.opengwas.io/profile/
    python scripts/00_fetch_opengwas.py --keaton path/to/keaton2024_supp.xlsx   # --keaton is optional

It writes small files to ./data/raw/ in the folder you run it from (a few MB). Zip that folder and send it back.

What it does
1. Resolves an OpenGWAS study ID for each trait in STUDY_SPECS (pinned ID if it
   matches, otherwise a keyword search that prefers the named author/year,
   European ancestry and the largest sample). Every choice is written to
   data/raw/manifest.csv so it can be checked.
2. Builds the variant list: SBP genome-wide hits (Evangelou 2018, ieu-b-38),
   plus the Keaton 2024 rsIDs if the supplementary table is supplied, plus
   instruments for BMI, LDL and T2D (needed for multivariable MR).
3. Looks every variant up in every study (LD proxies r2 >= 0.8 allowed).
4. LD-clumps the BP variants at r2 < 0.001 and r2 < 0.05 (1000G EUR, server side).

It caches each API call under data/raw/cache/, so if it stops (rate limit,
network), just run it again and it carries on where it left off.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://api.opengwas.io/api/"
OUT = Path.cwd() / "data" / "raw"  # created in the folder you run the script from
CACHE = OUT / "cache"

# --------------------------------------------------------------------------- #
# Trait specification.
#   role:     exposure | cluster_trait | outcome | replication | mvmr
#   pin:      OpenGWAS id to use if its trait name matches `all`
#   all:      every keyword must appear in the trait name (case-insensitive)
#   any:      at least one must appear (optional)
#   exclude:  none may appear
#   author:   preferred first-author surname (optional)
#   year:     preferred year (optional)
#   prefix:   restrict to ids starting with this (optional)
#   binary:   rank candidates by number of cases instead of sample size
# Stroke/CAD/MI/HF/AF are deliberately NOT cluster traits (avoids circularity).
# --------------------------------------------------------------------------- #
STUDY_SPECS = [
    # BP itself (exposure + part of the clustering profile)
    dict(key="sbp", role="exposure", pin="ieu-b-38", all=["systolic blood pressure"], author="Evangelou"),
    dict(key="dbp", role="cluster_trait", pin="ieu-b-39", all=["diastolic blood pressure"], author="Evangelou"),
    dict(key="pp", role="cluster_trait", all=["pulse pressure"], author="Evangelou"),
    # adiposity / anthropometry
    dict(key="bmi", role="cluster_trait", pin="ieu-b-40", all=["body mass index"],
         exclude=["adjusted", "childhood", "change", "age"], author="Yengo"),
    dict(key="whradjbmi", role="cluster_trait", all=["waist", "hip"], any=["bmi", "body mass"], author="Pulit"),
    dict(key="bodyfat", role="cluster_trait", all=["body fat percentage"], exclude=["arm", "leg", "trunk"]),
    dict(key="height", role="cluster_trait", all=["height"],
         exclude=["sitting", "standing", "child", "birth", "age 10", "comparative"], author="Yengo"),
    # lipids
    dict(key="ldl", role="cluster_trait", pin="ieu-b-110", all=["ldl cholesterol"]),
    dict(key="hdl", role="cluster_trait", pin="ieu-b-109", all=["hdl cholesterol"]),
    dict(key="tg", role="cluster_trait", pin="ieu-b-111", all=["triglycerides"]),
    # glycaemic
    dict(key="t2d", role="cluster_trait", all=["type 2 diabetes"],
         exclude=["adjusted", "bmi", "family", "age", "neph", "retin", "neuro"], author="Mahajan", binary=True),
    dict(key="fg", role="cluster_trait", all=["fasting glucose"], exclude=["interaction", "adjusted for bmi"]),
    dict(key="fi", role="cluster_trait", all=["fasting insulin"], exclude=["interaction", "adjusted for bmi"]),
    dict(key="hba1c", role="cluster_trait", any=["hba1c", "glycated haemoglobin", "glycated hemoglobin"]),
    # kidney / electrolytes
    dict(key="egfr", role="cluster_trait", any=["glomerular filtration", "egfr"],
         exclude=["cystatin", "change", "decline", "diabet"], author="Wuttke"),
    dict(key="urate", role="cluster_trait", any=["urate", "uric acid"], exclude=["gout", "excretion"]),
    dict(key="uacr", role="cluster_trait", any=["albumin-to-creatinine", "albumin to creatinine",
                                                "albumin/creatinine", "microalbumin in urine"]),
    dict(key="na_urine", role="cluster_trait", all=["sodium in urine"]),
    dict(key="k_urine", role="cluster_trait", all=["potassium in urine"]),
    # cardiac / inflammation / liver / behaviour
    dict(key="hr", role="cluster_trait", all=["heart rate"],
         exclude=["variability", "recovery", "increase", "during", "maximum"]),
    dict(key="crp", role="cluster_trait", all=["c-reactive protein"]),
    dict(key="alt", role="cluster_trait", all=["alanine aminotransferase"]),
    dict(key="alcohol", role="cluster_trait", any=["drinks per week", "alcohol consumption", "alcoholic drinks per week"]),
    dict(key="smoking", role="cluster_trait", all=["smoking initiation"], author="Liu"),
    # outcomes
    dict(key="stroke_any", role="outcome", all=["stroke"],
         exclude=["ischemic", "ischaemic", "large", "small", "cardio", "lacunar", "hemorrh",
                  "haemorrh", "family", "father", "mother", "sibling", "age", "outcome", "recurr"],
         author="Mishra", year=2022, binary=True),
    dict(key="stroke_is", role="outcome", any=["ischemic stroke", "ischaemic stroke"],
         exclude=["large", "small", "cardio", "lacunar", "outcome", "recurr", "excluding"],
         author="Mishra", year=2022, binary=True),
    dict(key="stroke_las", role="outcome", all=["large artery"], author="Mishra", year=2022, binary=True),
    dict(key="stroke_ces", role="outcome", all=["cardioembolic"], author="Mishra", year=2022, binary=True),
    dict(key="stroke_svs", role="outcome", all=["small vessel", "stroke"], author="Mishra", year=2022, binary=True),
    dict(key="ich", role="outcome", all=["intracerebral"], any=["hemorrhage", "haemorrhage"], binary=True),
    dict(key="cad", role="outcome", all=["coronary artery disease"],
         exclude=["family", "age", "father", "mother", "sibling", "self-reported"],
         author="Aragam", year=2022, binary=True),
    # replication outcomes (FinnGen release in OpenGWAS)
    dict(key="fg_stroke", role="replication", prefix="finn-b-", all=["stroke"], any=["excluding sah", "excl sah", "excluding subarachnoid"], binary=True),
    dict(key="fg_is", role="replication", prefix="finn-b-", any=["ischemic stroke", "ischaemic stroke"],
         exclude=["large", "small", "cardio", "lacunar", "ichd"], binary=True),
    dict(key="fg_chd", role="replication", prefix="finn-b-", any=["coronary heart disease", "major coronary"], binary=True),
]

MVMR_KEYS = ["bmi", "ldl", "t2d"]  # extra instruments pulled for multivariable MR


# --------------------------------------------------------------------------- #
# HTTP helpers
# --------------------------------------------------------------------------- #
class Api:
    def __init__(self, jwt: str):
        self.s = requests.Session()
        self.s.headers.update({"Authorization": f"Bearer {jwt}", "X-API-SOURCE": "bp-mechanisms-mr/0.1"})

    def call(self, path: str, payload: dict | None = None, cache_key: str | None = None, tries: int = 8):
        if cache_key:
            f = CACHE / f"{cache_key}.json"
            if f.exists():
                return json.loads(f.read_text())
        url = API + path
        for attempt in range(tries):
            try:
                if payload is None:
                    r = self.s.get(url, timeout=300)
                else:
                    r = self.s.post(url, json=payload, timeout=300)
                    if r.status_code in (400, 415, 422):  # older endpoints accept form data
                        r = self.s.post(url, data=payload, timeout=300)
            except requests.RequestException as e:
                wait = 10 * (attempt + 1)
                log(f"  network error ({e}); retry in {wait}s")
                time.sleep(wait)
                continue
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 60 * (attempt + 1)))
                log(f"  rate limited; waiting {wait}s (OpenGWAS allowance resets periodically)")
                time.sleep(min(wait, 900))
                continue
            if r.status_code >= 500:
                time.sleep(15 * (attempt + 1))
                continue
            if r.status_code == 401:
                sys.exit("401 Unauthorized: check OPENGWAS_JWT (get one at https://api.opengwas.io/profile/)")
            r.raise_for_status()
            out = r.json()
            if cache_key:
                (CACHE / f"{cache_key}.json").write_text(json.dumps(out))
            return out
        raise RuntimeError(f"Giving up on {path} after {tries} attempts")


LOGF = None


def log(msg: str):
    print(msg, flush=True)
    if LOGF:
        LOGF.write(msg + "\n")
        LOGF.flush()


# --------------------------------------------------------------------------- #
# Study resolution
# --------------------------------------------------------------------------- #
def _matches(spec, trait: str) -> bool:
    t = (trait or "").lower()
    if any(k.lower() not in t for k in spec.get("all", [])):
        return False
    if spec.get("any") and not any(k.lower() in t for k in spec["any"]):
        return False
    if any(k.lower() in t for k in spec.get("exclude", [])):
        return False
    return True


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def resolve(spec, catalog: pd.DataFrame):
    pin = spec.get("pin")
    if pin and pin in catalog.index and _matches(spec, catalog.loc[pin, "trait"]):
        return pin, "pinned", catalog.loc[[pin]]
    c = catalog[catalog["trait"].apply(lambda t: _matches(spec, t))]
    if spec.get("prefix"):
        c = c[c.index.str.startswith(spec["prefix"])]
    elif spec["role"] == "outcome":  # main outcomes come from consortium GWAS, not biobank-only releases
        c = c[~c.index.str.startswith(("finn-b-", "ukb-"))]
    if c.empty:
        return None, "no match", c
    c = c.copy()
    pop = c["population"].fillna("").str.lower()
    c["_eur"] = pop.str.contains("european").astype(int)
    c["_mixed"] = pop.str.contains("mixed").astype(int)
    c = c[(c["_eur"] == 1) | (c["_mixed"] == 1)] if ((c["_eur"] + c["_mixed"]) > 0).any() else c
    c["_author"] = c["author"].fillna("").str.contains(spec.get("author", "@@none@@"), case=False).astype(int)
    c["_year"] = (c["year"].apply(_num) == spec.get("year", -1)).astype(int)
    c["_size"] = c["ncase"].apply(_num) if spec.get("binary") else c["sample_size"].apply(_num)
    c["_size"] = c["_size"].where(c["_size"] > 0, c["sample_size"].apply(_num))
    c = c.sort_values(["_author", "_year", "_eur", "_size"], ascending=False)
    how = "searched" + (" (author match)" if c["_author"].iloc[0] else "")
    if pin:
        how += f"; pin {pin} did not match"
    return c.index[0], how, c


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def distance_prune(df: pd.DataFrame, window_bp: int) -> list[str]:
    """Keep the most significant variant, drop others within window_bp on the same chromosome, repeat."""
    d = df.dropna(subset=["chr", "position"]).sort_values("p")
    kept, taken = [], {}
    for r in d.itertuples():
        pos = taken.setdefault(str(r.chr), [])
        if all(abs(int(r.position) - q) > window_bp for q in pos):
            kept.append(r.rsid)
            pos.append(int(r.position))
    return kept


def keaton_rsids(path: Path) -> list[str]:
    """Pull every rsID from the sheet that contains the most of them."""
    sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=str)
    best = []
    for name, df in sheets.items():
        ids = set()
        for col in df.columns:
            vals = df[col].dropna().astype(str).str.strip()
            hits = vals[vals.str.match(r"^rs\d+$")]
            ids.update(hits)
        if len(best) < len(ids) <= 5000:  # the signal table has ~2,103 rows; ignore any huge weight tables
            best, best_name = sorted(ids), name
    log(f"Keaton file: {len(best)} rsIDs from sheet '{best_name}'" if best else "Keaton file: no rsIDs found")
    return best


def main():
    global LOGF
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--keaton", type=Path, help="Keaton et al. 2024 supplementary tables (.xlsx), optional")
    ap.add_argument("--chunk", type=int, default=100, help="variants per association call")
    ap.add_argument("--no-proxies", action="store_true")
    args = ap.parse_args()

    jwt = os.environ.get("OPENGWAS_JWT", "").strip()
    if not jwt:
        sys.exit("Set OPENGWAS_JWT first (free token: https://api.opengwas.io/profile/)")
    CACHE.mkdir(parents=True, exist_ok=True)
    LOGF = open(OUT / "fetch_log.txt", "a")
    log(f"=== fetch started {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    api = Api(jwt)
    log(f"API status: {api.call('status')}")

    # 1. study catalogue + resolution
    raw = api.call("gwasinfo", cache_key="gwasinfo_all")
    recs = list(raw.values()) if isinstance(raw, dict) else raw
    catalog = pd.DataFrame(recs).set_index("id")
    for col in ["trait", "author", "year", "population", "sample_size", "ncase", "ncontrol", "consortium", "pmid", "nsnp", "unit", "sd"]:
        if col not in catalog:
            catalog[col] = None
    log(f"{len(catalog)} studies in OpenGWAS")

    man_rows, alt_rows = [], []
    for spec in STUDY_SPECS:
        sid, how, cands = resolve(spec, catalog)
        if sid is None:
            log(f"  [{spec['key']}] NOT FOUND - skipped")
            man_rows.append(dict(key=spec["key"], role=spec["role"], id=None, how=how))
            continue
        r = catalog.loc[sid]
        man_rows.append(dict(key=spec["key"], role=spec["role"], id=sid, how=how, trait=r["trait"],
                             author=r["author"], year=r["year"], population=r["population"],
                             sample_size=r["sample_size"], ncase=r["ncase"], ncontrol=r["ncontrol"],
                             consortium=r["consortium"], pmid=r["pmid"], unit=r["unit"], sd=r["sd"]))
        for alt_id, a in cands.head(6).iterrows():
            alt_rows.append(dict(key=spec["key"], id=alt_id, trait=a["trait"], author=a["author"], year=a["year"],
                                 population=a["population"], sample_size=a["sample_size"], ncase=a["ncase"]))
        log(f"  [{spec['key']}] {sid}  {r['trait']} | {r['author']} {r['year']} | n={r['sample_size']} ({how})")
    manifest = pd.DataFrame(man_rows)
    manifest.to_csv(OUT / "manifest.csv", index=False)
    pd.DataFrame(alt_rows).to_csv(OUT / "manifest_alternatives.csv", index=False)
    ids = {r.key: r.id for r in manifest.itertuples() if isinstance(r.id, str)}
    if "sbp" not in ids:
        sys.exit("Could not resolve the SBP study; stopping.")

    # 2. variant list
    inst_rows = []
    th = api.call("tophits", dict(id=[ids["sbp"]], pval=5e-8, preclumped=1, clump=1, r2=0.001, kb=10000, pop="EUR"),
                  cache_key=f"tophits_{ids['sbp']}_strict")
    for h in th:
        inst_rows.append(dict(rsid=h["rsid"], source="sbp_tophits_r2_0.001"))
    try:  # server-side clumping of all SBP hits at a relaxed threshold; can be slow or refused
        th2 = api.call("tophits", dict(id=[ids["sbp"]], pval=5e-8, preclumped=0, clump=1, r2=0.05, kb=500, pop="EUR"),
                       cache_key=f"tophits_{ids['sbp']}_r2_0.05", tries=3)
    except Exception as e:  # noqa: BLE001
        log(f"  relaxed SBP clumping failed ({e}); continuing with strict hits + Keaton list")
        th2 = []
    for h in th2:
        inst_rows.append(dict(rsid=h["rsid"], source="sbp_tophits_r2_0.05"))
    log(f"SBP tophits: {len(th)} (r2<0.001), {len(th2)} (r2<0.05, 500kb)")
    if args.keaton:
        for rs in keaton_rsids(args.keaton):
            inst_rows.append(dict(rsid=rs, source="keaton2024"))
    for k in MVMR_KEYS:
        if k in ids:
            t = api.call("tophits", dict(id=[ids[k]], pval=5e-8, preclumped=1, clump=1, r2=0.001, kb=10000, pop="EUR"),
                         cache_key=f"tophits_{ids[k]}_strict")
            for h in t:
                inst_rows.append(dict(rsid=h["rsid"], source=f"mvmr_{k}"))
            log(f"MVMR instruments {k}: {len(t)}")
    inst = pd.DataFrame(inst_rows).drop_duplicates()
    snps = sorted(inst["rsid"].unique())
    log(f"{len(snps)} unique variants to look up")

    # 3. look-ups
    long_rows = []
    study_ids = sorted(set(ids.values()))
    for i, sid in enumerate(study_ids, 1):
        n_before = len(long_rows)
        for j, ch in enumerate(chunks(snps, args.chunk)):
            key = "assoc_" + sid + "_" + hashlib.md5(",".join(ch).encode()).hexdigest()[:12] + ("_np" if args.no_proxies else "")
            res = api.call("associations", dict(variant=ch, id=[sid], proxies=0 if args.no_proxies else 1, r2=0.8,
                                                 align_alleles=1, palindromes=1, maf_threshold=0.3), cache_key=key)
            if isinstance(res, list):
                long_rows.extend(res)
        log(f"  [{i}/{len(study_ids)}] {sid}: {len(long_rows) - n_before} associations")
    assoc = pd.DataFrame(long_rows)
    keep = [c for c in ["id", "trait", "rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p", "n",
                        "proxy", "target_snp", "proxy_snp", "rsq"] if c in assoc.columns]
    assoc = assoc[keep]
    assoc.to_csv(OUT / "associations.csv.gz", index=False, compression="gzip")
    log(f"associations.csv.gz: {len(assoc)} rows")

    # 4. LD clumping of BP variants (server side, 1000G EUR)
    bp = assoc[(assoc["id"] == ids["sbp"]) & assoc["rsid"].isin(inst.loc[~inst["source"].str.startswith("mvmr"), "rsid"])]
    bp = bp.sort_values("p").drop_duplicates("rsid")
    for r2, kb in [(0.001, 10000), (0.05, 500)]:
        try:
            kept = api.call("ld/clump", dict(rsid=bp["rsid"].tolist(), pval=bp["p"].astype(float).tolist(),
                                             pthresh=1.0, r2=r2, kb=kb, pop="EUR"),
                            cache_key=f"clump_{r2}_{kb}_{hashlib.md5(','.join(bp['rsid']).encode()).hexdigest()[:12]}",
                            tries=4)
            how = "LD clumping (1000G EUR)"
        except Exception as e:  # noqa: BLE001 - fall back to distance-based pruning (conservative)
            log(f"  ld/clump failed ({e}); falling back to distance pruning")
            kept = distance_prune(bp, kb * 1000 if r2 >= 0.01 else 10_000_000)
            how = "distance pruning fallback"
        inst[f"clump_r2_{r2}"] = inst["rsid"].isin(set(kept))
        log(f"clump r2<{r2}, {kb}kb: {len(kept)} kept of {len(bp)} [{how}]")
    inst.to_csv(OUT / "instruments.csv", index=False)
    log("Done. Zip data/raw (you can leave out data/raw/cache) and send it back.")


if __name__ == "__main__":
    main()
