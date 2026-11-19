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
network, a runtime that dies), just run it again and it carries on where it left off.
On Colab, keep the repo on Google Drive so the cache survives runtime restarts.
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
#   author:   preferred first-author surname, or a list in order of preference (optional)
#   year:     preferred year (optional)
#   prefix:   restrict to ids starting with this (optional)
#   binary:   rank candidates by number of cases instead of sample size
# Stroke/CAD/MI/HF/AF are deliberately NOT cluster traits (avoids circularity).
# --------------------------------------------------------------------------- #
STUDY_SPECS = [
    # BP itself (exposure + part of the clustering profile)
    dict(key="sbp", role="exposure", pin="ieu-b-38", all=["systolic blood pressure"], author="Evangelou"),
    dict(key="dbp", role="cluster_trait", pin="ieu-b-39", all=["diastolic blood pressure"], author="Evangelou"),
    dict(key="pp", role="cluster_trait", all=["pulse pressure"], author=["Evangelou", "Keaton", "Warren"]),
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
    dict(key="hr", role="cluster_trait", any=["heart rate", "pulse rate"],
         exclude=["variability", "recovery", "increase", "during", "maximum", "foetal", "fetal", "labour", "labor",
                  "delivery", "diagnos", "icd", "complicat", "pregnan", "change", "response", "exercise"]),
    dict(key="crp", role="cluster_trait", all=["c-reactive protein"]),
    dict(key="alt", role="cluster_trait", all=["alanine aminotransferase"]),
    dict(key="alcohol", role="cluster_trait", any=["drinks per week", "alcohol consumption", "alcoholic drinks per week"]),
    dict(key="smoking", role="cluster_trait", all=["smoking initiation"], author="Liu"),
    # outcomes
    dict(key="stroke_any", role="outcome", all=["stroke"],
         exclude=["ischemic", "ischaemic", "large", "small", "cardio", "lacunar", "hemorrh",
                  "haemorrh", "family", "father", "mother", "sibling", "age", "outcome", "recurr"],
         author=["Mishra", "Malik"], binary=True),
    dict(key="stroke_is", role="outcome", any=["ischemic stroke", "ischaemic stroke"],
         exclude=["large", "small", "cardio", "lacunar", "outcome", "recurr", "excluding"],
         author=["Mishra", "Malik"], binary=True),
    # GIGASTROKE (Mishra 2022) is preferred if OpenGWAS has it; otherwise MEGASTROKE (Malik 2018).
    # MEGASTROKE writes subtypes with hyphens ("small-vessel"), so match both spellings.
    dict(key="stroke_las", role="outcome", all=["stroke"], any=["large artery", "large-artery"],
         author=["Mishra", "Malik"], binary=True),
    dict(key="stroke_ces", role="outcome", all=["stroke"], any=["cardioembolic", "cardio-embolic"],
         author=["Mishra", "Malik"], binary=True),
    dict(key="stroke_svs", role="outcome", all=["stroke"], any=["small vessel", "small-vessel"],
         author=["Mishra", "Malik"], binary=True),
    dict(key="ich", role="outcome", all=["intracerebral"], any=["hemorrhage", "haemorrhage"], binary=True),
    dict(key="cad", role="outcome", all=["coronary artery disease"],
         exclude=["family", "age", "father", "mother", "sibling", "self-reported"],
         author=["Aragam", "van der Harst", "Nikpay"], binary=True),  # consortium meta-analyses, newest first
    # replication outcomes (FinnGen release in OpenGWAS)
    dict(key="fg_stroke", role="replication", prefix="finn-b-", all=["stroke"], any=["excluding sah", "excl sah", "excluding subarachnoid"], binary=True),
    dict(key="fg_is", role="replication", prefix="finn-b-", any=["ischemic stroke", "ischaemic stroke"],
         exclude=["large", "small", "cardio", "lacunar", "ichd"], binary=True),
    dict(key="fg_chd", role="replication", prefix="finn-b-", any=["coronary heart disease", "major coronary"], binary=True),
]

MVMR_KEYS = ["bmi", "ldl", "t2d"]  # extra instruments pulled for multivariable MR
MIN_NSNP = 1_000_000                 # a study must cover at least this many variants to be used
MAX_ASSOC_PER_REQUEST = 64           # OpenGWAS API v4 limit: N(id) * N(variant) <= 64 per request


# --------------------------------------------------------------------------- #
# HTTP helpers
# --------------------------------------------------------------------------- #
class ClientError(Exception):
    """4xx other than 401/429: the request itself was rejected (bad id, bad variant, bad parameter)."""

    def __init__(self, status: int, body: str):
        super().__init__(f"HTTP {status}: {body[:300]}")
        self.status, self.body = status, body


class Api:
    def __init__(self, jwt: str):
        self.s = requests.Session()
        self.s.headers.update({"Authorization": f"Bearer {jwt}", "X-API-SOURCE": "bp-mechanisms-mr/0.1"})

    def call(self, path: str, payload: dict | None = None, cache_key: str | None = None, tries: int = 8):
        if cache_key:
            f = CACHE / f"{cache_key}.json"
            if f.exists():
                try:
                    return json.loads(f.read_text())
                except (json.JSONDecodeError, OSError):  # half-written by a runtime that died: fetch again
                    f.unlink(missing_ok=True)
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
                exp = (f" Your token expired at {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(TOKEN_EXP))}."
                       if TOKEN_EXP and TOKEN_EXP < time.time() else "")
                log(f"401 Unauthorized.{exp} Tokens last 14 days: make a new one at https://api.opengwas.io/profile/, "
                    "update OPENGWAS_JWT and re-run. Progress so far is cached and will be reused.")
                sys.exit(1)
            if 400 <= r.status_code < 500:
                raise ClientError(r.status_code, r.text)
            r.raise_for_status()
            out = r.json()
            if cache_key:  # atomic: write to a temp file, then rename, so a crash never leaves a partial file
                tmp = CACHE / f".{cache_key}.json.tmp"
                tmp.write_text(json.dumps(out))
                os.replace(tmp, CACHE / f"{cache_key}.json")
            return out
        raise RuntimeError(f"Giving up on {path} after {tries} attempts")


def jwt_expiry(token: str):
    """Expiry time (UTC epoch seconds) from the token's payload, or None. No verification needed for this."""
    import base64
    try:
        part = token.split(".")[1]
        payload = json.loads(base64.urlsafe_b64decode(part + "=" * (-len(part) % 4)))
        return float(payload["exp"])
    except Exception:  # noqa: BLE001
        return None


TOKEN_EXP = None
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
    # genome-wide studies only: exome/targeted arrays (e.g. Surendran 2020 pulse pressure) miss most variants
    nsnp = c["nsnp"].apply(_num)
    gw = c[(nsnp == 0) | (nsnp >= MIN_NSNP)]
    if not gw.empty:
        c = gw.copy()
    authors = spec.get("author", [])
    authors = [authors] if isinstance(authors, str) else authors
    c["_author"] = 0
    for rank, a in enumerate(reversed(authors), 1):  # earlier in the list = higher score
        c.loc[c["author"].fillna("").str.contains(a, case=False, regex=False), "_author"] = rank
    c["_year"] = (c["year"].apply(_num) == spec.get("year", -1)).astype(int)
    c["_size"] = c["ncase"].apply(_num) if spec.get("binary") else c["sample_size"].apply(_num)
    c["_size"] = c["_size"].where(c["_size"] > 0, c["sample_size"].apply(_num))
    c = c.sort_values(["_author", "_year", "_eur", "_size"], ascending=False)
    how = "searched" + (f" (author match: {c['author'].iloc[0]})" if c["_author"].iloc[0] else "")
    if pin:
        how += f"; pin {pin} did not match"
    return c.index[0], how, c


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def _record_done(sid: str, variants: list[str]):
    with open(CACHE / f"done_{sid}.txt", "a") as fh:
        fh.write("\n".join(variants) + "\n")


def _cached_study(sid: str):
    """Rows and finished variants already on disk for a study, whatever chunk size produced them."""
    rows, done = [], set()
    for f in CACHE.glob(f"assoc_{sid}_*.json"):
        try:
            res = json.loads(f.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if isinstance(res, list):
            rows.extend(res)
            for r in res:
                is_proxy = str(r.get("proxy", "")).lower() in ("true", "1")
                done.add(r.get("target_snp") if is_proxy and r.get("target_snp") else r.get("rsid"))
    df = CACHE / f"done_{sid}.txt"
    if df.exists():
        done.update(x for x in df.read_text().split() if x)
    return rows, done


def lookup(api: "Api", sid: str, variants: list[str], proxies: int, failed: list, depth: int = 0) -> list:
    """
    Associations for one study. On a rejected request (HTTP 4xx): log the server's reason, retry without
    LD proxies, then split the chunk in half recursively so one bad variant cannot sink the rest.
    """
    key = "assoc_" + sid + "_" + hashlib.md5(",".join(variants).encode()).hexdigest()[:12] + ("" if proxies else "_np")
    payload = dict(variant=variants, id=[sid], proxies=proxies, r2=0.8, align_alleles=1, palindromes=1,
                   maf_threshold=0.3)
    try:
        res = api.call("associations", payload, cache_key=key)
        _record_done(sid, variants)
        return res if isinstance(res, list) else []
    except ClientError as e:
        if depth == 0:
            log(f"  {sid}: request rejected ({e}); retrying")
        if proxies:
            try:
                res = api.call("associations", dict(payload, proxies=0), cache_key=key + "_np")
                _record_done(sid, variants)
                if depth == 0:
                    log(f"  {sid}: OK without LD proxies for this chunk")
                return res if isinstance(res, list) else []
            except ClientError:
                pass
        if len(variants) == 1:
            failed.append(dict(study=sid, variant=variants[0], reason=e.body[:200]))
            return []
        mid = len(variants) // 2
        return (lookup(api, sid, variants[:mid], proxies, failed, depth + 1)
                + lookup(api, sid, variants[mid:], proxies, failed, depth + 1))


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
    ap.add_argument("--chunk", type=int, default=MAX_ASSOC_PER_REQUEST,
                    help=f"variants per association call (API maximum {MAX_ASSOC_PER_REQUEST})")
    ap.add_argument("--no-proxies", action="store_true")
    args = ap.parse_args()
    args.chunk = max(1, min(args.chunk, MAX_ASSOC_PER_REQUEST))

    global TOKEN_EXP
    jwt = os.environ.get("OPENGWAS_JWT", "").strip()
    if not jwt:
        sys.exit("Set OPENGWAS_JWT first (free token: https://api.opengwas.io/profile/)")
    TOKEN_EXP = jwt_expiry(jwt)
    if TOKEN_EXP is not None:
        left_h = (TOKEN_EXP - time.time()) / 3600
        if left_h <= 0:
            sys.exit(f"OPENGWAS_JWT expired at {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(TOKEN_EXP))}. "
                     "Make a new one at https://api.opengwas.io/profile/ and re-run.")
    CACHE.mkdir(parents=True, exist_ok=True)
    LOGF = open(OUT / "fetch_log.txt", "a")
    log(f"=== fetch started {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    if TOKEN_EXP is not None:
        log(f"token valid until {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(TOKEN_EXP))} "
            f"({(TOKEN_EXP - time.time()) / 86400:.1f} days left)")
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
    snps = sorted({str(x).strip() for x in inst["rsid"].dropna()})
    bad = [x for x in snps if not re.match(r"^(rs\d+|\d{1,2}:\d+)$", x)]
    if bad:
        log(f"  dropping {len(bad)} malformed variant ids, e.g. {bad[:5]}")
    snps = [x for x in snps if x not in set(bad)]
    log(f"{len(snps)} unique variants to look up")

    # 3. look-ups
    long_rows, failed, dead_studies = [], [], []
    study_ids = sorted(set(ids.values()))
    probe = th[0]["rsid"] if th else snps[0]  # an SBP lead variant: valid in any study
    for i, sid in enumerate(study_ids, 1):
        n_before, n_failed = len(long_rows), len(failed)
        try:  # is the study itself usable? one cheap request before 30+ chunks
            api.call("associations", dict(variant=[probe], id=[sid], proxies=0, r2=0.8, align_alleles=1,
                                          palindromes=1, maf_threshold=0.3), cache_key=f"probe_{sid}_{probe}")
        except ClientError as e:
            log(f"  [{i}/{len(study_ids)}] {sid}: SKIPPED - the API rejects this study ({e})")
            dead_studies.append(dict(study=sid, reason=e.body[:300]))
            continue
        cached_rows, done = _cached_study(sid)
        long_rows.extend(cached_rows)
        todo = [v for v in snps if v not in done]
        if done:
            log(f"  [{i}/{len(study_ids)}] {sid}: {len(snps) - len(todo)} variants already cached, {len(todo)} to fetch")
        n_chunks = (len(todo) + args.chunk - 1) // args.chunk
        for j, ch in enumerate(chunks(todo, args.chunk), 1):
            long_rows.extend(lookup(api, sid, ch, 0 if args.no_proxies else 1, failed))
            if j % 10 == 0 and j < n_chunks:
                print(f"    {sid}: chunk {j}/{n_chunks}", flush=True)
        extra = f", {len(failed) - n_failed} variants rejected" if len(failed) > n_failed else ""
        log(f"  [{i}/{len(study_ids)}] {sid}: {len(long_rows) - n_before} associations{extra}")
    if failed:
        pd.DataFrame(failed).to_csv(OUT / "rejected_variants.csv", index=False)
    if dead_studies:
        pd.DataFrame(dead_studies).to_csv(OUT / "rejected_studies.csv", index=False)
        dead = {d["study"] for d in dead_studies}
        manifest.loc[manifest["id"].isin(dead), "how"] += " | REJECTED by API at look-up"
        manifest.to_csv(OUT / "manifest.csv", index=False)
        if ids["sbp"] in dead:
            sys.exit(f"The SBP study {ids['sbp']} was rejected by the API; see rejected_studies.csv")
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
