"""Load the fetched OpenGWAS files and harmonise every association to the SBP-raising allele."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass
class RawData:
    manifest: pd.DataFrame      # one row per trait key: key, role, id, trait, unit, ncase, ncontrol, ...
    assoc: pd.DataFrame         # long table: id, rsid, ea, nea, eaf, beta, se, p, n
    instruments: pd.DataFrame   # rsid, source, clump flags

    def study(self, key: str) -> str | None:
        row = self.manifest.loc[self.manifest["key"] == key, "id"]
        return row.iloc[0] if len(row) and isinstance(row.iloc[0], str) else None

    def keys(self, role: str) -> list[str]:
        m = self.manifest
        return [k for k, i in zip(m.loc[m["role"] == role, "key"], m.loc[m["role"] == role, "id"]) if isinstance(i, str)]


def load_raw(raw_dir: str | Path) -> RawData:
    raw_dir = Path(raw_dir)
    manifest = pd.read_csv(raw_dir / "manifest.csv")
    assoc = pd.read_csv(raw_dir / "associations.csv.gz")
    inst = pd.read_csv(raw_dir / "instruments.csv")
    for c in ["beta", "se", "p", "eaf", "n"]:
        if c in assoc:
            assoc[c] = pd.to_numeric(assoc[c], errors="coerce")
    assoc["ea"] = assoc["ea"].astype(str).str.upper()
    assoc["nea"] = assoc["nea"].astype(str).str.upper()
    assoc = assoc.dropna(subset=["beta", "se"])
    assoc = assoc[assoc["se"] > 0]
    # a proxy hit is reported for the variant that was asked for (target_snp), not the proxy's own rsid
    if "proxy" in assoc and "target_snp" in assoc:
        is_proxy = assoc["proxy"].astype(str).str.lower().isin(["true", "1"]) & assoc["target_snp"].notna()
        assoc.loc[is_proxy, "rsid"] = assoc.loc[is_proxy, "target_snp"]
    # if a variant was found both directly and via a proxy, keep the direct hit
    if "proxy" in assoc:
        assoc["_proxy"] = assoc["proxy"].astype(str).str.lower().isin(["true", "1"])
        assoc = assoc.sort_values("_proxy").drop_duplicates(["id", "rsid"]).drop(columns="_proxy")
    else:
        assoc = assoc.drop_duplicates(["id", "rsid"])
    return RawData(manifest, assoc, inst)


def _is_palindromic(a1: pd.Series, a2: pd.Series) -> pd.Series:
    return a1.map(COMPLEMENT) == a2


def harmonise(raw: RawData, ref_key: str = "sbp", maf_ambiguous: float = 0.42) -> pd.DataFrame:
    """
    Return a long table with every association re-expressed for the SBP-raising allele.

    Rules (TwoSampleMR action = 2 equivalent):
      * reference allele = allele that raises SBP in the reference study
      * same alleles -> keep sign; swapped -> flip sign and eaf
      * strand-flipped (complement) alleles -> treated like the above after complementing
      * palindromic SNPs: kept only if eaf is informative (outside [0.42, 0.58]) and consistent
    Output columns: rsid, key, beta, se, p, n, eaf (all aligned), proxy flag.
    """
    sid = raw.study(ref_key)
    ref = raw.assoc[raw.assoc["id"] == sid][["rsid", "ea", "nea", "eaf", "beta"]].copy()
    flip = ref["beta"] < 0
    ref.loc[flip, ["ea", "nea"]] = ref.loc[flip, ["nea", "ea"]].values
    ref.loc[flip, "eaf"] = 1 - ref.loc[flip, "eaf"]
    ref = ref.rename(columns={"ea": "ref_ea", "nea": "ref_nea", "eaf": "ref_eaf"})[["rsid", "ref_ea", "ref_nea", "ref_eaf"]]

    key_of = dict(zip(raw.manifest["id"], raw.manifest["key"]))
    a = raw.assoc.merge(ref, on="rsid", how="inner")
    a["key"] = a["id"].map(key_of)
    # a study id can serve several keys (e.g. the same FinnGen file) -> expand
    multi = raw.manifest.dropna(subset=["id"]).groupby("id")["key"].apply(list)
    a = a.drop(columns="key").merge(multi.rename("key").reset_index().explode("key"), on="id")

    same = (a["ea"] == a["ref_ea"]) & (a["nea"] == a["ref_nea"])
    swap = (a["ea"] == a["ref_nea"]) & (a["nea"] == a["ref_ea"])
    ea_c, nea_c = a["ea"].map(COMPLEMENT), a["nea"].map(COMPLEMENT)
    same_c = (ea_c == a["ref_ea"]) & (nea_c == a["ref_nea"]) & ~same & ~swap
    swap_c = (ea_c == a["ref_nea"]) & (nea_c == a["ref_ea"]) & ~same & ~swap
    sign = np.select([same | same_c, swap | swap_c], [1.0, -1.0], default=np.nan)

    pal = _is_palindromic(a["ref_ea"], a["ref_nea"])
    eaf = a["eaf"].where(sign == 1, 1 - a["eaf"])
    ambiguous = pal & (a["ref_eaf"].between(maf_ambiguous, 1 - maf_ambiguous) | eaf.between(maf_ambiguous, 1 - maf_ambiguous))
    # palindromic and informative: decide orientation by allele frequency
    pal_ok = pal & ~ambiguous & a["eaf"].notna() & a["ref_eaf"].notna()
    pal_sign = np.where((a["eaf"] > 0.5) == (a["ref_eaf"] > 0.5), 1.0, -1.0)
    sign = np.where(pal_ok, pal_sign, sign)
    sign = np.where(ambiguous, np.nan, sign)

    a["beta"] = a["beta"] * sign
    a["eaf"] = np.where(sign == 1, a["eaf"], 1 - a["eaf"])
    a = a[np.isfinite(sign)]
    out = a[["rsid", "key", "id", "beta", "se", "p", "n", "eaf"] + (["proxy"] if "proxy" in a else [])].copy()
    out["z"] = out["beta"] / out["se"]
    return out.reset_index(drop=True)


def wide(h: pd.DataFrame, value: str, keys: list[str], rsids=None) -> pd.DataFrame:
    w = h[h["key"].isin(keys)].pivot_table(index="rsid", columns="key", values=value, aggfunc="first")
    w = w.reindex(columns=[k for k in keys if k in w.columns])
    if rsids is not None:
        w = w.reindex(rsids)
    return w
