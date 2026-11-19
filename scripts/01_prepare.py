#!/usr/bin/env python3
"""Step 1: harmonise all associations to the SBP-raising allele and build the variant x trait matrix."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr.data import harmonise, load_raw, wide  # noqa: E402
from bpmr.matrix import build_matrix, profile_matrix  # noqa: E402


SBP_SD_MMHG = 19.0  # population SD of SBP in UK Biobank-scale cohorts, used only if effects are SD-scaled


def infer_sbp_unit(bx: pd.Series, sx: pd.Series) -> str:
    """
    OpenGWAS leaves the unit blank for several BP studies (Evangelou 2018 included).
    Genome-wide significant SBP effects are ~0.2-1 mmHg per allele, i.e. ~0.01-0.05 SD (SD ~19 mmHg),
    so the median |beta| of the significant variants separates the two scales by an order of magnitude.
    """
    sig = (bx / sx).abs() > 5.45  # p < 5e-8
    med = bx[sig].abs().median() if sig.sum() >= 20 else np.nan
    if not np.isfinite(med):
        return "unknown"
    return "mmHg" if med > 0.08 else "SD"


def exposure_scale(manifest: pd.DataFrame, bx: pd.Series | None = None, sx: pd.Series | None = None):
    row = manifest[manifest["key"] == "sbp"].iloc[0]
    unit = str(row.get("unit", "")).lower()
    sd = pd.to_numeric(row.get("sd", np.nan), errors="coerce")
    if "mmhg" in unit:
        return 10.0, "per 10 mmHg higher SBP"
    if np.isfinite(sd) and sd > 0:
        return 10.0 / sd, f"per 10 mmHg higher SBP (converted from SD units, SD = {sd:g} mmHg)"
    if bx is not None and unit in ("", "nan", "none", "na"):
        guess = infer_sbp_unit(bx, sx)
        med = bx[(bx / sx).abs() > 5.45].abs().median()
        if guess == "mmHg":
            return 10.0, f"per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = {med:.3f})"
        if guess == "SD":
            return 10.0 / SBP_SD_MMHG, (f"per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred SD units from "
                                        f"median |beta| = {med:.3f}; assumed SD = {SBP_SD_MMHG:g} mmHg)")
    return 1.0, f"per unit of the SBP GWAS ({row.get('unit', 'unknown unit')})"


def main():
    s = config.parse(__doc__)
    out = s.step("01_prepare")
    raw = load_raw(s.raw)
    log = [f"raw: {s.raw}", f"{len(raw.assoc)} association rows, {raw.assoc['rsid'].nunique()} variants"]

    h = harmonise(raw)
    log.append(f"harmonised rows: {len(h)} ({len(raw.assoc) - len(h)} dropped: allele mismatch or ambiguous palindromes)")
    h.to_csv(out / "harmonised.csv.gz", index=False, compression="gzip")

    inst = raw.instruments
    bp_src = ~inst["source"].astype(str).str.startswith("mvmr")
    bp = inst[bp_src]
    clump05 = bp.groupby("rsid")["clump_r2_0.05"].any() if "clump_r2_0.05" in bp else pd.Series(True, index=bp["rsid"].unique())
    clump001 = bp.groupby("rsid")["clump_r2_0.001"].any() if "clump_r2_0.001" in bp else clump05
    cluster_set = sorted(clump05[clump05].index)
    zs = h[(h["key"] == "sbp")].drop_duplicates("rsid").set_index("rsid")["z"]
    weak = [r for r in cluster_set if not (abs(zs.get(r, 0.0)) >= s.min_sbp_z_cluster)]
    cluster_set = [r for r in cluster_set if r not in set(weak)]
    log.append(f"{len(weak)} variants dropped: |z_SBP| < {s.min_sbp_z_cluster} or absent in the SBP GWAS")
    log.append(f"BP variants: {bp['rsid'].nunique()} candidates; {len(cluster_set)} independent at r2<0.05 (clustering set)")

    traits = [k for k in raw.keys("exposure") + raw.keys("cluster_trait")]
    zm = build_matrix(h, cluster_set, traits, s.max_missing_variant, s.max_missing_trait, s.corr_cutoff)
    log += zm.log
    zm.z.to_csv(out / "z_scaled.csv")
    zm.z_raw.to_csv(out / "z_raw.csv")
    zm.noise_sd.to_csv(out / "noise_sd.csv")
    zc = profile_matrix(zm.z, s.profile)
    zc.to_csv(out / "z_cluster.csv")
    log.append(f"clustering input: {s.profile} ({zc.shape[1]} traits)")

    # exposure table for MR
    bx = wide(h, "beta", ["sbp"]).iloc[:, 0]
    sx = wide(h, "se", ["sbp"]).iloc[:, 0]
    nx = wide(h, "n", ["sbp"]).iloc[:, 0]
    eaf = wide(h, "eaf", ["sbp"]).iloc[:, 0]
    var = pd.DataFrame(dict(bx=bx, sx=sx, nx=nx, eaf=eaf)).dropna(subset=["bx", "sx"])
    var["F"] = (var["bx"] / var["sx"]) ** 2
    var["in_cluster_set"] = var.index.isin(zm.z.index)
    var["mr_independent"] = var.index.isin(clump001[clump001].index)
    var["mr_instrument"] = var["in_cluster_set"] & var["mr_independent"] & (var["F"] >= s.f_min)
    var.index.name = "rsid"
    var.to_csv(out / "variants.csv")
    log.append(f"MR instruments (in clustering set, r2<0.001, F>={s.f_min:g}): {int(var['mr_instrument'].sum())}")

    scale, label = exposure_scale(raw.manifest, var["bx"], var["sx"])
    info = dict(scale=scale, label=label, n_cluster_variants=int(len(zm.z)), traits=list(zm.z.columns),
                n_mr_instruments=int(var["mr_instrument"].sum()))
    (out / "info.json").write_text(json.dumps(info, indent=2))
    raw.manifest.to_csv(out / "manifest_used.csv", index=False)
    (out / "log.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))
    print(f"exposure scaling: {label}")


if __name__ == "__main__":
    main()
