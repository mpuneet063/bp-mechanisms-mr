#!/usr/bin/env python3
"""Step 4: Mendelian randomisation of SBP on each outcome, for all BP variants and for each mechanism cluster."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr.data import wide  # noqa: E402
from bpmr.mr import run_all, steiger_keep, ivw  # noqa: E402


def outcome_table(h, key, rsids):
    t = h[h["key"] == key].drop_duplicates("rsid").set_index("rsid")
    return t.reindex(rsids)[["beta", "se", "eaf"]].rename(columns={"beta": "by", "se": "sy", "eaf": "eaf_y"})


def main():
    s = config.parse(__doc__)
    prep, cl, stab, out = s.outp / "01_prepare", s.outp / "02_cluster", s.outp / "03_stability", s.step("04_mr")
    h = pd.read_csv(prep / "harmonised.csv.gz")
    var = pd.read_csv(prep / "variants.csv", index_col=0)
    info = json.loads((prep / "info.json").read_text())
    labels = pd.read_csv(cl / "labels.csv", index_col=0)
    meta = json.loads((cl / "clusters.json").read_text())
    man = pd.read_csv(prep / "manifest_used.csv").set_index("key")
    scale = info["scale"]

    inst = var[var["mr_instrument"]].join(labels[["hclust", "bnmf"]], how="inner")
    outcomes = [o for o in s.outcomes if o in man.index and isinstance(man.loc[o, "id"], str)]
    rows, steiger_rows = [], []
    for partition in ["hclust", "bnmf"]:
        gate = pd.read_csv(stab / f"{partition}_stability_summary.csv").set_index("cluster")["goes_to_mr"]
        groups = [("all", inst.index)] + [(f"c{c}", inst.index[inst[partition] == c]) for c in sorted(inst[partition].unique())]
        for o in outcomes:
            oy = outcome_table(h, o, inst.index)
            d = inst.join(oy).dropna(subset=["by", "sy"])
            binary = pd.notna(man.loc[o, "ncase"]) and float(man.loc[o, "ncase"] or 0) > 0
            for g, idx in groups:
                if partition == "bnmf" and g == "all":
                    continue  # identical to the hclust "all" row
                dd = d.loc[d.index.intersection(idx)]
                if len(dd) < 3:
                    continue
                bx, sx, by, sy = (dd[c].to_numpy(float) for c in ["bx", "sx", "by", "sy"])
                res = run_all(bx, sx, by, sy, nboot=s.wm_boot, nsim=s.presso_sim, seed=s.seed)
                # Steiger-filtered IVW (sensitivity)
                keep = steiger_keep(bx, sx, dd["nx"].to_numpy(float), by, sy, dd["eaf"].fillna(dd["eaf_y"]).fillna(0.3).to_numpy(float),
                                    float(man.loc[o, "ncase"] or 0), float(man.loc[o, "ncontrol"] or 0), binary)
                if keep.sum() >= 3:
                    r = ivw(bx[keep], sx[keep], by[keep], sy[keep])
                    r["method"] = "IVW (Steiger-filtered)"
                    res = pd.concat([res, pd.DataFrame([r])], ignore_index=True)
                steiger_rows.append(dict(partition=partition, cluster=g, outcome=o, n=len(keep), steiger_removed=int((~keep).sum())))
                res.insert(0, "outcome", o)
                res.insert(0, "cluster", g)
                res.insert(0, "partition", partition)
                res["name"] = "all BP variants" if g == "all" else meta["names"][partition].get(g[1:], "")
                res["passes_stability"] = True if g == "all" else bool(gate.get(g, False))
                rows.append(res)
    mr = pd.concat(rows, ignore_index=True)
    mr["b_scaled"] = mr["b"] * scale
    mr["se_scaled"] = mr["se"] * scale
    mr["OR"] = np.exp(mr["b_scaled"])
    mr["OR_lo"] = np.exp(mr["b_scaled"] - 1.96 * mr["se_scaled"])
    mr["OR_hi"] = np.exp(mr["b_scaled"] + 1.96 * mr["se_scaled"])
    mr["scale_label"] = info["label"]
    mr.drop(columns=[c for c in ["outliers"] if c in mr]).to_csv(out / "mr_results.csv", index=False)
    pd.DataFrame(steiger_rows).to_csv(out / "steiger.csv", index=False)

    main_ivw = mr[(mr["method"] == "IVW") & (mr["partition"] == s.primary)]
    piv = main_ivw.pivot_table(index=["cluster", "name"], columns="outcome", values="OR")
    piv = piv[[o for o in outcomes if o in piv.columns]]
    print(f"IVW odds ratios {info['label']} ({s.primary} partition):")
    print(piv.round(2).to_string())
    piv.to_csv(out / f"ivw_OR_matrix_{s.primary}.csv")


if __name__ == "__main__":
    main()
