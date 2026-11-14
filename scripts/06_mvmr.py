#!/usr/bin/env python3
"""
Step 6: multivariable MR. Each stable cluster's SBP effect is a separate exposure (SBP association for that
cluster's variants, 0 elsewhere), fitted jointly with BMI, LDL and T2D liability. This separates "SBP acting
through mechanism k" from the cluster's other (pleiotropic) associations.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr.compare import conditional_f, mvmr_ivw  # noqa: E402
from bpmr.data import wide  # noqa: E402


def main():
    s = config.parse(__doc__)
    prep, cl, stab, mrd, out = (s.outp / "01_prepare", s.outp / "02_cluster", s.outp / "03_stability",
                                s.outp / "04_mr", s.step("06_mvmr"))
    h = pd.read_csv(prep / "harmonised.csv.gz")
    var = pd.read_csv(prep / "variants.csv", index_col=0)
    info = json.loads((prep / "info.json").read_text())
    labels = pd.read_csv(cl / "labels.csv", index_col=0)
    meta = json.loads((cl / "clusters.json").read_text())
    mr = pd.read_csv(mrd / "mr_results.csv")
    inst_tab = pd.read_csv(Path(s.raw) / "instruments.csv")
    part = s.primary
    gate = pd.read_csv(stab / f"{part}_stability_summary.csv").set_index("cluster")["goes_to_mr"]
    ok = [int(c[1:]) for c, g in gate.items() if g]
    mediators = [m for m in config.MVMR_MEDIATORS if m in set(h["key"])]

    bp = var[var["mr_instrument"]].join(labels[[part]], how="inner")
    bp = bp[bp[part].isin(ok)]
    med_snps = inst_tab.loc[inst_tab["source"].isin([f"mvmr_{m}" for m in mediators]), "rsid"].unique().tolist()
    snps = sorted(set(bp.index) | set(med_snps))
    B = wide(h, "beta", ["sbp"] + mediators, snps)
    S = wide(h, "se", ["sbp"] + mediators, snps)
    B, S = B.dropna(), S.loc[B.dropna().index]
    snps = B.index

    cols, BX, SX = [], [], []
    member = labels[part].reindex(snps)
    for c in ok:
        m = (member == c).to_numpy() & np.isin(snps, bp.index)
        BX.append(np.where(m, B["sbp"], 0.0)), SX.append(np.where(m, S["sbp"], 0.0)), cols.append(f"SBP via c{c}")
    other = ~np.isin(snps, bp.index)
    BX.append(np.where(other, B["sbp"], 0.0)), SX.append(np.where(other, S["sbp"], 0.0)), cols.append("SBP (other variants)")
    for m in mediators:
        BX.append(B[m].to_numpy()), SX.append(S[m].to_numpy()), cols.append(m.upper())
    BX, SX = np.column_stack(BX), np.column_stack(SX)
    cf = conditional_f(BX, SX, cols)
    cf.to_csv(out / "conditional_F.csv")

    rows = []
    for o in mr["outcome"].unique():
        y = h[h["key"] == o].drop_duplicates("rsid").set_index("rsid").reindex(snps)
        keep = y["beta"].notna().to_numpy()
        res = mvmr_ivw(BX[keep], y["beta"].to_numpy()[keep], y["se"].to_numpy()[keep], cols)
        res.insert(0, "outcome", o)
        rows.append(res)
    mv = pd.concat(rows, ignore_index=True)
    sbp_rows = mv["exposure"].str.startswith("SBP")
    mv["b_scaled"] = np.where(sbp_rows, mv["b"] * info["scale"], mv["b"])
    mv["se_scaled"] = np.where(sbp_rows, mv["se"] * info["scale"], mv["se"])
    mv["OR"] = np.exp(mv["b_scaled"])
    mv["OR_lo"] = np.exp(mv["b_scaled"] - 1.96 * mv["se_scaled"])
    mv["OR_hi"] = np.exp(mv["b_scaled"] + 1.96 * mv["se_scaled"])
    mv["name"] = mv["exposure"].map(lambda e: meta["names"][part].get(e.split("c")[-1], "") if "via c" in e else "")
    mv.to_csv(out / "mvmr_results.csv", index=False)

    # side-by-side: univariable vs multivariable for each cluster
    uni = mr[(mr["partition"] == part) & (mr["method"] == "IVW")].set_index(["cluster", "outcome"])
    comp = []
    for r in mv[mv["exposure"].str.contains("via c")].itertuples():
        c = r.exposure.split("via ")[1]
        if (c, r.outcome) in uni.index:
            u = uni.loc[(c, r.outcome)]
            comp.append(dict(cluster=c, name=r.name, outcome=r.outcome, OR_univariable=u["OR"], OR_mvmr=r.OR,
                             OR_mvmr_lo=r.OR_lo, OR_mvmr_hi=r.OR_hi, p_mvmr=r.p))
    pd.DataFrame(comp).to_csv(out / "univariable_vs_mvmr.csv", index=False)
    print("Conditional F:", cf.round(1).to_dict())
    print(pd.DataFrame(comp).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
