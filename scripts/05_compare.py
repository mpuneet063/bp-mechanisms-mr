#!/usr/bin/env python3
"""Step 5: do mechanism clusters differ in their effect on each outcome? Q_between, permutation null, pairwise."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr.compare import mr_clusters, pairwise, permutation_test, q_between  # noqa: E402


def main():
    s = config.parse(__doc__)
    prep, cl, stab, mrd, out = (s.outp / "01_prepare", s.outp / "02_cluster", s.outp / "03_stability",
                                s.outp / "04_mr", s.step("05_compare"))
    h = pd.read_csv(prep / "harmonised.csv.gz")
    var = pd.read_csv(prep / "variants.csv", index_col=0)
    info = json.loads((prep / "info.json").read_text())
    labels = pd.read_csv(cl / "labels.csv", index_col=0)
    meta = json.loads((cl / "clusters.json").read_text())
    mr = pd.read_csv(mrd / "mr_results.csv")
    scale = info["scale"]

    rows, pairs = [], []
    for partition in config.PARTITIONS:
        ss = pd.read_csv(stab / f"{partition}_stability_summary.csv").set_index("cluster")
        ok, exploratory = mr_clusters(ss["goes_to_mr"], ss["n_variants"], s.min_cluster_snps)
        inst = var[var["mr_instrument"]].join(labels[[partition]], how="inner")
        inst = inst[inst[partition].isin(ok)]
        names = {c: meta["names"][partition][str(c)] for c in ok}
        for o in mr["outcome"].unique():
            y = h[h["key"] == o].drop_duplicates("rsid").set_index("rsid")[["beta", "se"]]
            d = inst.join(y, how="inner").dropna(subset=["beta", "se"])
            if d[partition].nunique() < 2:
                continue
            bx, sx, by, sy = (d[c].to_numpy(float) for c in ["bx", "sx", "beta", "se"])
            lab = d[partition].to_numpy()
            sub = mr[(mr["partition"] == partition) & (mr["outcome"] == o) & (mr["method"] == "IVW")
                     & mr["cluster"].isin([f"c{c}" for c in ok])]
            Q, Qp = q_between(sub["b"].to_numpy(), sub["se"].to_numpy())
            perm = permutation_test(bx, sx, by, sy, lab, B=s.perm_B, seed=s.seed)
            rows.append(dict(partition=partition, outcome=o, n_clusters=len(sub), Q_between=Q, Q_df=len(sub) - 1,
                             Q_p=Qp, perm_p=perm["perm_p"], perm_null_median_Q=perm["null_median"],
                             exploratory=exploratory))
            pw = pairwise(bx, sx, by, sy, lab, names, B=s.pair_B, seed=s.seed)
            pw.insert(0, "outcome", o)
            pw.insert(0, "partition", partition)
            pw["exploratory"] = exploratory
            pairs.append(pw)
    het = pd.DataFrame(rows)
    het.to_csv(out / "between_cluster_heterogeneity.csv", index=False)
    pw = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame(
        columns=["partition", "outcome", "name_a", "name_b", "diff_logOR", "diff_p", "ratio_logOR", "ratio_lo", "ratio_hi",
                 "ratio_interpretable", "exploratory"])
    pw["diff_logOR_scaled"] = pw["diff_logOR"] * scale
    pw["ratio_OR_scaled"] = np.exp(pw["diff_logOR_scaled"])
    pw.to_csv(out / "pairwise_contrasts.csv", index=False)

    # stroke-vs-CAD specificity: does a cluster's stroke effect differ from its CAD effect more than average?
    spec_rows = []
    for partition in config.PARTITIONS:
        iv = mr[(mr["partition"] == partition) & (mr["method"] == "IVW")]
        for o in [x for x in iv["outcome"].unique() if x.startswith("stroke") or x == "ich"]:
            a = iv[iv["outcome"] == o].set_index("cluster")
            c = iv[iv["outcome"] == "cad"].set_index("cluster")
            for cl_ in a.index.intersection(c.index):
                d = a.loc[cl_, "b"] - c.loc[cl_, "b"]
                se = np.sqrt(a.loc[cl_, "se"] ** 2 + c.loc[cl_, "se"] ** 2)  # outcome GWAS treated as independent
                spec_rows.append(dict(partition=partition, cluster=cl_, name=a.loc[cl_, "name"], stroke_outcome=o,
                                      logOR_stroke=a.loc[cl_, "b"] * scale, logOR_cad=c.loc[cl_, "b"] * scale,
                                      diff=d * scale, diff_se=se * scale))
    pd.DataFrame(spec_rows).to_csv(out / "stroke_vs_cad_by_cluster.csv", index=False)

    if het.empty:
        print("no partition had two or more clusters to compare")
    else:
        print(het[het["partition"] == s.primary].round(4).to_string(index=False))
    top = pw[(pw["partition"] == s.primary) & pw["ratio_interpretable"]].sort_values("diff_p").head(10)
    print("\nStrongest pairwise contrasts (" + info["label"] + "):")
    for r in top.itertuples():
        print(f"  {r.outcome}: {r.name_a} vs {r.name_b}: logOR ratio {r.ratio_logOR:.2f} "
              f"(95% CI {r.ratio_lo:.2f}-{r.ratio_hi:.2f}), diff p = {r.diff_p:.2g}")


if __name__ == "__main__":
    main()
