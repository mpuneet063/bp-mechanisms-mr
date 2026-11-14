#!/usr/bin/env python3
"""
Score a pipeline run on SIMULATED data against the planted truth (data/sim_raw/truth.csv):
  * cluster recovery: adjusted Rand index of hclust / bNMF labels vs the planted mechanisms
  * MR recovery: each found cluster is matched to its majority planted mechanism; IVW log-OR vs planted theta
  * MVMR: does adjusting for BMI remove the planted BMI->CAD pleiotropy in the adiposity cluster?
Writes <out>/sim_validation.md.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr.cluster import ari  # noqa: E402
from bpmr.simulate import THETA  # noqa: E402


def main():
    s = config.parse(__doc__)
    truth = pd.read_csv(Path(s.raw) / "truth.csv", index_col=0)["mech"]
    labels = pd.read_csv(s.outp / "02_cluster" / "labels.csv", index_col=0)
    mr = pd.read_csv(s.outp / "04_mr" / "mr_results.csv")
    info = json.loads((s.outp / "01_prepare" / "info.json").read_text())
    t = truth.reindex(labels.index)
    L = ["# Simulation validation", "", f"Planted mechanisms: {t.value_counts().to_dict()}", ""]
    L += ["| partition | ARI vs truth | K found |", "|---|---|---|"]
    for p in ["hclust", "bnmf"]:
        L.append(f"| {p} | {ari(t.astype('category').cat.codes, labels[p]):.3f} | {labels[p].nunique()} |")

    rows = []
    for p in ["hclust", "bnmf"]:
        for c in sorted(labels[p].unique()):
            m = t[labels[p] == c]
            mech, purity = m.value_counts().idxmax(), m.value_counts(normalize=True).max()
            for o in config.MAIN_OUTCOMES:
                r = mr[(mr["partition"] == p) & (mr["cluster"] == f"c{c}") & (mr["outcome"] == o) & (mr["method"] == "IVW")]
                if r.empty:
                    continue
                est, se = r["b_scaled"].iloc[0], r["se_scaled"].iloc[0]
                # expected = mixture of planted thetas over the cluster's true members
                exp = np.mean([THETA[x][o] for x in m])
                rows.append(dict(partition=p, cluster=f"c{c}", majority=mech, purity=purity, outcome=o,
                                 est=est, se=se, expected=exp, z=(est - exp) / se))
    d = pd.DataFrame(rows)
    d.to_csv(s.outp / "sim_validation_mr.csv", index=False)
    for p in ["hclust", "bnmf"]:
        dd = d[d["partition"] == p]
        cover = (dd["z"].abs() < 1.96).mean()
        corr = np.corrcoef(dd["est"], dd["expected"])[0, 1]
        L += ["", f"## {p}: MR estimates vs planted effects ({info['label']})", "",
              f"- correlation(estimate, planted) = {corr:.3f}; 95% CI covers the planted value in {cover:.0%} of "
              f"{len(dd)} cluster x outcome estimates (expect ~95%, lower for the adiposity x CAD cells with planted pleiotropy)."]
        for c, g in dd.groupby("cluster"):
            L.append(f"- {c}: majority '{g['majority'].iloc[0]}' (purity {g['purity'].iloc[0]:.0%})")
    mv = s.outp / "06_mvmr" / "univariable_vs_mvmr.csv"
    if mv.exists():
        u = pd.read_csv(mv)
        adip = d[(d["partition"] == s.primary) & (d["majority"] == "adiposity")]["cluster"].unique()
        L += ["", "## MVMR check (planted BMI->CAD pleiotropy in the adiposity cluster)", ""]
        for c in adip:
            r = u[(u["cluster"] == c) & (u["outcome"] == "cad")]
            if not r.empty:
                L.append(f"- {c}: CAD OR univariable {r['OR_univariable'].iloc[0]:.3f} -> MVMR {r['OR_mvmr'].iloc[0]:.3f} "
                         f"(planted SBP effect OR {np.exp(THETA['adiposity']['cad']):.3f})")
    (s.outp / "sim_validation.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
