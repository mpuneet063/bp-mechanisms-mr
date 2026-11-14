#!/usr/bin/env python3
"""Step 7: figures (PNG + PDF) and a markdown summary of every headline number."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402

INK, INK2, MUTED, SURF = "#0b0b0b", "#52514e", "#8a8984", "#fcfcfb"
BLUE, RED, NEUTRAL = "#2a78d6", "#e34948", "#f0efec"
DIVERGING = LinearSegmentedColormap.from_list("div", ["#104281", BLUE, NEUTRAL, RED, "#9e1f1f"])
SEQ = LinearSegmentedColormap.from_list("seq", [SURF, "#cde2fb", "#6da7ec", BLUE, "#104281"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 200})


from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator  # noqa: E402

OR_TICKS = [0.25, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 3, 4, 6, 8]


def or_axis(ax):
    """Log-scaled odds-ratio axis with plain-number ticks."""
    lo, hi = ax.get_xlim()
    ticks = [t for t in OR_TICKS if lo <= t <= hi] or [1]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig_heatmap(s, out, part, meta):
    prof = pd.read_csv(s.outp / "02_cluster" / f"profile_{part}.csv", index_col=0)
    prof.index = [f"c{c}" for c in prof.index]
    sizes = meta["sizes"][part]
    lim = np.nanpercentile(np.abs(prof.values), 98) or 1
    fig, ax = plt.subplots(figsize=(0.42 * prof.shape[1] + 2.5, 0.42 * prof.shape[0] + 1.4))
    im = ax.imshow(prof.values, cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(prof.shape[1]), prof.columns, rotation=60, ha="right")
    ax.set_yticks(range(prof.shape[0]), [f"{c} (n={sizes[c[1:]]})" for c in prof.index])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-.5, prof.shape[1]), minor=True)
    ax.set_yticks(np.arange(-.5, prof.shape[0]), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("mean trait z per unit SBP z" if s.profile == "per_sbp" else "mean scaled z (SBP-raising allele)")
    cb.outline.set_visible(False)
    ax.set_title(f"Mechanism clusters of BP variants ({part}): mean trait profile (SBP-raising allele)", loc="left", color=INK)
    save(fig, out, f"fig1_cluster_heatmap_{part}")


def fig_stability(s, out, meta):
    st = pd.read_csv(s.outp / "03_stability" / "hclust_stability_summary.csv")
    J = pd.read_csv(s.outp / "03_stability" / "hclust_bootstrap_jaccard.csv")
    kscan = pd.read_csv(s.outp / "03_stability" / "hclust_stability_by_k.csv")
    kc = pd.read_csv(s.outp / "02_cluster" / "bnmf_K_counts.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    ax = axes[0]
    y = np.arange(len(st))[::-1]
    lo, hi = J.quantile(0.05), J.quantile(0.95)
    for yi, r in zip(y, st.itertuples()):
        ax.plot([lo[r.cluster], hi[r.cluster]], [yi, yi], color=MUTED, lw=2, solid_capstyle="round")
        ax.plot(r.jaccard_mean, yi, "o", ms=8, color=BLUE if r.goes_to_mr else MUTED, mec="white", mew=1.5)
    ax.axvline(s.jaccard_min, color=INK2, lw=1, ls=":")
    ax.axvline(s.jaccard_stable, color=INK2, lw=1, ls="--")
    ax.text(s.jaccard_min, -0.9, f" gate {s.jaccard_min}", color=INK2, fontsize=7, va="center")
    ax.text(s.jaccard_stable, -0.9, f" stable {s.jaccard_stable}", color=INK2, fontsize=7, va="center")
    ax.set_ylim(-1.3, y.max() + 0.5)
    ax.set_yticks(y, [f"{r.cluster} (n={r.n_variants})" for r in st.itertuples()])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("bootstrap Jaccard (mean, 5th-95th pct)")
    ax.set_title("A  Per-cluster stability (Ward)", loc="left", color=INK)
    ax = axes[1]
    ax.plot(kscan["k"], kscan["mean_jaccard"], "-o", color=BLUE, lw=2, ms=6, label="mean over clusters")
    ax.plot(kscan["k"], kscan["min_cluster_jaccard"], "-o", color=MUTED, lw=2, ms=6, label="least stable cluster")
    ax.axvline(meta["k"], color=INK2, lw=1, ls=":")
    ax.set_xlabel("number of clusters K")
    ax.set_ylabel("bootstrap Jaccard")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_title("B  Stability across K", loc="left", color=INK)
    ax = axes[2]
    ax.bar(kc["K"], kc["runs"], color=BLUE, width=0.6)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("K selected by bNMF (ARD)")
    ax.set_ylabel("runs")
    ax.set_title("C  bNMF choice of K", loc="left", color=INK)
    fig.tight_layout()
    save(fig, out, "fig2_stability")

    co = np.load(s.outp / "03_stability" / "hclust_noise_coassign.npy")
    lab = pd.read_csv(s.outp / "02_cluster" / "labels.csv", index_col=0)["hclust"].to_numpy()
    o = np.argsort(lab, kind="stable")
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    im = ax.imshow(co[np.ix_(o, o)], cmap=SEQ, vmin=0, vmax=1, interpolation="nearest")
    b = np.cumsum(np.bincount(lab))[:-1]
    for x in b:
        ax.axhline(x - .5, color="white", lw=1.5)
        ax.axvline(x - .5, color="white", lw=1.5)
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("P(same cluster) under GWAS sampling noise")
    cb.outline.set_visible(False)
    ax.set_title("Co-assignment matrix, variants ordered by cluster", loc="left", color=INK, fontsize=9)
    save(fig, out, "fig2b_coassignment")


def fig_forest(s, out, part, info, meta):
    mr = pd.read_csv(s.outp / "04_mr" / "mr_results.csv")
    d = mr[(mr["method"] == "IVW") & (mr["partition"].isin([part, "hclust"]))]
    d = d[(d["partition"] == part) | (d["cluster"] == "all")]
    outs = [o for o in config.MAIN_OUTCOMES if o in d["outcome"].unique()]
    clusters = ["all"] + sorted([c for c in d["cluster"].unique() if c != "all"], key=lambda c: int(c[1:]))
    ncol = min(4, len(outs))
    nrow = int(np.ceil(len(outs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol + 1.8, 0.33 * len(clusters) * nrow + 1.2 * nrow),
                             sharey=True, squeeze=False)
    xs = d[["OR_lo", "OR_hi"]].replace([np.inf, -np.inf], np.nan)
    xmin, xmax = max(0.25, np.nanmin(xs.values) * 0.9), min(8, np.nanmax(xs.values) * 1.1)
    y = np.arange(len(clusters))[::-1]
    for i, o in enumerate(outs):
        ax = axes[i // ncol][i % ncol]
        sub = d[d["outcome"] == o].set_index("cluster")
        for yi, c in zip(y, clusters):
            if c not in sub.index:
                continue
            r = sub.loc[c]
            col = INK if c == "all" else (BLUE if r["passes_stability"] else MUTED)
            ax.plot([r["OR_lo"], r["OR_hi"]], [yi, yi], color=col, lw=2, solid_capstyle="round")
            ax.plot(r["OR"], yi, "s" if c == "all" else "o", ms=7, color=col, mec="white", mew=1.2)
        ax.axvline(1, color=MUTED, lw=1)
        ax.set_xscale("log")
        ax.set_xlim(xmin, xmax)
        or_axis(ax)
        ax.set_title(config.OUTCOME_LABELS[o], loc="left", color=INK, fontsize=9)
        ax.grid(axis="x", color="#e8e7e3", lw=0.6)
    for j in range(len(outs), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    names = meta["names"][part]
    labels = ["All BP variants"] + [f"{c}: {names[c[1:]]}" for c in clusters[1:]]
    for r in range(nrow):
        axes[r][0].set_yticks(y, labels)
    fig.supxlabel(f"Odds ratio {info['label']} (IVW, 95% CI, log scale)", color=INK2, fontsize=9)
    fig.suptitle(f"Effect of SBP on stroke subtypes and CAD, by mechanism cluster ({part}); grey = failed stability gate",
                 x=0.01, ha="left", color=INK, fontsize=10)
    fig.tight_layout()
    save(fig, out, f"fig3_forest_{part}")


def fig_mvmr(s, out):
    p = s.outp / "06_mvmr" / "univariable_vs_mvmr.csv"
    if not p.exists():
        return
    d = pd.read_csv(p)
    d = d[d["outcome"].isin(["cad", "stroke_svs", "stroke_las", "stroke_any"])]
    if d.empty:
        return
    outs = [o for o in ["stroke_any", "stroke_las", "stroke_svs", "cad"] if o in d["outcome"].unique()]
    fig, axes = plt.subplots(1, len(outs), figsize=(3.2 * len(outs) + 1.5, 0.4 * d["cluster"].nunique() + 1.4), sharey=True)
    axes = np.atleast_1d(axes)
    cl = sorted(d["cluster"].unique(), key=lambda c: int(c[1:]))
    y = np.arange(len(cl))[::-1]
    for ax, o in zip(axes, outs):
        sub = d[d["outcome"] == o].set_index("cluster")
        for yi, c in zip(y, cl):
            if c not in sub.index:
                continue
            r = sub.loc[c]
            ax.plot([r["OR_univariable"], r["OR_mvmr"]], [yi, yi], color="#d8d7d2", lw=3, solid_capstyle="round")
            ax.plot(r["OR_univariable"], yi, "o", ms=7, color=MUTED, mec="white", mew=1.2)
            ax.plot(r["OR_mvmr"], yi, "o", ms=7, color=BLUE, mec="white", mew=1.2)
        ax.axvline(1, color=MUTED, lw=1)
        ax.set_xscale("log")
        or_axis(ax)
        ax.set_title(config.OUTCOME_LABELS[o], loc="left", color=INK, fontsize=9)
    axes[0].set_yticks(y, [f"{c}: {sub['name'].get(c, '') if c in sub.index else ''}" for c in cl])
    axes[0].plot([], [], "o", color=MUTED, label="univariable IVW")
    axes[0].plot([], [], "o", color=BLUE, label="MVMR (adj. BMI, LDL, T2D)")
    fig.legend(loc="upper right", frameon=False, ncol=2, fontsize=8)
    fig.supxlabel("Odds ratio per 10 mmHg SBP (log scale)", color=INK2, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, out, "fig4_mvmr")


def summary(s, out, info, meta):
    st = pd.read_csv(s.outp / "03_stability" / "hclust_stability_summary.csv")
    mr = pd.read_csv(s.outp / "04_mr" / "mr_results.csv")
    het = pd.read_csv(s.outp / "05_compare" / "between_cluster_heterogeneity.csv")
    pw = pd.read_csv(s.outp / "05_compare" / "pairwise_contrasts.csv")
    part = s.primary
    L = [f"# Results summary ({'SIMULATED DATA' if 'sim' in s.raw else 'real data'})", "",
         f"Exposure scaling: {info['label']}. Clustering variants: {info['n_cluster_variants']}; "
         f"MR instruments: {info['n_mr_instruments']}. Traits: {', '.join(info['traits'])}.", "",
         f"## Clusters (K = {meta['k']}, bNMF K counts {meta['k_counts']}; hclust vs bNMF ARI = {meta['ari']:.2f})", "",
         "| cluster | n | top traits | bootstrap Jaccard | noise keep prob | stability | to MR |", "|---|---|---|---|---|---|---|"]
    for r in st.itertuples():
        L.append(f"| {r.cluster} | {r.n_variants} | {r.name} | {r.jaccard_mean:.2f} | {r.noise_keep_prob:.2f} | {r.stability} | {'yes' if r.goes_to_mr else 'no'} |")
    L += ["", f"## IVW odds ratios {info['label']} ({part})", ""]
    iv = mr[(mr["method"] == "IVW") & ((mr["partition"] == part) | (mr["cluster"] == "all"))]
    outs = [o for o in s.outcomes if o in iv["outcome"].unique()]
    L.append("| cluster | " + " | ".join(config.OUTCOME_LABELS[o] for o in outs) + " |")
    L.append("|---" * (len(outs) + 1) + "|")
    for c in ["all"] + sorted([c for c in iv["cluster"].unique() if c != "all"], key=lambda c: int(c[1:])):
        cells = []
        for o in outs:
            r = iv[(iv["cluster"] == c) & (iv["outcome"] == o)]
            cells.append("" if r.empty else f"{r['OR'].iloc[0]:.2f} ({r['OR_lo'].iloc[0]:.2f}-{r['OR_hi'].iloc[0]:.2f})")
        L.append(f"| {c} | " + " | ".join(cells) + " |")
    L += ["", "## Do clusters differ? (stable clusters only)", "", "| outcome | Q_between | df | p (chi2) | p (permutation) |", "|---|---|---|---|---|"]
    for r in het[het["partition"] == part].itertuples():
        L.append(f"| {config.OUTCOME_LABELS.get(r.outcome, r.outcome)} | {r.Q_between:.1f} | {r.Q_df} | {r.Q_p:.2g} | {r.perm_p:.3g} |")
    top = pw[(pw["partition"] == part) & pw["ratio_interpretable"]].sort_values("diff_p").head(8)
    L += ["", "## Strongest pairwise contrasts", ""]
    for r in top.itertuples():
        L.append(f"- {config.OUTCOME_LABELS.get(r.outcome, r.outcome)}: c{r.cluster_a} ({r.name_a}) vs c{r.cluster_b} ({r.name_b}): "
                 f"log-OR ratio {r.ratio_logOR:.2f} (bootstrap 95% CI {r.ratio_lo:.2f} to {r.ratio_hi:.2f}); "
                 f"OR ratio {r.ratio_OR_scaled:.2f} {info['label']}; difference p = {r.diff_p:.2g}")
    egg = mr[(mr["method"] == "MR-Egger") & (mr["partition"] == part)]
    L += ["", "## Pleiotropy checks", "",
          f"- MR-Egger intercept p < 0.05 in {int((egg['int_p'] < 0.05).sum())} of {len(egg)} cluster x outcome analyses.",
          f"- MR-PRESSO global test p < 0.05 in {int((mr.loc[mr['method'].str.startswith('MR-PRESSO') & (mr['partition'] == part), 'global_p'] < 0.05).sum())} analyses."]
    (out / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    s = config.parse(__doc__)
    out = s.step("07_report")
    info = json.loads((s.outp / "01_prepare" / "info.json").read_text())
    meta = json.loads((s.outp / "02_cluster" / "clusters.json").read_text())
    for part in ["hclust", "bnmf"]:
        fig_heatmap(s, out, part, meta)
        fig_forest(s, out, part, info, meta)
    fig_stability(s, out, meta)
    fig_mvmr(s, out)
    summary(s, out, info, meta)


if __name__ == "__main__":
    main()
