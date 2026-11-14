"""
Cluster stability: the part the February 2026 paper does not report.

* variant bootstrap  (Hennig 2007 clusterboot): resample variants, recluster, best-match Jaccard per cluster
* noise bootstrap    : redraw every z-score from its sampling distribution, recluster all variants,
                       per-variant assignment probability + co-assignment matrix
* leave-one-trait-out and K +/- 1: adjusted Rand index against the main clustering
Thresholds (fixed before any MR is run): Jaccard >= 0.75 stable, 0.6-0.75 moderately stable, < 0.6 unstable.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

from .cluster import align_labels, ari


def jaccard_best(orig_members: np.ndarray, boot_labels: np.ndarray, boot_idx: np.ndarray) -> float:
    """Best Jaccard between an original cluster (restricted to sampled variants) and any bootstrap cluster."""
    uniq, first = np.unique(boot_idx, return_index=True)
    lab = boot_labels[first]
    a = set(uniq[orig_members[uniq]])
    if not a:
        return np.nan
    best = 0.0
    for c in np.unique(lab):
        b = set(uniq[lab == c])
        best = max(best, len(a & b) / len(a | b))
    return best


def variant_bootstrap(x: np.ndarray, labels: np.ndarray, cluster_fn: Callable[[np.ndarray], np.ndarray],
                      B: int = 500, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    k = labels.max() + 1
    J = np.full((B, k), np.nan)
    for b in range(B):
        idx = rng.integers(0, len(x), len(x))
        bl = cluster_fn(x[idx])
        for c in range(k):
            J[b, c] = jaccard_best(labels == c, bl, idx)
    return pd.DataFrame(J, columns=[f"c{c}" for c in range(k)])


def noise_bootstrap(x: np.ndarray, noise_sd: np.ndarray, labels: np.ndarray,
                    cluster_fn: Callable[[np.ndarray], np.ndarray], B: int = 200, seed: int = 2):
    """Returns (per-variant probability of keeping its cluster, co-assignment matrix, per-run ARI)."""
    rng = np.random.default_rng(seed)
    n, k = len(x), labels.max() + 1
    keep = np.zeros(n)
    co = np.zeros((n, n), dtype=np.float32)
    aris = []
    for _ in range(B):
        xb = x + rng.normal(0, 1, x.shape) * noise_sd
        bl = align_labels(labels, cluster_fn(xb), k_ref=k)
        keep += bl == labels
        co += (bl[:, None] == bl[None, :])
        aris.append(ari(labels, bl))
    return keep / B, co / B, np.array(aris)


def summarise(labels: np.ndarray, J: pd.DataFrame, keep_prob: np.ndarray, co: np.ndarray) -> pd.DataFrame:
    rows = []
    for c in range(labels.max() + 1):
        m = labels == c
        within = co[np.ix_(m, m)]
        rows.append(dict(cluster=f"c{c}", n_variants=int(m.sum()),
                         jaccard_mean=J[f"c{c}"].mean(), jaccard_q05=J[f"c{c}"].quantile(0.05),
                         dissolved_pct=(J[f"c{c}"] < 0.5).mean() * 100,
                         noise_keep_prob=keep_prob[m].mean(),
                         coassign_within=within[np.triu_indices(m.sum(), 1)].mean() if m.sum() > 1 else np.nan))
    s = pd.DataFrame(rows)
    s["stability"] = pd.cut(s["jaccard_mean"], [-1, 0.6, 0.75, 1.01], labels=["unstable", "moderate", "stable"])
    return s


def leave_one_trait_out(z: pd.DataFrame, labels: np.ndarray, make_x: Callable, cluster_fn: Callable) -> pd.Series:
    out = {}
    for t in z.columns:
        if t == "sbp":
            continue
        out[t] = ari(labels, cluster_fn(make_x(z.drop(columns=t))))
    return pd.Series(out, name="ARI_without_trait").sort_values()


def k_sensitivity(x: np.ndarray, labels: np.ndarray, fn_k: Callable[[np.ndarray, int], np.ndarray], k: int) -> pd.Series:
    return pd.Series({kk: ari(labels, fn_k(x, kk)) for kk in [k - 1, k + 1] if kk >= 2}, name="ARI_vs_main")


def stability_by_k(x: np.ndarray, fn_k: Callable[[np.ndarray, int], np.ndarray], ks=range(2, 11), B: int = 100,
                   seed: int = 3) -> pd.DataFrame:
    """Mean bootstrap Jaccard (averaged over clusters) for each K: a data-driven check on the choice of K."""
    rows = []
    for k in ks:
        lab = fn_k(x, k)
        J = variant_bootstrap(x, lab, lambda xx: fn_k(xx, k), B=B, seed=seed)
        rows.append(dict(k=k, mean_jaccard=J.mean().mean(), min_cluster_jaccard=J.mean().min()))
    return pd.DataFrame(rows)
