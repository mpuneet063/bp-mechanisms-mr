"""Clustering: Ward hierarchical clustering (Pascat et al. 2026) and Bayesian NMF with ARD (Udler et al. 2018)."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import adjusted_rand_score, silhouette_score


# --------------------------------------------------------------------------- #
# Hierarchical clustering
# --------------------------------------------------------------------------- #
def hclust_input(z: pd.DataFrame, clip_sd: float | None = 2.0) -> np.ndarray:
    """Pascat et al.: Euclidean distance on z-scores truncated at +/- clip_sd standard deviations per trait."""
    x = z.to_numpy(dtype=float)
    if clip_sd is not None:
        mu, sd = x.mean(0), x.std(0) + 1e-12
        x = np.clip(x, mu - clip_sd * sd, mu + clip_sd * sd)
    return x


def ward(x: np.ndarray, k: int) -> np.ndarray:
    lab = fcluster(linkage(x, method="ward", metric="euclidean"), t=k, criterion="maxclust")
    return lab - 1


def silhouette_by_k(x: np.ndarray, ks=range(2, 11)) -> pd.Series:
    Z = linkage(x, method="ward")
    out = {}
    for k in ks:
        lab = fcluster(Z, t=k, criterion="maxclust")
        if len(np.unique(lab)) > 1:
            out[k] = silhouette_score(x, lab)
    return pd.Series(out, name="silhouette")


# --------------------------------------------------------------------------- #
# Bayesian NMF, L2/Euclidean, automatic relevance determination
# Direct port of BayesNMF.L2EU from github.com/gwas-partitioning/bnmf-clustering
# (Tan & Fevotte 2013; Udler et al. 2018 PLoS Med).
# --------------------------------------------------------------------------- #
@dataclass
class BNMFRun:
    W: np.ndarray          # variants x K_active
    H: np.ndarray          # K_active x traits(2M)
    evidence: float
    k: int
    iters: int


def bnmf_l2eu(V0: np.ndarray, rng: np.random.Generator, n_iter: int = 10000, a0: float = 10, tol: float = 1e-7,
              K: int = 20, K0: int = 10, phi: float = 1.0, min_iter: int = 100, window: int = 25) -> BNMFRun:
    eps = 1e-50
    V = V0 - V0.min()
    Vmax = V.max()
    N, M = V.shape
    W = rng.random((N, K)) * Vmax
    H = rng.random((K, M)) * Vmax
    V_ap = W @ H + eps
    phi = V.std(ddof=1) ** 2 * phi
    C = (N + M) / 2 + a0 + 1
    b0 = 3.14 * (a0 - 1) * V.mean() / (2 * K0)
    lambda_bound = b0 / C
    lam = (0.5 * (W ** 2).sum(0) + 0.5 * (H ** 2).sum(1) + b0) / C
    lambda_cut = lambda_bound * 1.5
    errors = []
    it = 1
    delta = np.inf
    evid = np.nan
    while it < n_iter:
        it += 1
        inv = 1.0 / lam
        H = H * (W.T @ V) / (W.T @ V_ap + phi * H * inv[:, None] + eps)
        V_ap = W @ H + eps
        W = W * (V @ H.T) / (V_ap @ H.T + phi * W * inv[None, :] + eps)
        V_ap = W @ H + eps
        new = (0.5 * (W ** 2).sum(0) + 0.5 * (H ** 2).sum(1) + b0) / C
        new[~np.isfinite(new)] = 1e-6
        delta = np.max(np.abs(new - lam) / (lam + 1e-10))
        lam = new
        err = ((V - V_ap) ** 2).sum()
        errors.append(err)
        if it > min_iter and delta < tol:
            break
        if it > min_iter and it % 20 == 0 and len(errors) >= window:
            rec = np.array(errors[-window:])
            if abs(np.mean(np.diff(rec)) / np.mean(rec)) < tol:
                break
    like = err / 2
    evid = like + phi * np.sum((0.5 * (W ** 2).sum(0) + 0.5 * (H ** 2).sum(1) + b0) / lam + C * np.log(lam))
    active = np.where(lam >= lambda_cut)[0]
    if len(active) == 0:
        active = np.array([int(np.argmax(lam))])
    return BNMFRun(W=W[:, active], H=H[active, :], evidence=float(evid), k=len(active), iters=it)


def run_bnmf(V: pd.DataFrame, n_reps: int = 100, seed: int = 1, **kw) -> list[BNMFRun]:
    rng = np.random.default_rng(seed)
    return [bnmf_l2eu(V.to_numpy(float), rng, **kw) for _ in range(n_reps)]


def summarise_bnmf(runs: list[BNMFRun], V: pd.DataFrame, k_fixed: int | None = None):
    """Pick K (modal across runs) and the best-evidence run at that K. Return normalised W/H and hard labels."""
    ks = Counter(r.k for r in runs)
    k = k_fixed if k_fixed is not None else max(ks.items(), key=lambda kv: (kv[1], -kv[0]))[0]
    cands = [r for r in runs if r.k == k]
    best = min(cands, key=lambda r: r.evidence)
    W, H = best.W, best.H
    # scale so H rows sum to 1 -> W carries the magnitude; then express each variant as shares across clusters
    s = H.sum(1, keepdims=True) + 1e-12
    H_n, W_n = H / s, W * s.T
    share = W_n / (W_n.sum(1, keepdims=True) + 1e-12)
    labels = share.argmax(1)
    Wdf = pd.DataFrame(W_n, index=V.index, columns=[f"c{i}" for i in range(k)])
    Hdf = pd.DataFrame(H_n, index=Wdf.columns, columns=V.columns)
    sharedf = pd.DataFrame(share, index=V.index, columns=Wdf.columns)
    return dict(k=k, k_counts=dict(sorted(ks.items())), best=best, W=Wdf, H=Hdf, share=sharedf,
                labels=pd.Series(labels, index=V.index, name="bnmf"), runs_at_k=cands)


def h_signed(H: pd.DataFrame) -> pd.DataFrame:
    """Collapse the 2M non-negative columns back to M signed trait loadings (pos - neg)."""
    traits = sorted({c.rsplit("_", 1)[0] for c in H.columns}, key=lambda t: list(H.columns).index(f"{t}_pos"))
    return pd.DataFrame({t: H[f"{t}_pos"] - H[f"{t}_neg"] for t in traits})


def cluster_profiles(z: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Mean scaled z per cluster x trait (used for naming and the heatmap)."""
    return z.groupby(labels.reindex(z.index)).mean()


def name_clusters(profile: pd.DataFrame, top: int = 3, exclude=("sbp",), overall: pd.Series | None = None) -> dict:
    """Name each cluster by the traits where it departs most from the average BP variant (`overall` = mean profile
    over all variants; if None, the unweighted mean of the cluster profiles)."""
    names = {}
    centre = overall if overall is not None else profile.mean()
    for c, row in (profile - centre).iterrows():
        r = row.drop([e for e in exclude if e in row.index])
        top_t = r.abs().sort_values(ascending=False).head(top)
        names[c] = ", ".join(f"{t}{'+' if r[t] > 0 else '-'}" for t in top_t.index)
    return names


def align_labels(ref: np.ndarray, other: np.ndarray, k_ref: int | None = None, k_other: int | None = None) -> np.ndarray:
    """Relabel `other` to best match `ref` (Hungarian on overlap counts)."""
    from scipy.optimize import linear_sum_assignment
    k_ref = k_ref or ref.max() + 1
    k_other = k_other or other.max() + 1
    M = np.zeros((k_other, k_ref))
    np.add.at(M, (other, ref), 1)
    r, c = linear_sum_assignment(-M)
    mapping = {o: rr for o, rr in zip(r, c)}
    nxt = k_ref
    out = np.empty_like(other)
    for i, o in enumerate(other):
        if o not in mapping:
            mapping[o] = nxt
            nxt += 1
        out[i] = mapping[o]
    return out


def ari(a, b) -> float:
    return float(adjusted_rand_score(a, b))
