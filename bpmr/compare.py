"""Between-cluster comparison of MR estimates: heterogeneity, permutation null, pairwise contrasts, and MVMR."""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy import stats

from .mr import ivw


def _ivw_b(bx, by, sy):
    return np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2)


def q_between(b: np.ndarray, se: np.ndarray) -> tuple[float, float]:
    """Cochran's Q across cluster estimates (fixed-effect), df = K - 1."""
    w = 1 / se ** 2
    bbar = np.sum(w * b) / np.sum(w)
    Q = float(np.sum(w * (b - bbar) ** 2))
    return Q, float(stats.chi2.sf(Q, len(b) - 1))


def permutation_test(bx, sx, by, sy, labels: np.ndarray, B: int = 2000, seed: int = 1) -> dict:
    """
    Null: cluster membership is unrelated to the causal estimate.
    Shuffle labels across variants (cluster sizes fixed), recompute Q_between of the IVW estimates.
    """
    ks = np.unique(labels)

    def stat(lab):
        b, se = [], []
        for k in ks:
            m = lab == k
            r = ivw(bx[m], sx[m], by[m], sy[m])
            b.append(r["b"]), se.append(r["se"])
        return q_between(np.array(b), np.array(se))[0]

    obs = stat(labels)
    rng = np.random.default_rng(seed)
    null = np.array([stat(rng.permutation(labels)) for _ in range(B)])
    return dict(Q_between=obs, perm_p=(np.sum(null >= obs) + 1) / (B + 1), null_median=float(np.median(null)))


def pairwise(bx, sx, by, sy, labels: np.ndarray, names: dict, B: int = 2000, seed: int = 1) -> pd.DataFrame:
    """
    For every pair of clusters: difference in log OR (per exposure unit) with a Wald test, plus the ratio of
    log ORs with a percentile bootstrap CI (variants resampled within cluster). The ratio is only reported as
    interpretable when the denominator's 95% CI excludes zero.
    """
    rng = np.random.default_rng(seed)
    ks = np.unique(labels)
    est, boot = {}, {}
    for k in ks:
        m = np.where(labels == k)[0]
        r = ivw(bx[m], sx[m], by[m], sy[m])
        est[k] = r
        bb = np.empty(B)
        for i in range(B):
            s = rng.choice(m, len(m), replace=True)
            bb[i] = _ivw_b(bx[s], by[s], sy[s])
        boot[k] = bb
    rows = []
    for a, b in itertools.combinations(ks, 2):
        if est[a]["b"] < est[b]["b"]:  # orient each pair so cluster_a has the larger effect
            a, b = b, a
        d = est[a]["b"] - est[b]["b"]
        se_d = np.sqrt(est[a]["se"] ** 2 + est[b]["se"] ** 2)
        ratio_boot = boot[a] / boot[b]
        denom_excl0 = (est[b]["b"] - 1.96 * est[b]["se"] > 0) or (est[b]["b"] + 1.96 * est[b]["se"] < 0)
        rows.append(dict(cluster_a=a, cluster_b=b, name_a=names.get(a, a), name_b=names.get(b, b),
                         diff_logOR=d, diff_se=se_d, diff_p=2 * stats.norm.sf(abs(d / se_d)),
                         ratio_logOR=est[a]["b"] / est[b]["b"],
                         ratio_lo=np.nanpercentile(ratio_boot, 2.5), ratio_hi=np.nanpercentile(ratio_boot, 97.5),
                         ratio_interpretable=bool(denom_excl0 and est[b]["b"] > 0)))
    return pd.DataFrame(rows)


def mvmr_ivw(BX: np.ndarray, by: np.ndarray, sy: np.ndarray, names: list[str]) -> pd.DataFrame:
    """Multivariable IVW: weighted regression of by on the exposure columns, no intercept, random-effects SE."""
    w = 1 / sy ** 2
    XtW = BX.T * w
    A = XtW @ BX
    beta = np.linalg.solve(A, XtW @ by)
    resid = by - BX @ beta
    n, p = BX.shape
    phi = max(1.0, np.sum(w * resid ** 2) / (n - p))
    se = np.sqrt(np.diag(np.linalg.inv(A)) * phi)
    return pd.DataFrame(dict(exposure=names, b=beta, se=se, p=2 * stats.t.sf(np.abs(beta / se), df=n - p), nsnp=n))


def conditional_f(BX: np.ndarray, SX: np.ndarray, names: list[str]) -> pd.Series:
    """Sanderson et al. 2019 conditional F-statistic for each exposure (covariance between exposures ignored)."""
    out = {}
    for j in range(BX.shape[1]):
        others = [i for i in range(BX.shape[1]) if i != j]
        y, X = BX[:, j], BX[:, others]
        gamma = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ gamma
        var = SX[:, j] ** 2 + (SX[:, others] ** 2) @ (gamma ** 2)
        Q = np.sum(resid ** 2 / var)
        out[names[j]] = Q / (BX.shape[0] - BX.shape[1] + 1)
    return pd.Series(out, name="conditional_F")


def mr_clusters(gate: pd.Series, sizes: pd.Series, min_snps: int):
    """Clusters that passed the stability gate; if none did, every cluster with enough variants, flagged exploratory."""
    ok = [int(c[1:]) for c, g in gate.items() if g]
    if ok:
        return ok, False
    return [int(c[1:]) for c, n in sizes.items() if n >= min_snps], True
