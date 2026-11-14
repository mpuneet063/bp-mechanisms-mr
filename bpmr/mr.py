"""Two-sample Mendelian randomisation estimators (summary data). All inputs are 1-D numpy arrays aligned to the
exposure-raising allele: bx, sx (exposure), by, sy (outcome)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def ivw(bx, sx, by, sy):
    """Inverse-variance weighted, multiplicative random effects (SE inflated when Q/(n-1) > 1)."""
    w = bx ** 2 / sy ** 2
    b = np.sum(bx * by / sy ** 2) / np.sum(w)
    se_fe = 1 / np.sqrt(np.sum(w))
    n = len(bx)
    resid = (by - b * bx) / sy
    Q = float(np.sum(resid ** 2))
    phi = Q / (n - 1) if n > 1 else 1.0
    se = se_fe * np.sqrt(max(1.0, phi))
    p = 2 * stats.t.sf(abs(b / se), df=max(n - 1, 1)) if n > 1 else np.nan
    return dict(method="IVW", b=b, se=se, p=p, nsnp=n, Q=Q, Q_df=n - 1,
                Q_p=stats.chi2.sf(Q, n - 1) if n > 1 else np.nan,
                I2=max(0.0, (Q - (n - 1)) / Q) if Q > 0 else 0.0)


def egger(bx, sx, by, sy):
    """MR-Egger: weighted regression of by on bx with intercept (bx oriented positive)."""
    n = len(bx)
    if n < 3:
        return dict(method="MR-Egger", b=np.nan, se=np.nan, p=np.nan, nsnp=n, intercept=np.nan, int_se=np.nan, int_p=np.nan)
    s = np.sign(bx)
    bx, by = bx * s, by * s
    w = 1 / sy ** 2
    X = np.column_stack([np.ones(n), bx])
    XtW = X.T * w
    beta = np.linalg.solve(XtW @ X, XtW @ by)
    resid = by - X @ beta
    sigma = np.sqrt(np.sum(w * resid ** 2) / (n - 2))
    cov = np.linalg.inv(XtW @ X) * (max(1.0, sigma) ** 2)  # residual SE floored at 1 (as TwoSampleMR)
    se = np.sqrt(np.diag(cov))
    p = 2 * stats.t.sf(np.abs(beta / se), df=n - 2)
    return dict(method="MR-Egger", b=beta[1], se=se[1], p=p[1], nsnp=n, intercept=beta[0], int_se=se[0], int_p=p[0])


def _wm(b, w):
    o = np.argsort(b)
    b, w = b[o], w[o] / w.sum()
    cum = np.cumsum(w) - 0.5 * w
    i = np.searchsorted(cum, 0.5)
    if i == 0:
        return b[0]
    if i >= len(b):
        return b[-1]
    return b[i - 1] + (b[i] - b[i - 1]) * (0.5 - cum[i - 1]) / (cum[i] - cum[i - 1])


def weighted_median(bx, sx, by, sy, nboot: int = 1000, seed: int = 1):
    n = len(bx)
    if n < 3:
        return dict(method="Weighted median", b=np.nan, se=np.nan, p=np.nan, nsnp=n)
    br = by / bx
    se_r = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)  # second-order weights
    b = _wm(br, 1 / se_r ** 2)
    rng = np.random.default_rng(seed)
    boots = np.empty(nboot)
    for i in range(nboot):
        bxb, byb = rng.normal(bx, sx), rng.normal(by, sy)
        boots[i] = _wm(byb / bxb, 1 / se_r ** 2)
    se = boots.std(ddof=1)
    return dict(method="Weighted median", b=b, se=se, p=2 * stats.norm.sf(abs(b / se)), nsnp=n)


def presso(bx, sx, by, sy, nsim: int = 2000, alpha: float = 0.05, seed: int = 1):
    """MR-PRESSO global test and outlier-corrected IVW (Verbanck et al. 2018), compact re-implementation."""
    n = len(bx)
    if n < 4:
        return dict(method="MR-PRESSO (outlier-corrected)", b=np.nan, se=np.nan, p=np.nan, nsnp=n,
                    global_p=np.nan, n_outliers=0, outliers=[])
    w = 1 / sy ** 2

    def loo_b(bxv, byv):
        tot_xy, tot_xx = np.sum(w * bxv * byv), np.sum(w * bxv ** 2)
        return (tot_xy - w * bxv * byv) / (tot_xx - w * bxv ** 2)

    b_loo = loo_b(bx, by)
    rss_obs_i = (by - b_loo * bx) ** 2 / sy ** 2
    rss_obs = rss_obs_i.sum()
    # the Bonferroni outlier test needs at least n/alpha simulations to be able to reject (as in MR-PRESSO)
    nsim = int(max(nsim, np.ceil(2 * n / alpha)))
    rng = np.random.default_rng(seed)
    sim_i = np.empty((nsim, n))
    for start in range(0, nsim, 2000):  # vectorised in blocks
        m = min(2000, nsim - start)
        bx_s = rng.normal(bx, sx, size=(m, n))
        by_s = rng.normal(b_loo * bx, sy, size=(m, n))
        tot_xy = np.sum(w * bx_s * by_s, axis=1, keepdims=True)
        tot_xx = np.sum(w * bx_s ** 2, axis=1, keepdims=True)
        bl = (tot_xy - w * bx_s * by_s) / (tot_xx - w * bx_s ** 2)
        sim_i[start:start + m] = (by_s - bl * bx_s) ** 2 / sy ** 2
    rss_sim = sim_i.sum(1)
    global_p = (np.sum(rss_sim >= rss_obs) + 1) / (nsim + 1)
    outliers = np.array([], dtype=int)
    if global_p < alpha:
        p_i = (np.sum(sim_i >= rss_obs_i[None, :], axis=0) + 1) / (nsim + 1)
        outliers = np.where(p_i * n < alpha)[0]
    keep = np.setdiff1d(np.arange(n), outliers)
    r = ivw(bx[keep], sx[keep], by[keep], sy[keep])
    r.update(method="MR-PRESSO (outlier-corrected)", global_p=global_p, n_outliers=len(outliers), outliers=outliers.tolist())
    return r


def r2_continuous(z, n):
    return z ** 2 / (z ** 2 + n - 2)


def r2_binary(b, eaf, ncase, ncontrol):
    """Variance explained in a binary outcome from a log OR on the logistic latent scale
    (r2 = var_x b^2 / (var_x b^2 + pi^2/3)), without the prevalence adjustment TwoSampleMR applies.
    An approximation: the Steiger filter is used only as a sensitivity analysis. ncase/ncontrol kept for API symmetry."""
    var_x = 2 * eaf * (1 - eaf)
    vl = var_x * b ** 2
    return vl / (vl + np.pi ** 2 / 3)


def steiger_keep(bx, sx, nx, by, sy, eaf, ncase, ncontrol, binary: bool):
    rx = r2_continuous(bx / sx, nx)
    ry = r2_binary(by, eaf, ncase, ncontrol) if binary else r2_continuous(by / sy, ncase + ncontrol)
    return rx > ry


def run_all(bx, sx, by, sy, nboot=1000, nsim=2000, seed=1) -> pd.DataFrame:
    rows = [ivw(bx, sx, by, sy), egger(bx, sx, by, sy), weighted_median(bx, sx, by, sy, nboot, seed),
            presso(bx, sx, by, sy, nsim, seed=seed)]
    df = pd.DataFrame(rows)
    df["F_mean"] = np.mean((bx / sx) ** 2)
    return df
