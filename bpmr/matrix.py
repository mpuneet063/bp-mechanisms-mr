"""Variant x trait z-score matrix: filtering, trait pruning, sample-size scaling (follows Udler/Smith bNMF prep)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import norm

from .data import wide


@dataclass
class ZMatrix:
    z: pd.DataFrame                 # scaled z, NaN filled with 0  (variants x traits)
    z_raw: pd.DataFrame             # unscaled z with NaN
    n: pd.DataFrame                 # sample sizes
    noise_sd: pd.DataFrame          # sd of sampling noise in scaled units (for the noise bootstrap)
    log: list = field(default_factory=list)

    def nonneg(self) -> pd.DataFrame:
        pos = self.z.clip(lower=0).add_suffix("_pos")
        neg = (-self.z).clip(lower=0).add_suffix("_neg")
        return pd.concat([pos, neg], axis=1)


def build_matrix(h: pd.DataFrame, rsids: list[str], traits: list[str], max_missing_variant: float = 0.2,
                 max_missing_trait: float = 0.5, corr_cutoff: float = 0.85, keep_traits=("sbp",),
                 prune_by_signal: bool = True) -> ZMatrix:
    log = []
    z = wide(h, "z", traits, rsids)
    n = wide(h, "n", traits, rsids)
    log.append(f"start: {z.shape[0]} variants x {z.shape[1]} traits")

    tmiss = z.isna().mean()
    drop_t = tmiss[tmiss > max_missing_trait].index.tolist()
    if drop_t:
        log.append(f"traits dropped for >{max_missing_trait:.0%} missing: {drop_t}")
    z, n = z.drop(columns=drop_t), n.drop(columns=drop_t)

    vmiss = z.isna().mean(axis=1)
    z, n = z[vmiss <= max_missing_variant], n[vmiss <= max_missing_variant]
    log.append(f"variants kept with <= {max_missing_variant:.0%} missing traits: {z.shape[0]}")

    if prune_by_signal:  # a trait no variant is associated with (Bonferroni) adds only noise
        minp = pd.DataFrame(2 * norm.sf(z.abs().to_numpy()), index=z.index, columns=z.columns).min()
        weak = minp[minp >= 0.05 / len(z)].index.difference(list(keep_traits)).tolist()
        if weak:
            log.append(f"traits dropped, no variant reaches p < 0.05/{len(z)}: {weak}")
        z, n = z.drop(columns=weak), n.drop(columns=weak)

    # correlation pruning, strongest-signal trait first (bNMF prep)
    cor = z.corr()
    order = z.abs().max().sort_values(ascending=False).index.tolist()
    order = [t for t in keep_traits if t in order] + [t for t in order if t not in keep_traits]
    kept, removed = [], {}
    remaining = order[:]
    while remaining:
        t = remaining.pop(0)
        kept.append(t)
        close = [u for u in remaining if abs(cor.loc[t, u]) >= corr_cutoff and u not in keep_traits]
        for u in close:
            removed[u] = t
        remaining = [u for u in remaining if u not in close]
    if removed:
        log.append("traits dropped for |r| >= %.2f with a kept trait: %s" % (corr_cutoff, removed))
    z, n = z[kept], n[kept]

    # sample-size adjustment: z / sqrt(N) * mean(sqrt(median N))
    n_filled = n.fillna(n.median())
    scale = np.sqrt(n_filled.median()).mean()
    z_scaled = z / np.sqrt(n_filled) * scale
    noise_sd = (1.0 / np.sqrt(n_filled)) * scale
    log.append(f"final: {z.shape[0]} variants x {z.shape[1]} traits; {int(z.isna().sum().sum())} missing cells set to 0")
    return ZMatrix(z=z_scaled.fillna(0.0), z_raw=z, n=n, noise_sd=noise_sd, log=log)


def profile_matrix(z: pd.DataFrame, mode: str = "per_sbp") -> pd.DataFrame:
    """
    Clustering input.
      per_sbp    : each trait's scaled z divided by the variant's SBP z (times the median SBP z), SBP column removed.
                   Every variant was selected for raising SBP, so this is "how much does the variant move trait X
                   per unit of SBP it raises" - the mechanism profile, without the effect-size magnitude.
      raw_no_sbp : scaled z with the SBP column removed
      raw        : scaled z as is (closest to Pascat et al. 2026)
    """
    if mode == "raw":
        return z
    nos = z.drop(columns=["sbp"])
    if mode == "raw_no_sbp":
        return nos
    if mode == "per_sbp":
        return nos.div(z["sbp"], axis=0) * z["sbp"].median()
    raise ValueError(mode)
