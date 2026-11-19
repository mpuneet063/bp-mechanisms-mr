#!/usr/bin/env python3
"""Step 2: k-means with K chosen by bootstrap stability (primary), Ward at the same K, Bayesian NMF with its own K;
concordance and cluster names."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr import stability as st  # noqa: E402
from bpmr.cluster import (align_labels, ari, choose_k, cluster_profiles, h_signed, hclust_input, kmeans,  # noqa: E402
                          name_clusters, run_bnmf, silhouette_by_k, summarise_bnmf, ward)
from bpmr.matrix import ZMatrix  # noqa: E402


def main():
    s = config.parse(__doc__)
    prep, out = s.outp / "01_prepare", s.step("02_cluster")
    z = pd.read_csv(prep / "z_cluster.csv", index_col=0)
    V = ZMatrix(z=z, z_raw=z, n=z, noise_sd=z).nonneg()
    log = [f"matrix: {z.shape[0]} variants x {z.shape[1]} traits"]

    # --- bNMF ---
    scan = []
    for phi in s.bnmf_phi_scan:
        rs = run_bnmf(V, n_reps=s.bnmf_phi_scan_runs, seed=s.seed + 7, K=s.bnmf_k_init, K0=s.bnmf_k0, phi=phi)
        for r in rs:
            scan.append(dict(phi=phi, K=r.k))
    scan = pd.DataFrame(scan).groupby(["phi", "K"]).size().rename("runs").reset_index()
    scan.to_csv(out / "bnmf_phi_scan.csv", index=False)
    log.append("bNMF phi scan (phi: K counts): " + "; ".join(
        f"{p}: " + ", ".join(f"K{r.K}x{r.runs}" for r in g.itertuples()) for p, g in scan.groupby("phi")))
    runs = run_bnmf(V, n_reps=s.bnmf_runs, seed=s.seed, K=s.bnmf_k_init, K0=s.bnmf_k0, phi=s.bnmf_phi)
    bn = summarise_bnmf(runs, V)
    kb = bn["k"]
    log.append(f"bNMF K across {s.bnmf_runs} runs: {bn['k_counts']} -> K_bNMF = {kb}")
    pd.Series(bn["k_counts"], name="runs").rename_axis("K").to_csv(out / "bnmf_K_counts.csv")
    # bNMF co-assignment across the runs that reached K (Udler's own stability check)
    labs = np.array([np.argmax(r.W / (r.W.sum(1, keepdims=True) + 1e-12) * r.H.sum(1)[None, :], axis=1) for r in bn["runs_at_k"]])
    co = (labs[:, :, None] == labs[:, None, :]).mean(0)
    np.save(out / "bnmf_run_coassign.npy", co.astype(np.float32))

    # --- K by stability (k-means), then k-means and Ward at that K ---
    x = hclust_input(z, s.hclust_clip_sd)
    km = lambda xx, kk: kmeans(xx, kk, n_init=s.kmeans_n_init, seed=s.seed % 2**31)  # noqa: E731
    kscan = st.stability_by_k(x, km, ks=s.k_scan, B=s.boot_k_select, seed=s.seed + 11)
    kscan.to_csv(out / "kmeans_k_selection.csv", index=False)
    k, why = (s.k_override, "forced with --k") if s.k_override else choose_k(kscan, s.jaccard_min, s.jaccard_stable)
    log.append("k-means stability by K (mean / weakest-cluster Jaccard): "
               + ", ".join(f"{r.k}:{r.mean_jaccard:.2f}/{r.min_cluster_jaccard:.2f}" for r in kscan.itertuples()))
    log.append(f"K = {k} ({why})")
    kl = km(x, k)
    hl = align_labels(kl, ward(x, k), k_ref=k)
    sil = silhouette_by_k(x)
    sil.rename_axis("k").to_csv(out / "hclust_silhouette.csv")
    bl = align_labels(kl, bn["labels"].to_numpy(), k_ref=k)
    log.append(f"Ward silhouette by k: " + ", ".join(f"{i}:{v:.3f}" for i, v in sil.items()))
    aris = {"kmeans_vs_hclust": ari(kl, hl), "kmeans_vs_bnmf": ari(kl, bl), "hclust_vs_bnmf": ari(hl, bl)}
    log.append("concordance (ARI): " + ", ".join(f"{a} {v:.3f}" for a, v in aris.items()))

    labels = pd.DataFrame(dict(kmeans=kl, hclust=hl, bnmf=bl), index=z.index)
    # rename bNMF share columns so c<i> matches the aligned bNMF label i
    raw_lab = bn["labels"].to_numpy()
    mapping = {int(r): int(a) for r, a in zip(raw_lab, bl)}
    nxt = int(max(bl.max() + 1, k))
    for c in range(bn["share"].shape[1]):
        if c not in mapping:
            mapping[c], nxt = nxt, nxt + 1
    share = bn["share"].copy()
    share.columns = [f"bnmf_share_c{mapping[int(c[1:])]}" for c in share.columns]
    labels = labels.join(share[sorted(share.columns)])
    labels.index.name = "rsid"
    labels.to_csv(out / "labels.csv")
    pd.crosstab(labels["kmeans"], labels["hclust"]).to_csv(out / "kmeans_vs_hclust_crosstab.csv")
    pd.crosstab(labels["kmeans"], labels["bnmf"]).to_csv(out / "kmeans_vs_bnmf_crosstab.csv")

    names = {}
    for part in config.PARTITIONS:
        prof = cluster_profiles(z, labels[part])
        prof.to_csv(out / f"profile_{part}.csv")
        names[part] = {int(c): n for c, n in name_clusters(prof, overall=z.mean()).items()}
    h_signed(bn["H"]).to_csv(out / "bnmf_H_signed.csv")
    sizes = {p: labels[p].value_counts().sort_index().to_dict() for p in config.PARTITIONS}
    meta = dict(k=int(k), k_reason=why, k_bnmf=int(kb), k_counts={int(a): int(b) for a, b in bn["k_counts"].items()}, ari=aris,
                names=names, sizes={p: {int(a): int(b) for a, b in d.items()} for p, d in sizes.items()},
                primary=s.primary)
    (out / "clusters.json").write_text(json.dumps(meta, indent=2))
    for part in config.PARTITIONS:
        log.append(f"{part} clusters:")
        for c in sorted(names[part]):
            log.append(f"   c{c} (n={sizes[part][c]}): {names[part][c]}")
    (out / "log.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
