#!/usr/bin/env python3
"""Step 3: cluster stability (variant bootstrap, sampling-noise bootstrap, leave-one-trait-out, K scan, bNMF)."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr import config  # noqa: E402
from bpmr import stability as st  # noqa: E402
from bpmr.cluster import align_labels, ari, bnmf_l2eu, hclust_input, kmeans, ward  # noqa: E402
from bpmr.matrix import ZMatrix, profile_matrix  # noqa: E402


def main():
    s = config.parse(__doc__)
    prep, cl, out = s.outp / "01_prepare", s.outp / "02_cluster", s.step("03_stability")
    zfull = pd.read_csv(prep / "z_scaled.csv", index_col=0)
    z = pd.read_csv(prep / "z_cluster.csv", index_col=0)
    zfull = zfull.reindex(z.index)
    noise = pd.read_csv(prep / "noise_sd.csv", index_col=0).reindex(index=zfull.index, columns=zfull.columns)
    labels = pd.read_csv(cl / "labels.csv", index_col=0).reindex(z.index)
    meta = json.loads((cl / "clusters.json").read_text())
    k = meta["k"]
    log = []

    x = hclust_input(z, s.hclust_clip_sd)
    methods = {
        "kmeans": lambda xx, kk: kmeans(xx, kk, n_init=s.kmeans_n_init, seed=s.seed % 2**31),
        "hclust": lambda xx, kk: ward(xx, kk),
    }
    summaries = {}
    for part, fk in methods.items():
        lab_main = labels[part].to_numpy()
        fn = lambda xx, fk=fk: fk(xx, k)  # noqa: E731

        def fn_from_z(zz, fk=fk):  # perturbed full z -> same profile + clipping -> same method
            zc = profile_matrix(pd.DataFrame(zz, columns=zfull.columns), s.profile)
            return fk(hclust_input(zc, s.hclust_clip_sd), k)

        J = st.variant_bootstrap(x, lab_main, fn, B=s.boot_variant, seed=s.seed)
        J.to_csv(out / f"{part}_bootstrap_jaccard.csv", index=False)
        keep, co, aris = st.noise_bootstrap(zfull.to_numpy(), noise.fillna(noise.mean()).to_numpy(), lab_main, fn_from_z,
                                            B=s.boot_noise, seed=s.seed + 1)
        np.save(out / f"{part}_noise_coassign.npy", co.astype(np.float32))
        summ = st.summarise(lab_main, J, keep, co)
        summ["name"] = summ["cluster"].map(lambda c: meta["names"][part][c[1:]])
        summ["goes_to_mr"] = (summ["jaccard_mean"] >= s.jaccard_min) & (summ["n_variants"] >= s.min_cluster_snps)
        summ.to_csv(out / f"{part}_stability_summary.csv", index=False)
        summaries[part] = summ
        pd.DataFrame(dict(rsid=z.index, keep_prob=keep, cluster=lab_main)).to_csv(out / f"{part}_variant_keep_prob.csv", index=False)
        log.append(f"{part} noise-bootstrap ARI vs main: median {np.median(aris):.3f} "
                   f"(IQR {np.percentile(aris, 25):.3f}-{np.percentile(aris, 75):.3f})")
        loto = st.leave_one_trait_out(z, lab_main, lambda zz: hclust_input(zz, s.hclust_clip_sd), fn)
        loto.to_csv(out / f"{part}_leave_one_trait_out_ARI.csv")
        st.k_sensitivity(x, lab_main, fk, k).to_csv(out / f"{part}_k_sensitivity_ARI.csv")
        kscan = st.stability_by_k(x, fk, ks=range(2, min(11, len(z) // 10)), B=s.boot_k_scan, seed=s.seed + 2)
        kscan.to_csv(out / f"{part}_stability_by_k.csv", index=False)
        log.append(f"{part} stability by K (mean Jaccard): " + ", ".join(f"{r.k}:{r.mean_jaccard:.2f}" for r in kscan.itertuples()))
        log.append(f"{part} leave-one-trait-out ARI (lowest 5): " + ", ".join(f"{t_}:{v:.2f}" for t_, v in loto.head(5).items()))

    # ---------- bNMF ----------
    co_runs = np.load(cl / "bnmf_run_coassign.npy")
    bl = labels["bnmf"].to_numpy()
    V = ZMatrix(z=z, z_raw=z, n=z, noise_sd=z).nonneg()
    rng = np.random.default_rng(s.seed + 3)
    bn_aris, bn_keep = [], np.zeros(len(z))
    for b in range(s.bnmf_noise_runs):
        zb = profile_matrix(zfull + rng.normal(0, 1, zfull.shape) * noise.fillna(noise.mean()).to_numpy(), s.profile)
        Vb = ZMatrix(z=zb, z_raw=zb, n=zb, noise_sd=zb).nonneg()
        runs = [bnmf_l2eu(Vb.to_numpy(float), rng, K=s.bnmf_k_init, K0=s.bnmf_k0, phi=s.bnmf_phi) for _ in range(3)]
        best = min(runs, key=lambda r: r.evidence)
        lab = np.argmax(best.W * best.H.sum(1)[None, :], axis=1)
        lab = align_labels(bl, lab, k_ref=int(bl.max()) + 1)
        bn_aris.append(ari(bl, lab))
        bn_keep += lab == bl
    bn_rows = []
    for c in sorted(np.unique(bl)):
        m = bl == c
        w = co_runs[np.ix_(m, m)]
        bn_rows.append(dict(cluster=f"c{c}", n_variants=int(m.sum()), name=meta["names"]["bnmf"].get(str(c), ""),
                            run_coassign_within=w[np.triu_indices(m.sum(), 1)].mean() if m.sum() > 1 else np.nan,
                            noise_keep_prob=(bn_keep[m] / max(1, s.bnmf_noise_runs)).mean()))
    bsum = pd.DataFrame(bn_rows)
    # bNMF has no variant bootstrap (too slow); its MR gate uses the noise-bootstrap keep probability instead
    bsum["goes_to_mr"] = (bsum["noise_keep_prob"] >= s.jaccard_min) & (bsum["n_variants"] >= s.min_cluster_snps)
    bsum.to_csv(out / "bnmf_stability_summary.csv", index=False)
    log.append(f"bNMF noise-bootstrap ARI vs main ({s.bnmf_noise_runs} reps): median {np.median(bn_aris):.3f}")

    for part, summ in summaries.items():
        log.append(f"{part} per-cluster stability (K = {k}):")
        for r in summ.itertuples():
            log.append(f"   {r.cluster} n={r.n_variants:4d} Jaccard={r.jaccard_mean:.2f} keep_prob={r.noise_keep_prob:.2f} "
                       f"[{r.stability}] MR={'yes' if r.goes_to_mr else 'flagged'} :: {r.name}")
    (out / "log.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
