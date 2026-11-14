"""
Synthetic data in exactly the format written by scripts/00_fetch_opengwas.py, with known truth.

Used to check the whole pipeline recovers (a) the planted mechanism clusters and (b) the planted
cluster-specific causal effects, and that MVMR removes planted BMI pleiotropy. Never used for results.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

TRAITS = ["sbp", "dbp", "pp", "bmi", "whradjbmi", "height", "ldl", "hdl", "tg", "t2d", "fg", "fi", "hba1c",
          "egfr", "urate", "uacr", "na_urine", "k_urine", "hr", "crp", "alt", "alcohol", "smoking"]
OUTCOMES = {  # key: (ncase, ncontrol)
    "stroke_any": (73652, 1234808), "stroke_is": (62100, 1234808), "stroke_las": (6399, 1234808),
    "stroke_ces": (10804, 1234808), "stroke_svs": (6811, 1234808), "ich": (1545, 1481),
    "cad": (181522, 984168),
}
REPL = {"fg_stroke": (10551, 202223), "fg_is": (8046, 202223), "fg_chd": (21012, 197780)}
BINARY_TRAITS = {"t2d": (74124, 824006)}

# planted mechanism profiles (expected z per unit of a variant's "strength"; sbp is the exposure itself)
PROFILES = {
    "vascular": dict(pp=1.0, dbp=-0.2, hr=-0.3),
    "adiposity": dict(bmi=1.0, whradjbmi=0.5, t2d=0.6, tg=0.5, hdl=-0.5, crp=0.5, fi=0.5),
    "kidney": dict(egfr=-0.9, urate=0.7, uacr=0.6, na_urine=-0.4),
    "lipid": dict(ldl=0.9, tg=0.4, hdl=-0.3, alt=0.3),
    "cardiac": dict(dbp=0.9, hr=0.8, pp=-0.6),
    "bp_only": dict(),
}
# planted causal effects: log OR per 10 mmHg SBP acting through each mechanism
THETA = {
    "vascular":  dict(stroke_any=0.40, stroke_is=0.38, stroke_las=0.30, stroke_ces=0.10, stroke_svs=0.70, ich=0.60, cad=0.25),
    "adiposity": dict(stroke_any=0.20, stroke_is=0.20, stroke_las=0.30, stroke_ces=0.20, stroke_svs=0.15, ich=0.10, cad=0.25),
    "kidney":    dict(stroke_any=0.25, stroke_is=0.25, stroke_las=0.25, stroke_ces=0.15, stroke_svs=0.35, ich=0.20, cad=0.45),
    "lipid":     dict(stroke_any=0.25, stroke_is=0.28, stroke_las=0.45, stroke_ces=0.10, stroke_svs=0.20, ich=0.10, cad=0.45),
    "cardiac":   dict(stroke_any=0.30, stroke_is=0.30, stroke_las=0.25, stroke_ces=0.45, stroke_svs=0.25, ich=0.30, cad=0.20),
    "bp_only":   dict(stroke_any=0.30, stroke_is=0.30, stroke_las=0.30, stroke_ces=0.20, stroke_svs=0.35, ich=0.35, cad=0.30),
}
GAMMA_BMI_CAD = 0.12   # planted pleiotropy: log OR for CAD per SD BMI (acts only through the BMI association)
SIZES = dict(vascular=120, adiposity=90, kidney=80, lipid=70, cardiac=90, bp_only=250)


def _alleles(rng, n, pal_frac=0.03):
    pairs = [("A", "C"), ("A", "G"), ("C", "T"), ("G", "T")]
    idx = rng.integers(0, 4, n)
    ea = np.array([pairs[i][0] for i in idx], dtype=object)
    nea = np.array([pairs[i][1] for i in idx], dtype=object)
    pal = rng.random(n) < pal_frac
    ea[pal], nea[pal] = "A", "T"
    return ea, nea


def make(out_dir: str | Path, seed: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ------------------- variants ------------------- #
    rows = []
    for mech, n in SIZES.items():
        for _ in range(n):
            rows.append(dict(mech=mech))
    for kind, n in [("mvmr_bmi", 250), ("mvmr_ldl", 150), ("mvmr_t2d", 120)]:
        for _ in range(n):
            rows.append(dict(mech=kind))
    v = pd.DataFrame(rows)
    v["rsid"] = [f"rs{1000000 + i}" for i in range(len(v))]
    v["ea"], v["nea"] = _alleles(rng, len(v))
    v["eaf"] = rng.uniform(0.05, 0.95, len(v))
    bp = ~v["mech"].str.startswith("mvmr")
    v["strength"] = rng.gamma(4, 2.0, len(v))                      # expected |z| scale of non-BP traits
    v["bx_true"] = np.where(bp, rng.uniform(0.18, 0.7, len(v)), rng.normal(0, 0.05, len(v)))  # mmHg

    n_trait = {t: int(rng.integers(150_000, 750_000)) for t in TRAITS}
    n_trait["sbp"] = 757_601
    se_sbp = 1 / np.sqrt(n_trait["sbp"]) * 19.0                    # mmHg scale (SD ~ 19 mmHg)

    true = {}
    for t in TRAITS:
        if t == "sbp":
            continue
        mu = np.zeros(len(v))
        for mech, prof in PROFILES.items():
            m = (v["mech"] == mech).to_numpy()
            mu[m] = prof.get(t, 0.0) * v.loc[m, "strength"]
        mu += rng.normal(0, 0.5, len(v)) * bp                       # background pleiotropy
        true[t] = mu
    # mediator instruments: strong on their own trait
    for kind, t in [("mvmr_bmi", "bmi"), ("mvmr_ldl", "ldl"), ("mvmr_t2d", "t2d")]:
        m = (v["mech"] == kind).to_numpy()
        true[t][m] = rng.choice([-1, 1], m.sum()) * rng.uniform(6, 15, m.sum())

    # ------------------- associations ------------------- #
    def emit(key, beta, se, n):
        flip = rng.random(len(v)) < 0.5            # report half on the other allele
        ea = np.where(flip, v["nea"], v["ea"])
        nea = np.where(flip, v["ea"], v["nea"])
        eaf = np.where(flip, 1 - v["eaf"], v["eaf"])
        b = np.where(flip, -beta, beta)
        miss = rng.random(len(v)) < (0.0 if key == "sbp" else 0.04)
        df = pd.DataFrame(dict(id=f"sim-{key}", trait=key, rsid=v["rsid"], chr=1, position=np.arange(len(v)),
                               ea=ea, nea=nea, eaf=eaf, beta=b, se=se, p=0.0, n=n, proxy=False))
        from scipy.stats import norm
        df["p"] = 2 * norm.sf(np.abs(df["beta"] / df["se"]))
        return df[~miss]

    parts = []
    bx_obs = v["bx_true"] + rng.normal(0, se_sbp, len(v))
    parts.append(emit("sbp", bx_obs.to_numpy(), np.full(len(v), se_sbp), n_trait["sbp"]))
    v["bx_obs"] = bx_obs
    trait_beta_true = {}
    for t in TRAITS:
        if t == "sbp":
            continue
        se = np.full(len(v), 1 / np.sqrt(n_trait[t]))
        bt = true[t] * se
        if t in ("dbp", "pp"):  # BP components also carry the SBP signal
            bt = bt + v["bx_true"].to_numpy() * (0.6 if t == "dbp" else 0.5) / 19.0
        trait_beta_true[t] = bt
        parts.append(emit(t, bt + rng.normal(0, 1, len(v)) * se, se, n_trait[t]))

    def outcome_se(nc, nn):
        f = v["eaf"].to_numpy()
        return np.sqrt((1 / nc + 1 / nn) / (2 * f * (1 - f)))

    bmi_sd_beta = trait_beta_true["bmi"]  # per SD BMI
    for key, (nc, nn) in {**OUTCOMES, **REPL}.items():
        base = key.replace("fg_", "")
        base = {"stroke": "stroke_any", "is": "stroke_is", "chd": "cad"}.get(base, base)
        theta = v["mech"].map(lambda m: THETA.get(m, {}).get(base, 0.0)).to_numpy() / 10.0
        by = theta * v["bx_true"].to_numpy()
        if base == "cad":
            by = by + GAMMA_BMI_CAD * bmi_sd_beta
        se = outcome_se(nc, nn)
        parts.append(emit(key, by + rng.normal(0, 1, len(v)) * se, se, nc + nn))

    assoc = pd.concat(parts, ignore_index=True)
    assoc.to_csv(out / "associations.csv.gz", index=False, compression="gzip")

    # ------------------- manifest / instruments ------------------- #
    man = [dict(key="sbp", role="exposure", id="sim-sbp", trait="systolic blood pressure", unit="mmHg",
                sample_size=n_trait["sbp"], ncase=None, ncontrol=None, how="simulated")]
    for t in TRAITS[1:]:
        nc = BINARY_TRAITS.get(t, (None, None))
        man.append(dict(key=t, role="cluster_trait", id=f"sim-{t}", trait=t, unit="SD", sample_size=n_trait[t],
                        ncase=nc[0], ncontrol=nc[1], how="simulated"))
    for k, (nc, nn) in OUTCOMES.items():
        man.append(dict(key=k, role="outcome", id=f"sim-{k}", trait=k, unit="log odds", sample_size=nc + nn,
                        ncase=nc, ncontrol=nn, how="simulated"))
    for k, (nc, nn) in REPL.items():
        man.append(dict(key=k, role="replication", id=f"sim-{k}", trait=k, unit="log odds", sample_size=nc + nn,
                        ncase=nc, ncontrol=nn, how="simulated"))
    pd.DataFrame(man).to_csv(out / "manifest.csv", index=False)

    src = np.where(bp, "sbp_tophits_r2_0.05", v["mech"])
    inst = pd.DataFrame(dict(rsid=v["rsid"], source=src))
    inst["clump_r2_0.05"] = bp
    inst["clump_r2_0.001"] = bp & (rng.random(len(v)) > 0.12)   # some BP variants "in LD" -> not MR instruments
    inst.to_csv(out / "instruments.csv", index=False)
    truth = v[["rsid", "mech", "bx_true", "strength"]]
    truth.to_csv(out / "truth.csv", index=False)
    return truth


if __name__ == "__main__":
    import sys
    make(sys.argv[1] if len(sys.argv) > 1 else "data/sim_raw")
