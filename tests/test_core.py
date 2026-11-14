"""Unit tests for harmonisation, MR estimators, bNMF and the comparison statistics. Run: pytest -q"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpmr.cluster import align_labels, bnmf_l2eu, ward  # noqa: E402
from bpmr.compare import mvmr_ivw, permutation_test, q_between  # noqa: E402
from bpmr.data import RawData, harmonise  # noqa: E402
from bpmr.matrix import profile_matrix  # noqa: E402
from bpmr.mr import egger, ivw, presso, weighted_median  # noqa: E402


def _sim_mr(theta=0.3, n=200, pleio=0.0, seed=0):
    rng = np.random.default_rng(seed)
    bx_t = rng.uniform(0.05, 0.3, n)
    sx = np.full(n, 0.01)
    sy = np.full(n, 0.02)
    bx = bx_t + rng.normal(0, sx)
    by = theta * bx_t + pleio + rng.normal(0, sy)
    return bx, sx, by, sy


# ---------------- harmonisation ---------------- #
def test_harmonise_flips_and_strands():
    man = pd.DataFrame(dict(key=["sbp", "bmi"], role=["exposure", "cluster_trait"], id=["e", "t"]))
    assoc = pd.DataFrame([
        # rs1: SBP effect allele A raises SBP; trait reported on swapped alleles -> sign flips
        dict(id="e", rsid="rs1", ea="A", nea="G", eaf=0.3, beta=0.5, se=0.1, p=1e-6, n=1),
        dict(id="t", rsid="rs1", ea="G", nea="A", eaf=0.7, beta=0.2, se=0.1, p=0.05, n=1),
        # rs2: SBP effect allele lowers SBP -> reference becomes the other allele; trait on opposite strand
        dict(id="e", rsid="rs2", ea="C", nea="T", eaf=0.2, beta=-0.4, se=0.1, p=1e-6, n=1),
        dict(id="t", rsid="rs2", ea="G", nea="A", eaf=0.2, beta=0.3, se=0.1, p=0.05, n=1),
        # rs3: ambiguous palindrome (A/T, eaf ~0.5) -> dropped
        dict(id="e", rsid="rs3", ea="A", nea="T", eaf=0.5, beta=0.4, se=0.1, p=1e-6, n=1),
        dict(id="t", rsid="rs3", ea="A", nea="T", eaf=0.5, beta=0.3, se=0.1, p=0.05, n=1),
    ])
    h = harmonise(RawData(man, assoc, pd.DataFrame()))
    get = lambda r, k: h[(h.rsid == r) & (h.key == k)]["beta"]  # noqa: E731
    assert get("rs1", "sbp").iloc[0] == pytest.approx(0.5)
    assert get("rs1", "bmi").iloc[0] == pytest.approx(-0.2)
    assert get("rs2", "sbp").iloc[0] == pytest.approx(0.4)      # re-oriented to the SBP-raising allele (T)
    assert get("rs2", "bmi").iloc[0] == pytest.approx(-0.3)     # G/A on the other strand = C/T; C is SBP-lowering
    assert get("rs3", "bmi").empty


# ---------------- MR estimators ---------------- #
def test_ivw_recovers_theta():
    bx, sx, by, sy = _sim_mr(0.3)
    r = ivw(bx, sx, by, sy)
    assert abs(r["b"] - 0.3) < 3 * r["se"]
    assert r["Q_p"] > 0.001


def test_ivw_matches_closed_form():
    bx, sx, by, sy = _sim_mr(0.2, n=50, seed=3)
    w = bx ** 2 / sy ** 2
    assert ivw(bx, sx, by, sy)["b"] == pytest.approx(np.sum(w * by / bx) / np.sum(w))


def test_egger_detects_directional_pleiotropy():
    bx, sx, by, sy = _sim_mr(0.3, pleio=0.02, n=400, seed=1)
    r = egger(bx, sx, by, sy)
    assert r["int_p"] < 0.01
    assert abs(r["b"] - 0.3) < 4 * r["se"]
    assert ivw(bx, sx, by, sy)["b"] > 0.3   # IVW is biased upward by the planted pleiotropy


def test_weighted_median_robust_to_invalid_minority():
    bx, sx, by, sy = _sim_mr(0.3, n=300, seed=2)
    by = by.copy()
    by[:60] += 0.05                          # 20% invalid instruments
    r = weighted_median(bx, sx, by, sy, nboot=200)
    assert abs(r["b"] - 0.3) < 0.05


def test_presso_flags_outliers():
    bx, sx, by, sy = _sim_mr(0.3, n=100, seed=4)
    by = by.copy()
    by[:3] += 0.4
    r = presso(bx, sx, by, sy, nsim=500)
    assert r["global_p"] < 0.05
    assert set(r["outliers"]) >= {0, 1, 2}


# ---------------- comparison ---------------- #
def test_q_between_null_and_alt():
    assert q_between(np.array([0.3, 0.3, 0.3]), np.array([0.05] * 3))[1] == pytest.approx(1.0)
    assert q_between(np.array([0.1, 0.5]), np.array([0.05, 0.05]))[1] < 1e-6


def test_permutation_detects_cluster_difference():
    a = _sim_mr(0.1, n=80, seed=5)
    b = _sim_mr(0.6, n=80, seed=6)
    bx, sx, by, sy = (np.concatenate([x, y]) for x, y in zip(a, b))
    lab = np.r_[np.zeros(80, int), np.ones(80, int)]
    assert permutation_test(bx, sx, by, sy, lab, B=200)["perm_p"] < 0.01


def test_mvmr_removes_mediated_effect():
    rng = np.random.default_rng(7)
    n = 400
    bx1 = rng.normal(0, 0.1, n)
    bx2 = rng.normal(0, 0.1, n) + 0.5 * bx1
    sy = np.full(n, 0.01)
    by = 0.2 * bx1 + 0.4 * bx2 + rng.normal(0, sy)
    r = mvmr_ivw(np.column_stack([bx1, bx2]), by, sy, ["x1", "x2"])
    assert r["b"].iloc[0] == pytest.approx(0.2, abs=0.03)
    assert r["b"].iloc[1] == pytest.approx(0.4, abs=0.03)


# ---------------- clustering ---------------- #
def test_bnmf_shapes_and_ward():
    rng = np.random.default_rng(0)
    V = np.abs(rng.normal(size=(60, 8)))
    r = bnmf_l2eu(V, rng, n_iter=300)
    assert r.W.shape == (60, r.k) and r.H.shape == (r.k, 8)
    lab = ward(np.vstack([rng.normal(0, 0.1, (20, 3)), rng.normal(5, 0.1, (20, 3))]), 2)
    assert len(set(lab[:20])) == 1 and len(set(lab[20:])) == 1 and lab[0] != lab[-1]


def test_align_labels():
    ref = np.array([0, 0, 1, 1, 2, 2])
    assert (align_labels(ref, np.array([2, 2, 0, 0, 1, 1])) == ref).all()


def test_profile_matrix_per_sbp():
    z = pd.DataFrame(dict(sbp=[10.0, 5.0], bmi=[2.0, 2.0]))
    pm = profile_matrix(z, "per_sbp")
    assert list(pm.columns) == ["bmi"]
    assert pm["bmi"].iloc[1] == pytest.approx(2 * pm["bmi"].iloc[0])
