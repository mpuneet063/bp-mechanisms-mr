"""Single place for paths and analysis settings. Every experiment script imports from here."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

MAIN_OUTCOMES = ["stroke_any", "stroke_is", "stroke_las", "stroke_ces", "stroke_svs", "ich", "cad"]
REPLICATION = ["fg_stroke", "fg_is", "fg_chd"]
MVMR_MEDIATORS = ["bmi", "ldl", "t2d"]
OUTCOME_LABELS = {
    "stroke_any": "Any stroke", "stroke_is": "Ischaemic stroke", "stroke_las": "Large-artery stroke",
    "stroke_ces": "Cardioembolic stroke", "stroke_svs": "Small-vessel stroke", "ich": "Intracerebral haemorrhage",
    "cad": "Coronary artery disease", "fg_stroke": "FinnGen stroke", "fg_is": "FinnGen ischaemic stroke",
    "fg_chd": "FinnGen CHD",
}


@dataclass
class Settings:
    raw: str = "data/raw"
    out: str = "results"
    seed: int = 20261114
    # matrix
    max_missing_variant: float = 0.2
    max_missing_trait: float = 0.5
    corr_cutoff: float = 0.85
    min_sbp_z_cluster: float = 4.0  # clustering variants need |z_SBP| >= 4 in the exposure GWAS (per-SBP profile stays stable)
    # clustering
    profile: str = "per_sbp"         # clustering input: per_sbp | raw_no_sbp | raw  (see bpmr.matrix.profile_matrix)
    bnmf_phi: float = 0.5            # bNMF regulariser (Udler's phi); scan reported in 02_cluster/bnmf_phi_scan.csv
    bnmf_phi_scan: list = field(default_factory=lambda: [0.25, 0.5, 1.0])
    bnmf_phi_scan_runs: int = 20
    primary: str = "hclust"          # partition used for the headline MR: "hclust" (Pascat 2026) or "bnmf"
    k_override: int = 0              # 0 = K from bNMF (modal K across runs)
    bnmf_runs: int = 100
    bnmf_k_init: int = 20
    bnmf_k0: int = 10
    hclust_clip_sd: float = 2.0
    # stability
    boot_variant: int = 500
    boot_noise: int = 200
    boot_k_scan: int = 100
    bnmf_noise_runs: int = 20
    jaccard_stable: float = 0.75
    jaccard_min: float = 0.6         # pre-specified: clusters below this are flagged unstable
    # MR
    min_cluster_snps: int = 10
    f_min: float = 10.0
    wm_boot: int = 1000
    presso_sim: int = 2000
    perm_B: int = 2000
    pair_B: int = 2000
    outcomes: list = field(default_factory=lambda: MAIN_OUTCOMES + REPLICATION)

    @property
    def outp(self) -> Path:
        return Path(self.out)

    def step(self, name: str) -> Path:
        p = self.outp / name
        p.mkdir(parents=True, exist_ok=True)
        return p


FAST = dict(bnmf_runs=20, bnmf_phi_scan_runs=5, boot_variant=60, boot_noise=40, boot_k_scan=20, bnmf_noise_runs=4, wm_boot=200,
            presso_sim=300, perm_B=200, pair_B=200)


def parse(description: str) -> Settings:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--raw", default=None, help="folder with manifest.csv, associations.csv.gz, instruments.csv")
    ap.add_argument("--out", default=None, help="results folder")
    ap.add_argument("--fast", action="store_true", help="small bootstrap/run counts for a quick end-to-end check")
    ap.add_argument("--primary", choices=["hclust", "bnmf"], default=None)
    ap.add_argument("--k", type=int, default=None, help="force the number of clusters")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--profile", choices=["per_sbp", "raw_no_sbp", "raw"], default=None)
    ap.add_argument("--phi", type=float, default=None, help="bNMF regulariser")
    a = ap.parse_args()
    s = Settings()
    if a.fast:
        for k, v in FAST.items():
            setattr(s, k, v)
    for k in ["raw", "out", "primary", "seed"]:
        if getattr(a, k) is not None:
            setattr(s, k, getattr(a, k))
    if a.profile is not None:
        s.profile = a.profile
    if a.phi is not None:
        s.bnmf_phi = a.phi
    if a.k is not None:
        s.k_override = a.k
    s.outp.mkdir(parents=True, exist_ok=True)
    (s.outp / "settings.json").write_text(json.dumps(asdict(s), indent=2))
    return s
