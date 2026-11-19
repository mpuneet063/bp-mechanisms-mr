"""Tests for the OpenGWAS fetch script: study choice, API batch limit, cache reuse, token expiry."""
import base64
import importlib.util
import json
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("fetch", ROOT / "scripts" / "00_fetch_opengwas.py")
fetch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch)

# A catalogue shaped like the real one, including the entries that were mis-picked on 2026-11-18
CATALOG = pd.DataFrame([
    dict(id="ukb-b-9270", trait="Diagnoses - main ICD10: O68.0 Labour and delivery complicated by foetal heart rate anomaly",
         author="Ben Elsworth", year=2018, population="European", sample_size=463010, ncase=500, nsnp=9851867),
    dict(id="ukb-b-18103", trait="Pulse rate", author="Ben Elsworth", year=2018, population="European",
         sample_size=461117, ncase=None, nsnp=9851867),
    dict(id="ebi-a-GCST006906", trait="Ischemic stroke (small-vessel)", author="Malik R", year=2018,
         population="European", sample_size=198048, ncase=4373, nsnp=7000000),
    dict(id="ebi-a-GCST90013868", trait="Coronary artery disease (SPA correction)", author="Mbatchou J", year=2021,
         population="European", sample_size=352063, ncase=None, nsnp=11000000),
    dict(id="ebi-a-GCST005195", trait="Coronary artery disease", author="van der Harst P", year=2017,
         population="European", sample_size=547261, ncase=122733, nsnp=7900000),
    dict(id="ieu-a-7", trait="Coronary heart disease / coronary artery disease", author="Nikpay", year=2015,
         population="Mixed", sample_size=184305, ncase=60801, nsnp=9455779),
    dict(id="ebi-a-GCST90000065", trait="Pulse pressure", author="Surendran P", year=2020, population="European",
         sample_size=810865, ncase=None, nsnp=250000),
    dict(id="ukb-b-pp", trait="Pulse pressure", author="Ben Elsworth", year=2018, population="European",
         sample_size=400000, ncase=None, nsnp=9851867),
]).set_index("id")


def pick(key):
    sp = next(s for s in fetch.STUDY_SPECS if s["key"] == key)
    return fetch.resolve(sp, CATALOG)[0]


def test_heart_rate_never_picks_foetal_diagnosis():
    assert pick("hr") == "ukb-b-18103"


def test_small_vessel_stroke_matches_hyphenated_megastroke_name():
    assert pick("stroke_svs") == "ebi-a-GCST006906"


def test_cad_prefers_consortium_meta_analysis():
    assert pick("cad") == "ebi-a-GCST005195"


def test_pulse_pressure_skips_exome_array():
    assert pick("pp") == "ukb-b-pp"


def test_chunk_respects_api_limit():
    assert fetch.MAX_ASSOC_PER_REQUEST == 64


def _token(exp):
    body = base64.urlsafe_b64encode(json.dumps({"exp": exp}).encode()).decode().rstrip("=")
    return f"xx.{body}.yy"


def test_token_expiry_decoded():
    t = time.time() + 3600
    assert fetch.jwt_expiry(_token(t)) == pytest.approx(t)
    assert fetch.jwt_expiry("not-a-jwt") is None


def test_cache_reused_across_chunk_sizes(tmp_path, monkeypatch):
    monkeypatch.setattr(fetch, "CACHE", tmp_path)
    # results cached by an older run (100-variant chunks bisected into 50s), one via a proxy
    rows = [dict(id="s1", rsid=f"rs{i}", proxy=False) for i in range(50)]
    rows.append(dict(id="s1", rsid="rs999", proxy=True, target_snp="rs50"))
    (tmp_path / "assoc_s1_abcdef123456.json").write_text(json.dumps(rows))
    fetch._record_done("s1", ["rs51", "rs52"])  # requested, but no association in the study
    got, done = fetch._cached_study("s1")
    assert len(got) == 51
    assert {f"rs{i}" for i in range(53)} <= done


def test_proxy_rows_keyed_by_target(tmp_path):
    from bpmr.data import load_raw
    pd.DataFrame([dict(key="sbp", role="exposure", id="e")]).to_csv(tmp_path / "manifest.csv", index=False)
    pd.DataFrame([dict(rsid="rs1", source="x")]).to_csv(tmp_path / "instruments.csv", index=False)
    pd.DataFrame([dict(id="e", rsid="rs777", ea="A", nea="G", eaf=0.3, beta=0.1, se=0.01, p=1e-5, n=1,
                       proxy=True, target_snp="rs1")]).to_csv(tmp_path / "associations.csv.gz", index=False)
    assert load_raw(tmp_path).assoc["rsid"].tolist() == ["rs1"]
