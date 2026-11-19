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


def test_parse_vcf_stream_keeps_only_wanted_and_reads_format():
    lines = [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.2:0.05:4:0.3:1000:rs1\n",
        "1\t200\trs2\tC\tT\t.\tPASS\t.\tSE:ES:LP:AF:ID\t0.01:-0.1:2:0.6:rs2\n",   # other FORMAT order, no SS
        "1\t300\trs3\tA\tG,T\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.01:5:0.2:900:rs3\n",  # multi-allelic: skipped
        "1\t400\trs4\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.01:5:0.2:900:rs4\n",   # not wanted
    ]
    rows = fetch.parse_vcf_stream(lines, {"rs1", "rs2", "rs3"}, "s1", default_n=555)
    assert [r["rsid"] for r in rows] == ["rs1", "rs2"]
    assert rows[0]["ea"] == "G" and rows[0]["nea"] == "A" and rows[0]["p"] == pytest.approx(1e-4)
    assert rows[1]["beta"] == pytest.approx(-0.1) and rows[1]["n"] == 555


def test_t2d_takes_largest_european_study():
    cat = pd.DataFrame([
        dict(id="m14", trait="Type 2 diabetes", author="Mahajan A", year=2014, population="Mixed", sample_size=110452, ncase=26488, nsnp=2.5e6),
        dict(id="xue", trait="Type 2 diabetes", author="Xue A", year=2018, population="European", sample_size=655666, ncase=62892, nsnp=5e6),
        dict(id="small", trait="Type 2 diabetes", author="X", year=2020, population="European", sample_size=50000, ncase=5000, nsnp=8e6),
    ]).set_index("id")
    sp = next(s for s in fetch.STUDY_SPECS if s["key"] == "t2d")
    assert fetch.resolve(sp, cat)[0] == "xue"


def test_t2d_case_control_beats_linear_biobank_release():
    # 2026-11-19: Loh 2018 (BOLT-LMM, no case count, linear scale) outranked Xue 2018 because its total N
    # was compared with Xue's case count
    cat = pd.DataFrame([
        dict(id="loh", trait="Type 2 diabetes", author="Loh PR", year=2018, population="European", sample_size=468298, ncase=None, nsnp=1.1e7),
        dict(id="xue", trait="Type 2 diabetes", author="Xue A", year=2018, population="European", sample_size=655666, ncase=61714, nsnp=5e6),
        dict(id="sak", trait="Type 2 diabetes", author="Sakaue S", year=2021, population="European", sample_size=490089, ncase=38841, nsnp=1.3e7),
    ]).set_index("id")
    sp = next(s for s in fetch.STUDY_SPECS if s["key"] == "t2d")
    sid, _, cands = fetch.resolve(sp, cat)
    assert sid == "xue" and list(cands.index) == ["xue", "sak", "loh"]


def test_pinned_study_keeps_ranked_candidates_for_fallback():
    cat = pd.DataFrame([
        dict(id="ieu-b-40", trait="body mass index", author="Yengo, L", year=2018, population="European", sample_size=681275, ncase=None, nsnp=2.3e6),
        dict(id="ukb-b-19953", trait="Body mass index (BMI)", author="Ben Elsworth", year=2018, population="European", sample_size=461460, ncase=None, nsnp=9.8e6),
    ]).set_index("id")
    sp = next(s for s in fetch.STUDY_SPECS if s["key"] == "bmi")
    sid, how, cands = fetch.resolve(sp, cat)
    assert (sid, how) == ("ukb-b-19953", "pinned")
    assert list(cands.index) == ["ukb-b-19953", "ieu-b-40"]


def test_whr_fallback_needs_bmi_adjustment_in_trait_or_note():
    cat = pd.DataFrame([
        dict(id="barton", trait="Waist-to-hip ratio adjusted for BMI", author="Barton AR", year=2021, population="European", sample_size=458349, ncase=None, nsnp=2e6, note=""),
        dict(id="giant_unadj", trait="Waist-to-hip ratio", author="Shungin D", year=2015, population="European", sample_size=212244, ncase=None, nsnp=2.5e6, note=""),
        dict(id="giant_adj", trait="Waist-to-hip ratio", author="Shungin D", year=2015, population="European", sample_size=210088, ncase=None, nsnp=2.5e6, note="Adjusted for BMI"),
        dict(id="mixed", trait="Waist-to-hip ratio adjusted for BMI", author="Fernandez-Rhodes L", year=2022, population="Mixed", sample_size=205581, ncase=None, nsnp=8e6, note=""),
    ]).set_index("id")
    sp = next(s for s in fetch.STUDY_SPECS if s["key"] == "whradjbmi")
    sid, _, cands = fetch.resolve(sp, cat)
    assert sid == "barton" and list(cands.index) == ["barton", "giant_adj"]


VCF = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"] + [
    f"1\t{100 * i}\trs{i}\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.01:5:0.2:900:rs{i}\n" for i in range(1, 11)]


def test_extract_reads_only_new_variants_and_honours_legacy_checkpoints(tmp_path, monkeypatch):
    import gzip
    monkeypatch.setattr(fetch, "CACHE", tmp_path / "cache")
    (tmp_path / "cache").mkdir()
    keep = tmp_path / "files"
    keep.mkdir()
    with gzip.open(keep / "s1.vcf.gz", "wt") as fh:
        fh.writelines(VCF)
    # a checkpoint written by the old script: rows for rs1-rs3, searched for rs1-rs5, no sidecar
    pd.DataFrame(fetch.parse_vcf_stream(VCF, {"rs1", "rs2", "rs3"}, "s1")).to_csv(
        tmp_path / "cache" / "extract_s1.csv.gz", index=False, compression="gzip")
    legacy = {f"rs{i}" for i in range(1, 6)}
    rows, status = fetch.extract_study(None, "s1", {"rs1", "rs2"}, keep, legacy_searched=legacy)
    assert status == "extract already saved" and len(rows) == 3
    # rs4/rs5 were searched before (legacy) -> not re-read; rs7 is new -> read and appended
    rows, status = fetch.extract_study(None, "s1", {"rs1", "rs4", "rs5", "rs7"}, keep, legacy_searched=legacy)
    assert "+1 new variants" in status
    assert sorted(r["rsid"] for r in rows) == ["rs1", "rs2", "rs3", "rs7"]
    assert fetch._searched_path("s1").exists()
    rows, status = fetch.extract_study(None, "s1", {"rs7"}, keep)  # sidecar now authoritative
    assert status == "extract already saved"


def test_fetch_by_files_measures_coverage_on_bp_variants_only(tmp_path, monkeypatch):
    import gzip
    monkeypatch.setattr(fetch, "CACHE", tmp_path / "cache")
    (tmp_path / "cache").mkdir()
    keep = tmp_path / "files"
    keep.mkdir()
    with gzip.open(keep / "s1.vcf.gz", "wt") as fh:
        fh.writelines(VCF[:4])  # rs1-rs3 only
    core = {"rs1", "rs2", "rs3", "rs4"}
    got, rep, dead = fetch.fetch_by_files(None, {"s1": core | {"rs99"}}, core, CATALOG, workers=1, keep_dir=keep)
    assert not dead and rep[0]["found"] == 3 and rep[0]["missing_pct"] == 25.0
