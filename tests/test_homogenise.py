"""Lake level homogenisation (lib/homogenise), generalised from callide-fsl-reinstate.

The engine is the callide package's own modules, copied; what is new is the job
layer that takes every input and setting from a job file. So the tests here are
of that layer, on small synthetic files, and one that holds the whole thing to
the original: with the callide-fsl-reinstate checkout beside this one, the
generalised job run on Callide's data must write the same files, byte for byte,
as ``callide.cli.run_scenario`` does.
"""

from __future__ import annotations

import filecmp
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from lib.homogenise import evaporation, gauges, job as jobs

BRYAN_ROOT = Path(__file__).resolve().parents[1]
CALLIDE = BRYAN_ROOT.parent / "callide-fsl-reinstate"


def wmip(path, rows, gauge="G1"):
    """A WMIP export: four header rows, data, a licence footer."""
    lines = [f"{gauge},,", "Level,,", "m,,", "Point,Value,Quality"]
    lines += [f"{stamp:%H:%M:%S %d/%m/%Y},{level},9" for stamp, level in rows]
    lines += ["", "Licence footer,,"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def hours(start, levels, step="1h"):
    stamps = pd.date_range(start, periods=len(levels), freq=step)
    return list(zip(stamps, levels))


# -- the gauge chain ---------------------------------------------------------------

def test_a_chain_of_files_hands_over_at_each_gauge_s_first_reading(tmp_path):
    first = wmip(tmp_path / "A.csv", hours("2000-01-01", [100.0] * 10), "A")
    second = wmip(tmp_path / "C.csv", hours("2000-01-01 05:00", [101.0] * 10), "C")
    record = gauges.read_chain([first, second])
    assert record.loc["2000-01-01 04:00", "Gauge"] == "A"
    assert record.loc["2000-01-01 05:00", "Gauge"] == "C"
    assert not record.index.has_duplicates


def test_a_chain_out_of_order_is_refused(tmp_path):
    first = wmip(tmp_path / "A.csv", hours("2001-01-01", [100.0] * 3), "A")
    second = wmip(tmp_path / "C.csv", hours("2000-01-01", [100.0] * 3), "C")
    with pytest.raises(ValueError, match="chronological"):
        gauges.read_chain([first, second])


def test_a_gauge_named_twice_is_refused(tmp_path):
    first = wmip(tmp_path / "A.csv", hours("2000-01-01", [100.0] * 3), "A")
    with pytest.raises(ValueError, match="twice"):
        gauges.read_chain([first, first])


# -- evaporation --------------------------------------------------------------------

def silo(path, days=40, pan=5.0):
    lines = ["preamble", "Date Day Date2 T.Max Smx T.Min Smn Rain Srn Evap Sev",
             "(yyyymmdd) () (ddmmyyyy) (oC) () (oC) () (mm) () (mm) ()"]
    for stamp in pd.date_range("2000-01-01", periods=days, freq="D"):
        lines.append(f"{stamp:%Y%m%d} 1 {stamp:%d-%m-%Y} 30 25 20 25 0 25 {pan} 25")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_pan_factors_are_a_setting(tmp_path):
    path = silo(tmp_path / "silo.txt")
    default = evaporation.read_evaporation(path)
    custom = evaporation.read_evaporation(path, {month: 1.0 for month in range(1, 13)})
    assert default.daily_mm.iloc[0] == pytest.approx(5.0 * 0.82)       # January, Callide
    assert custom.daily_mm.iloc[0] == pytest.approx(5.0)


def test_pan_factors_need_all_twelve_months(tmp_path):
    with pytest.raises(ValueError, match="each month"):
        evaporation.read_evaporation(silo(tmp_path / "silo.txt"), {1: 0.8})


# -- the job ------------------------------------------------------------------------

@pytest.mark.parametrize("missing, message", [
    ("gauges", "no gauge exports"), ("storage", "no storage file"),
    ("targets", "no target rating")])
def test_a_job_says_what_it_is_missing(missing, message):
    settings = {"gauges": ["a.csv"], "storage": "s.els", "register": "r.xlsx",
                "evaporation": "e.txt", "targets": [{"rating": "t.sq"}]}
    settings[missing] = [] if missing in ("gauges", "targets") else ""
    with pytest.raises(jobs.JobError, match=message):
        jobs.from_dict(settings, ".")


def test_relative_paths_resolve_against_the_job_file(tmp_path):
    (tmp_path / "job.json").write_text(json.dumps({
        "gauges": ["g/a.csv"], "storage": "s.els", "register": "r.xlsx",
        "evaporation": "e.txt", "targets": [{"rating": "t.sq"}]}), encoding="utf-8")
    job = jobs.load(tmp_path / "job.json")
    assert job.path("storage") == tmp_path / "s.els"
    assert job.resolve("g/a.csv") == tmp_path / "g" / "a.csv"
    assert job.out == tmp_path / "homogenised"


def test_a_missing_input_file_is_named(tmp_path):
    job = jobs.from_dict({"gauges": ["a.csv"], "storage": "nope.els", "register": "r.xlsx",
                          "evaporation": "e.txt", "targets": [{"rating": "t.sq"}]}, tmp_path)
    with pytest.raises(jobs.JobError, match="storage file not found"):
        jobs.load_inputs(job)


# -- held to the original ----------------------------------------------------------

@pytest.mark.skipif(not (CALLIDE / "callide" / "pipeline.py").is_file(),
                    reason="needs the callide-fsl-reinstate checkout beside Bryan")
def test_callide_homogenises_to_the_byte_as_the_original_did(tmp_path):
    data = CALLIDE / "data"
    job = {"gauges": [str(data / "gauge" / "130314A.csv"), str(data / "gauge" / "130314C.csv")],
           "overlay": {"file": str(data / "gauge" / "130314B_hydstra.csv"), "below": 200.90,
                       "reconnect_margin": 0.10},
           "storage": str(data / "storage" / "CALLIDE_STORAGE.els"),
           "register": str(data / "ratings" / "RatingCurves.xlsx"),
           "evaporation": str(data / "climate" / "silo_-24.35_150.65.txt"),
           "targets": [{"name": "rfsl", "rating": str(data / "ratings" / "CALLIDE_RFSL.sq"),
                        "fsl": 215.5}],
           "out": str(tmp_path / "bryan")}
    summary = jobs.run(jobs.from_dict(job, tmp_path))
    assert summary["record"]["steps"] == 508_179

    python = CALLIDE / ".venv" / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    subprocess.run([str(python if python.exists() else sys.executable), "-m",
                    "callide.cli.run_scenario", "--scenario", "rfsl-215-5-rapid",
                    "--outdir", str(tmp_path / "callide")],
                   cwd=CALLIDE, check=True, capture_output=True, timeout=900)
    original = tmp_path / "callide" / "rfsl-215-5-rapid"
    for name in ("ams.csv", "daily.csv", "peaks.csv"):
        assert filecmp.cmp(original / name, tmp_path / "bryan" / "rfsl" / name, shallow=False), name
    with gzip.open(original / "trace.csv.gz") as theirs, \
            gzip.open(tmp_path / "bryan" / "rfsl" / "trace.csv.gz") as ours:
        assert theirs.read() == ours.read()
