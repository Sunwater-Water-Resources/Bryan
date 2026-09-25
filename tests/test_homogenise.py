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

import numpy as np
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


def test_an_analysis_takes_its_own_water_year_or_the_homogenisation_s():
    """The antecedent storage and the inflow record may each override the water
    year the launcher shares across the study; with none given, they follow the
    homogenisation, as they always did."""
    homogenise = {"water_year_start": 10}
    assert jobs.water_year_start({}, homogenise) == 10
    assert jobs.water_year_start({"water_year_start": None}, homogenise) == 10
    assert jobs.water_year_start({"water_year_start": 7}, homogenise) == 7
    assert jobs.water_year_start({"water_year_start": "7"}, homogenise) == 7


# -- one rating for the whole record ------------------------------------------------------

from lib.homogenise import curves                                  # noqa: E402


def write_rat(path, pairs, *, fsl_note=None):
    lines = ["KROOMBIT OUTFLOW", "* written for a test"]
    if fsl_note:
        lines.append(f"* {fsl_note}")
    lines.append(f"{len(pairs)} PAIRS:")
    lines += [f"{level:.3f} {flow:.2f}" for level, flow in pairs]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


PAIRS = [(265.8, 0.0), (266.0, 60.0), (266.5, 350.0), (267.0, 900.0), (268.0, 2600.0)]


def test_a_rat_is_read_as_level_and_flow_with_the_header_s_full_supply(tmp_path):
    table, comments = curves.read_rat(write_rat(tmp_path / "k.rat", PAIRS,
                                                fsl_note="FSL is 265.8 mAHD"))
    assert list(table.columns) == ["level", "flow"] and len(table) == 5
    assert curves.declared_sq_fsl(comments) == 265.8
    rating = curves.load_rating(tmp_path / "k.rat", None)
    assert rating.fsl == 265.8 and rating.at(265.8) == 0.0
    assert rating.at(267.0) == pytest.approx(900.0)


def test_a_rat_s_full_supply_has_to_agree_with_what_is_given(tmp_path):
    path = write_rat(tmp_path / "k.rat", PAIRS, fsl_note="FSL is 265.8 mAHD")
    with pytest.raises(ValueError, match="declares a full supply level of 265.800"):
        curves.load_rating(path, None, fsl=266.0)
    starts_low = write_rat(tmp_path / "low.rat", [(265.0, 0.0)] + PAIRS[1:],
                           fsl_note="FSL is 265.8 mAHD")
    with pytest.raises(ValueError, match="starts at 265.000"):
        curves.load_rating(starts_low, None)


def test_one_rating_becomes_a_register_in_force_over_any_record(tmp_path):
    register, ratings = curves.read_rating_source(write_rat(tmp_path / "k.rat", PAIRS), None)
    assert list(register.columns) == ["from", "to", "FSL"] and len(register) == 1
    assert register["from"].iloc[0].year == 1800 and register["to"].iloc[0].year == 2200
    assert register["FSL"].iloc[0] == 265.8 and list(ratings) == [1]


def test_a_sq_without_its_full_supply_needs_one_given(tmp_path):
    sq = tmp_path / "k.sq"
    sq.write_text("KROOMBIT\n*\n2 PAIRS:\n0 0\n586 52\n", encoding="utf-8")
    volume, _ = curves.read_storage(write_els(tmp_path / "k.els"))
    with pytest.raises(ValueError, match="not stated in the header"):
        curves.read_rating_source(sq, volume)
    register, _ = curves.read_rating_source(sq, volume, fsl=265.8)
    assert register["FSL"].iloc[0] == 265.8


def write_els(path):
    levels = np.arange(250.0, 272.0, 0.5)
    frame = pd.DataFrame({"EL": levels, "A": 1000.0, "V": (levels - 250.0) * 10_000.0})
    frame.to_csv(path, index=False)
    return path


KROOMBIT_RECORD = BRYAN_ROOT.parent / "callide-fsl-reinstate" / "data" / "gauge" / "130360A_level_point.csv"
KROOMBIT_RATINGS = BRYAN_ROOT.parent / "callide-design-flood-hydrology" / "runs" / "Regional_E001" / "ratings"
SILO = BRYAN_ROOT.parent / "callide-fsl-reinstate" / "data" / "climate" / "silo_-24.35_150.65.txt"


@pytest.mark.skipif(not (KROOMBIT_RECORD.is_file() and (KROOMBIT_RATINGS / "kroombit.sq").is_file()
                         and SILO.is_file()),
                    reason="needs the callide-fsl-reinstate and callide-design-flood-hydrology "
                           "checkouts beside Bryan")
def test_kroombit_s_single_rat_and_a_one_row_register_give_the_same_inflow(tmp_path):
    """Phase 3's check: one rating read from a .rat is the same homogenisation input
    as the same rating in a register workbook, on Kroombit's own record."""
    fsl = 265.8
    els = pd.read_csv(KROOMBIT_RATINGS / "kroombit.els")
    sq, _ = curves.read_sq(KROOMBIT_RATINGS / "kroombit.sq")
    fsv = float(np.interp(fsl, els["EL"], els["V"]))
    levels = np.interp(fsv + sq["storage_ML"], els["V"], els["EL"])
    # Rounded as a .rat file holds it, so the workbook gets the very same rating.
    table = pd.DataFrame({"level": np.round(levels, 3), "flow": np.round(sq["flow_m3s"], 2)}) \
        .drop_duplicates("level").reset_index(drop=True)
    table.loc[0, "level"] = fsl
    rat = write_rat(tmp_path / "kroombit.rat", list(zip(table["level"], table["flow"])),
                    fsl_note=f"FSL is {fsl} mAHD")
    workbook = tmp_path / "register.xlsx"
    with pd.ExcelWriter(workbook) as writer:
        pd.DataFrame({"Rating": [1], "from": ["1/1/1990"], "to": ["1/1/2030"], "FSL": [fsl]}) \
            .to_excel(writer, sheet_name="Register", index=False)
        table.to_excel(writer, sheet_name="1", index=False)

    def derived(register):
        job = jobs.from_dict({"gauges": [str(KROOMBIT_RECORD)],
                              "storage": str(KROOMBIT_RATINGS / "kroombit.els"),
                              "register": str(register), "evaporation": str(SILO),
                              "targets": []}, tmp_path, routing=False)
        return jobs.load_inputs(job).dam().derive_inflow()

    single, workbook_rated = derived(rat), derived(workbook)
    assert len(single) == len(workbook_rated) > 10_000
    pd.testing.assert_series_equal(single["Inflow_ML"], workbook_rated["Inflow_ML"])
    pd.testing.assert_series_equal(single["Release_ML"], workbook_rated["Release_ML"])
    assert single["Release_ML"].sum() > 0                  # Kroombit spills in the record
