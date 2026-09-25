"""Check setup: what Bryan's interpreter imports, and whether the files are there."""

from __future__ import annotations

import json
import sys

import pytest

from core import setupcheck
from core.config import load_sims_config


def _answer(missing=()):
    """What the probe of an interpreter lacking ``missing`` would say."""
    def probe(python, modules):
        return {"python": "3.12.1", "executable": str(python),
                "packages": {name: ({"error": f"ModuleNotFoundError: No module named "
                                              f"'{name}'"} if name in missing
                                    else {"version": "1.0"}) for name in modules}}
    return probe


@pytest.fixture
def bryan(tmp_path):
    root = tmp_path / "Bryan"
    (root / "lib").mkdir(parents=True)
    (root / "Main.py").write_text("", encoding="utf-8")
    python = root / "python.exe"
    python.write_text("", encoding="utf-8")
    return {"bryan_python": str(python), "bryan_main": str(root / "Main.py")}


def _sims_config(folder, model_exe):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "model.json").write_text(json.dumps({"model_exe": model_exe}),
                                       encoding="utf-8")
    (folder / "sims.xlsx").write_text("", encoding="utf-8")
    config = folder / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx",
                                  "filepaths": {"model_config": "model.json"}}),
                      encoding="utf-8")
    return load_sims_config(config)


def _by_name(checks):
    return {(check.area, check.name): check for check in checks}


def test_a_missing_scipy_is_said_in_plain_words(bryan):
    checks = setupcheck.run_checks(**bryan, probe=_answer(missing={"scipy"}))
    scipy = _by_name(checks)[("Bryan", "scipy")]
    assert scipy.status == setupcheck.FAIL
    assert "Not installed in Bryan's interpreter" in scipy.detail
    assert "-m pip install scipy" in scipy.detail
    assert setupcheck.worst(checks) == setupcheck.FAIL
    assert "Not ready" in setupcheck.report_text(checks)


def test_pyarrow_is_only_a_warning(bryan):
    checks = setupcheck.run_checks(**bryan, probe=_answer(missing={"pyarrow"}))
    assert _by_name(checks)[("Bryan", "pyarrow")].status == setupcheck.WARN


def test_everything_present_reports_the_versions(bryan, tmp_path):
    exe = tmp_path / "urbs32.exe"
    exe.write_text("", encoding="utf-8")
    config = _sims_config(tmp_path / "model", str(exe))
    checks = setupcheck.run_checks(**bryan, sims_config=config, probe=_answer())
    found = _by_name(checks)
    assert found[("Bryan", "numpy")].detail == "1.0"
    assert found[("Model", "urbs32.exe")].status == setupcheck.OK
    assert "from model.json" in found[("Model", "urbs32.exe")].detail
    bryan_and_model = [check for check in checks if check.area != "Launcher"]
    assert setupcheck.worst(bryan_and_model) == setupcheck.OK


def test_a_missing_model_executable_fails(bryan, tmp_path):
    config = _sims_config(tmp_path / "model", r"C:\nowhere\urbs32.exe")
    found = _by_name(setupcheck.run_checks(**bryan, sims_config=config, probe=_answer()))
    assert found[("Model", "urbs32.exe")].status == setupcheck.FAIL


def test_without_a_sims_list_the_executable_is_not_judged(bryan):
    found = _by_name(setupcheck.run_checks(**bryan, probe=_answer()))
    assert found[("Model", "Model executable")].status == setupcheck.WARN


@pytest.mark.parametrize("change, expect", [
    ({"bryan_python": ""}, "Not set"),
    ({"bryan_python": "C:/nowhere/python.exe"}, "Not found"),
])
def test_an_interpreter_not_there_is_not_asked(bryan, change, expect):
    def probe(*_):
        raise AssertionError("should not be asked")
    found = _by_name(setupcheck.run_checks(**{**bryan, **change}, probe=probe))
    assert found[("Bryan", "Interpreter")].status == setupcheck.FAIL
    assert expect in found[("Bryan", "Interpreter")].detail


def test_main_py_without_lib_beside_it_is_questioned(bryan, tmp_path):
    other = tmp_path / "elsewhere" / "Main.py"
    other.parent.mkdir()
    other.write_text("", encoding="utf-8")
    found = _by_name(setupcheck.run_checks(**{**bryan, "bryan_main": str(other)},
                                           probe=_answer()))
    assert found[("Bryan", "Main.py")].status == setupcheck.WARN


def test_the_probe_really_asks_an_interpreter():
    answer = setupcheck.probe_interpreter(sys.executable, ["json", "no_such_module_x"])
    assert answer["packages"]["json"]["version"]
    assert "ModuleNotFoundError" in answer["packages"]["no_such_module_x"]["error"]


def test_a_broken_interpreter_is_a_probe_error(tmp_path):
    with pytest.raises(setupcheck.ProbeError):
        setupcheck.probe_interpreter(tmp_path / "not_python.exe", ["json"])


def test_the_report_is_written_to_a_file(bryan, tmp_path):
    checks = setupcheck.run_checks(**bryan, probe=_answer())
    path = setupcheck.write_report(checks, tmp_path)
    text = path.read_text(encoding="utf-8")
    assert path.name == setupcheck.REPORT_NAME
    assert "[Bryan]" in text and "[Launcher]" in text and "numpy: 1.0" in text
