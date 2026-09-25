"""Check setup: can this computer run Bryan? Said in words that can be sent on.

The launcher runs in its own environment and shells out to Bryan's interpreter,
so the two are checked apart:

- **Bryan's interpreter** is asked, in a subprocess, which of Bryan's packages
  it imports and at which version - the only way to know, since it is another
  Python. A missing scipy shows here, not an hour into a run.
- **Main.py**, and the **model executable** the open sims list's model config
  names (URBS's ``urbs32.exe``), are looked for on disk.
- **The launcher's own packages** are read from its installed metadata; the
  Lake record page's netCDF and shapefile readers are optional.

The result is also written to a text file, so a colleague setting up can send
it rather than describe it.
"""

from __future__ import annotations

import datetime as dt
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .paths import clean_path_text, normalise_sep

OK, WARN, FAIL = "ok", "warn", "fail"

# (module, required, what it is for). Required ones fail the check when missing.
BRYAN_PACKAGES = (
    ("numpy", True, ""),
    ("scipy", True, ""),
    ("pandas", True, ""),
    ("matplotlib", True, "the plots"),
    ("openpyxl", True, "reading the sims lists"),
    ("pyarrow", False, "only for .parquet results"),
)
# (distribution, required, what it is for), read from the launcher's metadata.
UI_PACKAGES = (
    ("nicegui", True, ""),
    ("pandas", True, ""),
    ("openpyxl", True, ""),
    ("psutil", True, "stopping runs"),
    ("netCDF4", False, "the Lake record page's catchment rainfall"),
    ("pyshp", False, "the Lake record page's catchment shapefile"),
    ("pyproj", False, "a catchment shapefile not in latitude and longitude"),
)

REPORT_NAME = "bryan_setup_check.txt"
PROBE_TIMEOUT = 120.0

# Run by Bryan's interpreter: its version, and each package's version or why it
# would not import. Standard library only, and one line of JSON out.
_PROBE = r"""
import importlib, json, sys
found = {}
for name in sys.argv[1:]:
    try:
        module = importlib.import_module(name)
        found[name] = {"version": str(getattr(module, "__version__", "") or "?")}
    except Exception as exc:
        found[name] = {"error": f"{type(exc).__name__}: {exc}"}
print(json.dumps({"python": sys.version.split()[0], "executable": sys.executable,
                  "packages": found}))
"""


@dataclass(frozen=True)
class Check:
    area: str        # "Bryan", "Model" or "Launcher"
    name: str
    status: str      # OK, WARN or FAIL
    detail: str


class ProbeError(Exception):
    """Bryan's interpreter could not be asked."""


def probe_interpreter(python, modules, timeout: float = PROBE_TIMEOUT) -> dict:
    """Ask an interpreter which of ``modules`` it imports, and at what version."""
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    try:
        done = subprocess.run([str(python), "-c", _PROBE, *modules],
                              capture_output=True, text=True, timeout=timeout,
                              creationflags=flags)
    except subprocess.TimeoutExpired:
        raise ProbeError(f"it did not answer within {timeout:.0f} s") from None
    except OSError as exc:
        raise ProbeError(f"it could not be started ({exc})") from exc
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    if done.returncode != 0 or not lines:
        tail = (done.stderr or done.stdout).strip().splitlines()[-3:]
        raise ProbeError("it stopped with an error" + (": " + " / ".join(tail) if tail else ""))
    return json.loads(lines[-1])


def run_checks(*, bryan_python="", bryan_main="", awap_folder="", sims_config=None,
               probe=None) -> list[Check]:
    """Everything a run needs, checked. ``sims_config`` is the open SimsConfig, if any."""
    checks = []
    checks += _interpreter_checks(bryan_python, probe or probe_interpreter)
    checks.append(_main_check(bryan_main))
    checks += _model_checks(sims_config)
    checks += _launcher_checks()
    if awap_folder:
        folder = Path(awap_folder)
        checks.append(Check("Launcher", "AWAP folder",
                            OK if folder.is_dir() else WARN,
                            str(folder) if folder.is_dir() else f"not found: {folder}"))
    return checks


def _interpreter_checks(python, probe) -> list[Check]:
    if not python:
        return [Check("Bryan", "Interpreter", FAIL,
                      "Not set. Give the python.exe Bryan runs with - the VENV_PY of "
                      "the model's batch files.")]
    if not Path(python).exists():
        return [Check("Bryan", "Interpreter", FAIL, f"Not found: {python}")]
    try:
        answer = probe(python, [name for name, _, _ in BRYAN_PACKAGES])
    except ProbeError as exc:
        return [Check("Bryan", "Interpreter", FAIL, f"{python}: {exc}")]
    checks = [Check("Bryan", "Interpreter", OK,
                    f"Python {answer.get('python', '?')}, {answer.get('executable') or python}")]
    packages = answer.get("packages") or {}
    for name, required, purpose in BRYAN_PACKAGES:
        found = packages.get(name) or {"error": "not asked"}
        if "version" in found:
            checks.append(Check("Bryan", name, OK, found["version"]))
            continue
        why = f" ({purpose})" if purpose else ""
        checks.append(Check(
            "Bryan", name, FAIL if required else WARN,
            f"Not installed in Bryan's interpreter{why}. Install it there with: "
            f"\"{python}\" -m pip install {name}.  [{found['error']}]"))
    return checks


def _main_check(main) -> Check:
    if not main:
        return Check("Bryan", "Main.py", FAIL, "Not set.")
    path = Path(main)
    if not path.is_file():
        return Check("Bryan", "Main.py", FAIL, f"Not found: {path}")
    if not (path.parent / "lib").is_dir():
        return Check("Bryan", "Main.py", WARN,
                     f"{path} has no lib folder beside it - is it Bryan's?")
    return Check("Bryan", "Main.py", OK, str(path))


def model_executables(sims_config) -> list[tuple[Path, Path]]:
    """(model config, executable) as the open sims list's model config names it."""
    if sims_config is None:
        return []
    config_path = (sims_config.filepaths or {}).get("model_config")
    if config_path is None or not Path(config_path).is_file():
        return []
    try:
        data = json.loads(Path(config_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    text = clean_path_text(data.get("model_exe") if isinstance(data, dict) else "")
    if not text:
        return []
    exe = Path(normalise_sep(text))
    if not exe.is_absolute():
        exe = Path(config_path).parent / exe
    return [(Path(config_path), exe)]


def _model_checks(sims_config) -> list[Check]:
    if sims_config is None:
        return [Check("Model", "Model executable", WARN,
                      "No sims list is open, so the model config that names it was "
                      "not read.")]
    found = model_executables(sims_config)
    if not found:
        return [Check("Model", "Model executable", WARN,
                      "The open sims list's model config names no model_exe.")]
    return [Check("Model", exe.name, OK if exe.is_file() else FAIL,
                  (str(exe) if exe.is_file() else f"Not found: {exe}")
                  + f"  (from {config.name})")
            for config, exe in found]


def _launcher_checks() -> list[Check]:
    checks = [Check("Launcher", "Interpreter", OK,
                    f"Python {platform.python_version()}, {sys.executable}")]
    for name, required, purpose in UI_PACKAGES:
        try:
            checks.append(Check("Launcher", name, OK, importlib.metadata.version(name)))
        except importlib.metadata.PackageNotFoundError:
            why = f" - needed for {purpose}" if purpose else ""
            checks.append(Check("Launcher", name, FAIL if required else WARN,
                                f"Not installed{why}. Install it with: "
                                f"\"{sys.executable}\" -m pip install {name}"))
    return checks


def worst(checks) -> str:
    statuses = {check.status for check in checks}
    return FAIL if FAIL in statuses else WARN if WARN in statuses else OK


def report_text(checks, when: dt.datetime | None = None) -> str:
    when = when or dt.datetime.now()
    verdict = {OK: "Ready to run.", WARN: "Ready to run, with warnings.",
               FAIL: "Not ready: fix the FAIL lines."}[worst(checks)]
    lines = [f"Bryan setup check, {when:%Y-%m-%d %H:%M}, on {platform.node()} "
             f"({platform.platform()})", verdict, ""]
    area = None
    for check in checks:
        if check.area != area:
            if area is not None:
                lines.append("")
            area = check.area
            lines.append(f"[{area}]")
        lines.append(f"  {check.status.upper():4}  {check.name}: {check.detail}")
    return "\n".join(lines) + "\n"


def write_report(checks, folder) -> Path:
    path = Path(folder) / REPORT_NAME
    path.write_text(report_text(checks), encoding="utf-8")
    return path
