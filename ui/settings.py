"""Where Bryan is, and how the UI behaves. Persisted to the user's home.

The UI runs in its own environment - Streamlit-style all-in-one installs would
mean adding dependencies to the frozen ``bryan29`` conda env that reproduces
study results. So the interpreter and Main.py are configuration, mirroring the
``PROJECT_ROOT`` / ``VENV_PY`` pair the existing per-model .bat files set.
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

from core.bryan import BRYAN_ROOT
from core.paths import atomic_write_json, read_json

SETTINGS_PATH = Path.home() / ".bryan_ui.json"
MAX_RECENT = 10


@dataclass
class UiSettings:
    bryan_python: str = ""
    bryan_main: str = ""
    max_parallel: int = 1          # see ui/README.md - parallelism rarely pays
    poll_seconds: float = 2.0
    show_consoles: bool = False
    keep_groups_together: bool = True
    recent_configs: list = field(default_factory=list)
    # The Downstream page's two config paths, per project. Keyed on the resolved absolute
    # path of the sims config, never its name: two projects routinely carry the same config
    # filename at different paths (one per model revision), and keying on the name would
    # hand one project's regional model to the other.
    downstream_configs: dict = field(default_factory=dict)

    @classmethod
    def load(cls) -> "UiSettings":
        data = read_json(SETTINGS_PATH, default={}) or {}
        known = {key: data.get(key, getattr(cls, key, None))
                 for key in cls.__dataclass_fields__}
        known["recent_configs"] = list(known.get("recent_configs") or [])
        known["downstream_configs"] = dict(known.get("downstream_configs") or {})
        settings = cls(**known)
        settings.fill_defaults()
        return settings

    def save(self) -> None:
        try:
            atomic_write_json(SETTINGS_PATH, asdict(self))
        except OSError:
            pass          # settings are a convenience, never worth an error

    def fill_defaults(self) -> None:
        """Guess Bryan's location, so the first launch usually just works."""
        if not self.bryan_main:
            candidate = BRYAN_ROOT / "Main.py"
            if candidate.is_file():
                self.bryan_main = str(candidate)
        if not self.bryan_python:
            self.bryan_python = str(_guess_interpreter())

    def remember(self, config_path) -> None:
        text = str(Path(config_path).resolve())
        self.recent_configs = [text] + [
            entry for entry in self.recent_configs if entry != text
        ][:MAX_RECENT - 1]
        self.save()

    @staticmethod
    def _project_key(config_path) -> str:
        return str(Path(config_path).resolve())

    def downstream_for(self, config_path) -> dict:
        """The Downstream page's remembered paths for this project, or blanks."""
        stored = self.downstream_configs.get(self._project_key(config_path)) or {}
        return {"config": str(stored.get("config", "")),
                "model": str(stored.get("model", ""))}

    def remember_downstream(self, config_path, config="", model="") -> None:
        """Keep the two paths against this project. Blanks are kept as blanks rather
        than merged over, so clearing a field in the page clears it here too."""
        key = self._project_key(config_path)
        if not config and not model:
            self.downstream_configs.pop(key, None)
        else:
            self.downstream_configs[key] = {"config": str(config or ""),
                                            "model": str(model or "")}
        self.save()

    def problems(self) -> list:
        """Why a launch would fail, before the user tries it."""
        issues = []
        if not self.bryan_main:
            issues.append("Bryan's Main.py has not been set.")
        elif not Path(self.bryan_main).is_file():
            issues.append(f"Main.py not found: {self.bryan_main}")
        if not self.bryan_python:
            issues.append("The Python interpreter for Bryan has not been set.")
        elif not Path(self.bryan_python).exists():
            issues.append(f"Interpreter not found: {self.bryan_python}")
        return issues


def _guess_interpreter() -> Path:
    """Prefer a venv next to Bryan over whatever is running the UI.

    The existing .bat files use ``%PROJECT_ROOT%\\env\\python.exe``, so look
    there first - the UI's own interpreter almost certainly lacks scipy.
    """
    names = ("env/Scripts/python.exe", "env/bin/python",
             ".venv/Scripts/python.exe", ".venv/bin/python")
    for name in names:
        candidate = BRYAN_ROOT / name
        if candidate.exists():
            return candidate
    return Path(sys.executable)
