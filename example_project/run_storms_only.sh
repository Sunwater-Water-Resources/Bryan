#!/usr/bin/env bash
# Juniper Creek Dam example - run_storms_only. POSIX companion to run_storms_only.bat.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

BRYAN_ROOT="${BRYAN_ROOT:-..}"
VENV_PY="${VENV_PY:-$(command -v python3 || command -v python)}"

exec "$VENV_PY" "$BRYAN_ROOT/Main.py" "sims_config_storms.json"
