#!/usr/bin/env bash
# Juniper Creek Dam example - run_urbs. POSIX companion to run_urbs.bat.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

BRYAN_ROOT="${BRYAN_ROOT:-..}"
VENV_PY="${VENV_PY:-$(command -v python3 || command -v python)}"

exec "$VENV_PY" "$BRYAN_ROOT/Main.py" "sims_config_urbs.json"
