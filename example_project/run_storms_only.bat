TITLE Juniper Creek Dam example - run_storms_only
cd /D "%~dp0"

:: Point these at your Bryan clone and its Python environment
set "BRYAN_ROOT=..\.."
set "VENV_PY=%BRYAN_ROOT%\env\python.exe"

"%VENV_PY%" "%BRYAN_ROOT%\Main.py" "sims_config_storms.json"

pause
