"""Has this row already run, and are its results still current?

Overwriting is destructive and irreversible: Bryan rewrites the results files
and ``shutil.rmtree``s the URBS working sub-folder on entry
(UrbsModel.__init__'s rmtree). There is no resume. So the UI works out, from the
filesystem, what state every row is in before anything is launched.

    NOT_RUN     no results database exists
    UP_TO_DATE  results exist and are newer than every input the row reads
    STALE       results exist but an input is newer  <- the one that matters
    INCOMPLETE  results exist but are truncated (a test_runs run, or a crash)
    NEEDS_PRIOR the inversion: Run models = no, and there is nothing to analyse

STALE is invisible today. Re-routing under an edited .sq, or re-running after
the storm config changed, leaves results that look perfectly fine.

**The Run models = no inversion.** A row with ``Run models`` = no and
``Analyse results`` = yes is a legitimate re-analysis of existing results - the
comment at lib/ReservoirRouting.py:790-792 says so. For those rows prior
results are a *requirement*, and their absence is an error, not 'not yet run'.

The run logs are deliberately NOT used to decide this. ``RunLog`` records
``Simulation`` = ``Output file`` (lib/RunLog.py:36) and no Duration, so in
CLD_RFSL_mc_sims_01.xlsx twenty-four rows share one name and the evidence is
ambiguous. The filesystem is the authority; the run log is provenance only
(``last_run_from_logs``).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .columns import normalise_method, path_columns_for, runs_models
from .paths import resolve_value, stat_or_none
from .outputs import find_primary, outputs_for

NOT_RUN = "not run"
UP_TO_DATE = "up to date"
STALE = "stale"
INCOMPLETE = "incomplete"
NEEDS_PRIOR = "needs prior results"
UNKNOWN = "unknown"

# Order for 'select everything that is not up to date'.
RERUN_WORTHY = (NOT_RUN, STALE, INCOMPLETE)


@dataclass(frozen=True)
class RowCompletion:
    state: str
    primary: Path | None = None
    newest_input: Path | None = None
    newest_input_mtime: float = 0.0
    result_mtime: float = 0.0
    detail: str = ""
    row_count: int | None = None
    expected_rows: int | None = None

    @property
    def has_results(self) -> bool:
        return self.primary is not None

    @property
    def needs_confirm_to_overwrite(self) -> bool:
        """Only UP_TO_DATE results are worth stopping the user over.

        Re-running a stale or incomplete row is the point, so it needs no
        extra confirmation - though the dialog still says what changed.
        """
        return self.state == UP_TO_DATE


def input_paths(row, config) -> list[Path]:
    """Every file this row reads, resolved absolute.

    The row's own path columns, plus the three global config files from
    sims_config.json - a changed climate or storm config makes every result
    in the project stale, and that is worth knowing.
    """
    method = normalise_method(row.get("Method"))
    found = []
    for column in path_columns_for(method):
        if column not in getattr(row, "index", ()):
            continue
        resolved = resolve_value(config.project_folder, row.get(column))
        if resolved is not None:
            found.append(resolved)
    found.extend(config.filepaths.values())
    return found


def assess(row, config, *, check_truncation: bool = False,
           mc_expected_rows: int | None = None) -> RowCompletion:
    """What state this row's results are in.

    ``check_truncation`` opens the results database to count its rows, so it is
    off by default - it is worth doing for the rows actually selected, not for
    all ninety-two on every render.
    """
    primary = find_primary(row, config.project_folder)
    analysis_only = not runs_models(row)

    if primary is None:
        outputs = outputs_for(row, config.project_folder)
        if analysis_only:
            return RowCompletion(
                NEEDS_PRIOR,
                detail=(
                    "'Run models' is no, so this row only re-analyses results "
                    "that already exist - and none were found. Expected one of: "
                    + ", ".join(str(p) for p in outputs.primary)
                ) if outputs.primary else outputs.note,
            )
        return RowCompletion(NOT_RUN, detail=outputs.note)

    result_stat = stat_or_none(primary)
    result_mtime = result_stat.st_mtime if result_stat else 0.0

    newest_path, newest_mtime = None, 0.0
    for candidate in input_paths(row, config):
        stat = stat_or_none(candidate)
        if stat and stat.st_mtime > newest_mtime:
            newest_path, newest_mtime = candidate, stat.st_mtime

    if check_truncation:
        counted, expected = _truncation(primary, mc_expected_rows)
        if counted is not None and expected and counted < expected:
            return RowCompletion(
                INCOMPLETE, primary, newest_path, newest_mtime, result_mtime,
                detail=(f"the results database holds {counted} rows, but the "
                        f"Monte Carlo config asks for {expected}. A run stopped "
                        f"early by 'test_runs', or one that crashed. Bryan's "
                        f"reservoir routing rejects an mcdf like this as input."),
                row_count=counted, expected_rows=expected,
            )

    if newest_mtime > result_mtime:
        return RowCompletion(
            STALE, primary, newest_path, newest_mtime, result_mtime,
            detail=f"{newest_path} is newer than the results",
        )

    return RowCompletion(UP_TO_DATE, primary, newest_path, newest_mtime,
                         result_mtime)


def _truncation(primary: Path, expected: int | None):
    """(rows in the database, rows expected). Either may be None."""
    if expected is None:
        return None, None
    try:
        if primary.suffix == ".parquet":
            return None, expected      # counting needs pyarrow; not worth it
        with primary.open("r", encoding="utf-8", errors="replace", newline="") as stream:
            count = sum(1 for _ in stream) - 1   # header
        return max(count, 0), expected
    except OSError:
        return None, expected


def expected_realisations(row, config) -> int | None:
    """main divisions x sub divisions from the row's Monte Carlo config.

    lib/ReservoirRouting.py reads these from ``scheme_config`` first, falling
    back to the top level (change_log, 31 July 2026) - a config written just
    for a routing run has them at the top. Mirror that lookup order.
    """
    import json

    config_file = resolve_value(config.project_folder, row.get("Config file"))
    if config_file is None or not config_file.is_file():
        return None
    try:
        data = json.loads(config_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

    scheme = data.get("scheme_config", {})
    main = scheme.get("number_of_main_divisions", data.get("number_of_main_divisions"))
    sub = scheme.get("number_of_sub_divisions", data.get("number_of_sub_divisions"))
    try:
        return int(main) * int(sub)
    except (TypeError, ValueError):
        return None


def assess_frame(frame, config, *, rows=None, check_truncation=False) -> dict:
    """Completion state for every row (or just ``rows``), keyed by frame index."""
    indices = list(frame.index) if rows is None else list(rows)
    out = {}
    for index in indices:
        row = frame.loc[index]
        expected = expected_realisations(row, config) if check_truncation else None
        out[index] = assess(row, config, check_truncation=check_truncation,
                            mc_expected_rows=expected)
    return out


def last_run_from_logs(config, extra_logs=()) -> dict:
    """When each ``Output file`` last completed, from the run logs.

    Provenance only - see the module docstring. Returns
    ``{Output file: (end time, status, computer)}``, last entry winning.
    """
    history: dict = {}
    candidates = [config.master_run_log, *extra_logs]
    run_root = config.run_root
    if run_root.is_dir():
        candidates.extend(sorted(run_root.glob("*/chunk_*_log.csv")))

    for path in candidates:
        try:
            with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as stream:
                for entry in csv.DictReader(stream):
                    name = (entry.get("Simulation") or "").strip()
                    if name:
                        history[name] = (
                            entry.get("End time", ""),
                            entry.get("Status", ""),
                            entry.get("Computer", ""),
                        )
        except (OSError, csv.Error):
            continue
    return history
