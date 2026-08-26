"""Where each method writes its results.

Taken from the code, not the manual, with the line each path came from. This is
what lets the UI say whether a row has already run without launching anything.

    monte carlo         <Output file>__mcdf.csv                MonteCarloSimulator
                        <Output file>_<type>.csv               MonteCarloSimulator
                        <Output file>_volume.csv               Simulator (volumes)
    ensemble            <Output file>.csv                      EnsembleSimulator
                        plots/ and csv/ beside it              lib/EnbAnalysis.py
    reservoir routing   <Results folder>/<Output file>__mcdf<suffix>.csv   (MC input)
                        <Results folder>/<Output file><suffix>.csv        (ensemble)
                        <Results folder>/<Output file>__<type>_quantiles<suffix>.csv
                        <Results folder>/<Output file>__inflow_volumes<suffix>.csv
                        <Hydrographs folder>/<Output file>_<series><suffix>.csv
                                                        lib/ReservoirRouting.py

``lib/ReservoirRouting.py:838-856`` (``_ensure_mcdf_loaded``) already probes the
reservoir-routing paths to find a previously-routed database, so this module
mirrors an existing convention rather than inventing one - the legacy candidate
below included.

**The pre-26-August-2026 Monte Carlo database carried no suffix.**
``_output_base`` tagged it ``__mcdf`` whatever ``Output suffix`` said, so every
suffix variant of one ``Output file`` overwrote the *same* file while
everything that identified the variant - the quantile tables, the volume table,
the hydrographs, the log - took the suffix. In
``TFD_SimsList_LongList_01.xlsx`` that is six ``Output file`` values times six
suffixes (FR-4B..FR-4J) over one results folder: thirty-six rows, six mcdf
files, and two hundred-odd suffixed quantile tables. Bryan now suffixes it too,
but those six files are still on disk and ``_ensure_mcdf_loaded`` will still
fall back to one. So the unsuffixed name stays a ``databases`` candidate - it
is readable - and is listed as ``shared``, never ``primary``: it proves that
*something* ran over this ``Output file``, not that this row did.

Bryan always *writes* csv - ``MCScheme.store_simulations`` is hard-coded to
``__mcdf.csv``. Parquet mcdfs in the wild are hand-converted, and
``ReservoirRouting._read_indexed`` (:85) falls back between the two extensions,
so the probes here accept either.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .columns import (ENSEMBLE, MONTE_CARLO, RESERVOIR_ROUTING,
                      analyses_results, analyses_volumes, normalise_method,
                      runs_models, stores_hydrographs)
from .paths import cell_text, resolve_value

RESULT_TYPES = ("inflow", "level", "outflow")
TABLE_SUFFIXES = (".csv", ".parquet")


@dataclass(frozen=True)
class OutputSet:
    """The files a row would write, split by how much they prove."""

    primary: tuple = ()      # this row's own results - absence means 'not run'
    secondary: tuple = ()    # volumes, plots - nice to have
    folders: tuple = ()      # directories that must exist for the row to run
    databases: tuple = ()    # the results database, in the order Bryan probes
                             # for it - what an analysis-only row reads back
    shared: tuple = ()       # written by this row but not owned by it: the
                             # unsuffixed mcdf, which every suffix variant of
                             # one Output file overwrites
    note: str = ""

    @property
    def all_paths(self) -> tuple:
        return tuple(self.primary) + tuple(self.secondary) + tuple(self.shared)


def _first_existing(base: Path, suffixes=TABLE_SUFFIXES) -> Path | None:
    for suffix in suffixes:
        candidate = base.with_suffix(suffix) if base.suffix else Path(str(base) + suffix)
        if candidate.is_file():
            return candidate
    return None


def existing_variant(base_without_extension: Path) -> Path | None:
    """The .csv or .parquet form of ``base``, whichever is on disk."""
    for suffix in TABLE_SUFFIXES:
        candidate = Path(str(base_without_extension) + suffix)
        if candidate.is_file():
            return candidate
    return None


def output_suffix(row) -> str:
    """lib/ReservoirRouting.py:737 - '_' + Output suffix, or ''."""
    suffix = cell_text(row.get("Output suffix"))
    return f"_{suffix}" if suffix else ""


def outputs_for(row, project_folder) -> OutputSet:
    """Every result path this row would write, resolved absolute."""
    method = normalise_method(row.get("Method"))
    output_file = cell_text(row.get("Output file"))
    if not output_file:
        return OutputSet(note="the row has no 'Output file', so its outputs "
                              "cannot be located")

    if method == RESERVOIR_ROUTING:
        return _reservoir_outputs(row, project_folder, output_file)
    if method == ENSEMBLE:
        return _ensemble_outputs(row, project_folder, output_file)
    if method == MONTE_CARLO:
        return _monte_carlo_outputs(row, project_folder, output_file)
    return OutputSet(note=f"unrecognised Method {row.get('Method')!r}")


def _monte_carlo_outputs(row, project_folder, output_file) -> OutputSet:
    base = resolve_value(project_folder, output_file)
    if base is None:
        return OutputSet()
    secondary = [Path(f"{base}_{kind}.csv") for kind in RESULT_TYPES]
    secondary.append(Path(f"{base}_volume.csv"))
    return OutputSet(
        primary=(Path(f"{base}__mcdf.csv"),),
        secondary=tuple(secondary),
        folders=(base.parent,),
        databases=(Path(f"{base}__mcdf.csv"),),
    )


def _ensemble_outputs(row, project_folder, output_file) -> OutputSet:
    base = resolve_value(project_folder, output_file)
    if base is None:
        return OutputSet()
    return OutputSet(
        primary=(Path(f"{base}.csv"),),
        secondary=(base.parent / "plots", base.parent / "csv"),
        folders=(base.parent,),
        databases=(Path(f"{base}.csv"),),
    )


def _reservoir_outputs(row, project_folder, output_file) -> OutputSet:
    results_folder = resolve_value(project_folder, row.get("Results folder"))
    if results_folder is None:
        return OutputSet(note="the row has no 'Results folder'")

    basename = Path(cell_text(output_file)).name
    suffix = output_suffix(row)

    # Which of the two the row writes depends on the scheme of its INPUT
    # database, which ReservoirRouting sniffs from that file's columns
    # (_detect_scheme). The UI does not open the input, so it accepts either -
    # the same candidates _ensure_mcdf_loaded probes, in its order.
    mc_database = results_folder / f"{basename}__mcdf{suffix}.csv"
    enb_database = results_folder / f"{basename}{suffix}.csv"
    legacy_database = results_folder / f"{basename}__mcdf.csv"
    databases = (mc_database, enb_database)
    if suffix:
        databases += (legacy_database,)

    # What the row leaves behind under its OWN name. The quantile tables are
    # the usual evidence; a row that routes without analysing leaves only the
    # hydrographs, and a row that does neither leaves nothing to find.
    quantiles = tuple(
        results_folder / f"{basename}__{kind}_quantiles{suffix}.csv"
        for kind in RESULT_TYPES
    ) if analyses_results(row) else ()

    hydrographs_folder = resolve_value(project_folder, row.get("Hydrographs folder"))
    hydrographs = ()
    if hydrographs_folder and runs_models(row) and stores_hydrographs(row):
        hydrographs = tuple(
            hydrographs_folder / f"{basename}_{series}{suffix}.csv"
            for series in ("outflows", "levels", "volumes")
        )

    secondary = ()
    if analyses_volumes(row):
        # One table of every event's volumes. The per-duration quantile files
        # beside it are named for the durations in the row's Config file, which
        # the UI does not open, so they are left out rather than guessed at.
        secondary = (results_folder / f"{basename}__inflow_volumes{suffix}.csv",)

    folders = (results_folder,) + ((hydrographs_folder,) if hydrographs_folder else ())

    # The legacy mcdf is a sibling's as much as it is this row's, so it never
    # counts as evidence - only as something an analysis could still read.
    primary = (mc_database, enb_database) + quantiles + hydrographs
    shared = (legacy_database,) if suffix else ()
    note = ("monte carlo input writes <name>__mcdf<suffix>.csv; ensemble input "
            "writes <name><suffix>.csv - either counts as having run")

    return OutputSet(
        primary=primary,
        secondary=secondary,
        folders=folders,
        databases=databases,
        shared=shared,
        note=note,
    )


def _first_on_disk(candidates, *, parquet_too: bool = False) -> Path | None:
    for candidate in candidates:
        if candidate.is_file():
            return candidate
        # accept a hand-converted parquet, as ReservoirRouting._read_indexed does
        if parquet_too:
            alternative = candidate.with_suffix(".parquet")
            if alternative.is_file():
                return alternative
    return None


def find_primary(row, project_folder) -> Path | None:
    """The output that proves THIS row has run, if any.

    Not simply the results database: a reservoir routing row with an
    ``Output suffix`` shares its mcdf with every other suffix over the same
    ``Output file``, so that file proves only that *something* ran.
    """
    outputs = outputs_for(row, project_folder)
    return _first_on_disk(outputs.primary, parquet_too=True)


def find_database(row, project_folder) -> Path | None:
    """The results database on disk, shared or not.

    What ``_ensure_mcdf_loaded`` will find for an analysis-only row, and the
    file worth counting rows in - a quantile table is ten rows by design.
    """
    outputs = outputs_for(row, project_folder)
    candidates = outputs.databases or outputs.primary
    return _first_on_disk(candidates, parquet_too=True)


def log_path_for(row, project_folder) -> Path | None:
    """The per-simulation log this row writes.

    The ``Log file`` column when it has one; for reservoir routing it may be
    blank, in which case lib/ReservoirRouting.py:440 puts it in the results
    folder as ``<Output file><suffix>_log.txt``.
    """
    declared = cell_text(row.get("Log file"))
    if declared:
        return resolve_value(project_folder, declared)

    if normalise_method(row.get("Method")) != RESERVOIR_ROUTING:
        return None
    results_folder = resolve_value(project_folder, row.get("Results folder"))
    output_file = cell_text(row.get("Output file"))
    if results_folder is None or not output_file:
        return None
    basename = Path(output_file).name
    return results_folder / f"{basename}{output_suffix(row)}_log.txt"


def collision_key(row) -> tuple:
    """What makes two rows write over each other.

    ``Output file`` alone is the URBS working-folder key
    (``Simulator.initialise_model`` -> ``UrbsModel.__init__``'s rmtree).
    Together with the suffix and results folder it
    is also the results-file key. In CLD_RFSL_mc_sims_01.xlsx twenty-four rows
    share one of these, differing only by Duration.
    """
    return (
        cell_text(row.get("Output file")),
        cell_text(row.get("Output suffix")),
        cell_text(row.get("Results folder")),
    )


def working_folder_key(row) -> str:
    """The URBS/RORB sub-folder name, from ``Simulator.initialise_model``.

    Two simulations sharing this must never run at the same time: the model
    wrapper rmtree's the folder on entry.
    """
    return Path(cell_text(row.get("Output file"))).name


# ---------------------------------------------------------------------------
# Quantile files - what the results viewer plots
# ---------------------------------------------------------------------------

# The volume analyses tag the FILE 'inflowVol24h' and the COLUMN inside it
# 'Vol24h' - Simulator.analyse_volumes and ReservoirRouting._tpt_result both
# say so, and util/MaxQuantiles.py reads both, so neither is free to change.
VOLUME_RE = re.compile(r"^(?P<kind>inflow|outflow)Vol(?P<window>\d+(?:[._]\d+)?)h$",
                       re.IGNORECASE)


@dataclass(frozen=True)
class QuantileFile:
    """One standard-AEP quantile table on disk.

    ``key`` is how the file is tagged; ``column`` is the column to read out of
    it. They differ only for the volumes, where a file called ``inflowVol24h``
    holds a column called ``Vol24h``.
    """

    key: str
    column: str
    path: Path

    @property
    def is_volume(self) -> bool:
        return VOLUME_RE.match(self.key) is not None

    @property
    def label(self) -> str:
        match = VOLUME_RE.match(self.key)
        if not match:
            return self.key
        window = match.group("window").replace("_", ".")
        return f"{match.group('kind').lower()} volume, {window} h window"


def quantile_files(row, project_folder) -> dict:
    """Every quantile table this row has actually written, by key.

    Monte carlo and reservoir routing name these differently but write the same
    three columns - ``aep (1 in x)``, ``probability``, ``<column>`` - from
    lib/MCScheme.py:353 and lib/ReservoirRouting.py:310 respectively. Only
    files that exist are returned: the point of this is what can be plotted.
    """
    method = normalise_method(row.get("Method"))
    if method == MONTE_CARLO:
        return _monte_carlo_quantiles(row, project_folder)
    if method == RESERVOIR_ROUTING:
        return _reservoir_quantiles(row, project_folder)
    return {}          # the ensemble method writes its own critical-duration
                       # analysis instead - see lib/EnbAnalysis.py


def _found(found: dict, key: str, path: Path) -> None:
    if not path.is_file():
        return
    match = VOLUME_RE.match(key)
    found[key] = QuantileFile(key=key,
                              column=f"Vol{match.group('window')}h" if match else key,
                              path=path)


def _monte_carlo_quantiles(row, project_folder) -> dict:
    base = resolve_value(project_folder, cell_text(row.get("Output file")))
    if base is None:
        return {}
    found: dict = {}
    for kind in RESULT_TYPES:
        _found(found, kind, Path(f"{base}_{kind}.csv"))
    for path in _safe_glob(base.parent, f"{base.name}_*Vol*h.csv"):
        _found(found, path.stem[len(base.name) + 1:], path)
    return found


def _reservoir_quantiles(row, project_folder) -> dict:
    results_folder = resolve_value(project_folder, row.get("Results folder"))
    output_file = cell_text(row.get("Output file"))
    if results_folder is None or not output_file:
        return {}
    basename = Path(output_file).name
    suffix = output_suffix(row)
    found: dict = {}
    for kind in RESULT_TYPES:
        _found(found, kind,
               results_folder / f"{basename}__{kind}_quantiles{suffix}.csv")
    tail = f"_quantiles{suffix}.csv"
    for path in _safe_glob(results_folder, f"{basename}__*Vol*h{tail}"):
        _found(found, path.name[len(basename) + 2:-len(tail)], path)
    return found


def _safe_glob(folder: Path, pattern: str) -> list:
    """Globbing a folder that may not exist, or may not be readable."""
    try:
        return sorted(folder.glob(pattern))
    except OSError:
        return []
