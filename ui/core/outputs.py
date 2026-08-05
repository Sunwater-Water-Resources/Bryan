"""Where each method writes its results.

Taken from the code, not the manual, with the line each path came from. This is
what lets the UI say whether a row has already run without launching anything.

    monte carlo         <Output file>__mcdf.csv                MonteCarloSimulator
                        <Output file>_<type>.csv               MonteCarloSimulator
                        <Output file>_volume.csv               Simulator (volumes)
    ensemble            <Output file>.csv                      EnsembleSimulator
                        plots/ and csv/ beside it              lib/EnbAnalysis.py
    reservoir routing   <Results folder>/<Output file>__mcdf.csv          (MC input)
                        <Results folder>/<Output file><suffix>.csv        (ensemble)
                        <Results folder>/<Output file>__<type>_quantiles<suffix>.csv
                                                        lib/ReservoirRouting.py:770,771,1019

``lib/ReservoirRouting.py:787-803`` (``_ensure_mcdf_loaded``) already probes the
two reservoir-routing paths to find a previously-routed database, so this module
mirrors an existing convention rather than inventing one.

Bryan always *writes* csv - ``MCScheme.store_simulations`` is hard-coded to
``__mcdf.csv``. Parquet mcdfs in the wild are hand-converted, and
``ReservoirRouting._read_indexed`` (:85) falls back between the two extensions,
so the probes here accept either.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .columns import ENSEMBLE, MONTE_CARLO, RESERVOIR_ROUTING, normalise_method
from .paths import cell_text, resolve_value

RESULT_TYPES = ("inflow", "level", "outflow")
TABLE_SUFFIXES = (".csv", ".parquet")


@dataclass(frozen=True)
class OutputSet:
    """The files a row would write, split by how much they prove."""

    primary: tuple = ()      # the results database - its absence means 'not run'
    secondary: tuple = ()    # quantiles, plots, volumes - nice to have
    folders: tuple = ()      # directories that must exist for the row to run
    note: str = ""

    @property
    def all_paths(self) -> tuple:
        return tuple(self.primary) + tuple(self.secondary)


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
    )


def _ensemble_outputs(row, project_folder, output_file) -> OutputSet:
    base = resolve_value(project_folder, output_file)
    if base is None:
        return OutputSet()
    return OutputSet(
        primary=(Path(f"{base}.csv"),),
        secondary=(base.parent / "plots", base.parent / "csv"),
        folders=(base.parent,),
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
    # the same two candidates _ensure_mcdf_loaded probes at :794.
    mc_base = results_folder / f"{basename}__mcdf"
    enb_base = results_folder / f"{basename}{suffix}"

    secondary = tuple(
        results_folder / f"{basename}__{kind}_quantiles{suffix}.csv"
        for kind in RESULT_TYPES
    )
    hydrographs = resolve_value(project_folder, row.get("Hydrographs folder"))
    folders = (results_folder,) + ((hydrographs,) if hydrographs else ())

    return OutputSet(
        primary=(Path(f"{mc_base}.csv"), Path(f"{enb_base}.csv")),
        secondary=secondary,
        folders=folders,
        note="monte carlo input writes __mcdf.csv; ensemble input writes "
             "<name><suffix>.csv - either counts as having run",
    )


def find_primary(row, project_folder) -> Path | None:
    """The results database this row has already written, if any."""
    outputs = outputs_for(row, project_folder)
    for candidate in outputs.primary:
        if candidate.is_file():
            return candidate
        # accept a hand-converted parquet, as ReservoirRouting._read_indexed does
        alternative = candidate.with_suffix(".parquet")
        if alternative.is_file():
            return alternative
    return None


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
