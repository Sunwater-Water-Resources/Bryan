"""What Bryan reads out of a sims-list row, and which columns hold paths.

This is the machine-readable version of ``Manual/SubDocs/sim_list.md``, but it
is written from the CODE, because the two disagree in places (see
``REPLICATE_FILE_ALIASES``). Entries name the function they came from rather
than a line number: lib/Simulator.py is under active development on other
branches, so its line numbers do not survive a merge.

The distinction that matters throughout is *when* a column is read:

- ``always`` - read for every included row of that method. A missing column is
  a KeyError that costs the whole simulation.
- ``when_running`` - read only when ``Run models`` is yes/storms only.
  ``Simulator.__init__`` returns at its ``if not self.do_runs`` guard before
  reading these, so an analysis-only row genuinely does not need them.

Pre-flight uses this to fail in a millisecond instead of an hour in.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .paths import cell_text, is_blank

MONTE_CARLO = "monte carlo"
ENSEMBLE = "ensemble"
RESERVOIR_ROUTING = "reservoir routing"
METHODS = (MONTE_CARLO, ENSEMBLE, RESERVOIR_ROUTING)

# Simulator.__init__ reads 'Run models'. 'storms only' generates storms but
# skips the hydrologic model; it does not exist on every branch, and treating
# an unknown value as 'runs' is the conservative way round.
RUN_MODES_THAT_RUN = ("yes", "storms only")

# 'Analyse volumes' values that switch the volume analysis OFF. Everything else
# switches it on: Simulator.__init__ and
# ReservoirRoutingSimulator._get_volume_setting accept different spellings
# ('yes'/'inflow'/'outflow'), and an unknown value is the conservative way
# round here too - it only ever adds an input to check.
VOLUME_MODES_THAT_SKIP = ("", "no", "none", "false")

# Read for every row of every method, before any method dispatch:
#   Main.py:77,81,82 and lib/RunLog.py:36,42,43,49
UNIVERSAL_ALWAYS = (
    "Include",
    "Method",
    "Output file",
    "Run models",
    "Analyse results",
    "Config file",
)

# Simulator.__init__, unconditionally - ahead of its `if not self.do_runs`
# guard, so an analysis-only row still needs every one of these.
SIMULATOR_ALWAYS = (
    "Run models",
    "Method",
    "Config file",
    "Duration",
    "Output file",
    "Analyse results",
    "Store hydrographs",
    "Mop up files",
)

# Simulator.__init__, after the do_runs guard.
SIMULATOR_WHEN_RUNNING = (
    "Log file",             # opened inside the do_runs branch
    "Exclusions",
    "Replicate file",       # see REPLICATE_FILE_ALIASES
    "Replicates",
    "IL",
    "CL",
    "ADV",                  # lib/Lake.py, via LakeConditions(parameters)
    "Focal subcatchments",  # read with .get(), but a falsy value raises
)

# ReservoirRoutingSimulator.__init__, whatever Run models says.
RESERVOIR_ALWAYS = (
    "Output file",
    "Hydrographs folder",
    "Results folder",
    "Run models",
    "Analyse results",
)

# lib/ReservoirRouting.py - needed to do the actual routing.
RESERVOIR_WHEN_RUNNING = (
    "ELS file",
    "SQ file",
    "FSL",
    "Inflow",
)

# Guarded by `in sim_row.index` / pd.notna - absent is fine.
RESERVOIR_OPTIONAL = (
    "Log file",
    "Analyse volumes",
    "ADV",
    "ADV source",
    "Output suffix",
    "Store hydrographs",
    "Config file",
    "Lake config",
    "Input MCDF",
    "Input database",
)

# ReservoirRouting accepts either spelling for the input database.
INPUT_DATABASE_ALIASES = ("Input MCDF", "Input database")

# Simulator.__init__ reads 'Replicate file'; Manual/SubDocs/sim_list.md
# documents 'Replication file'. A list built from the manual raises KeyError on
# every monte carlo row. Accept both, warn on the documented-but-wrong one.
REPLICATE_FILE_ALIASES = ("Replicate file", "Replication file")

# Simulator.__init__'s climate block looks GWL up with `in parameters.keys()`;
# without it, Year and SSP are read unconditionally.
CLIMATE_GWL = ("GWL",)
CLIMATE_SSP = ("Year", "SSP")

# Columns whose value is a path, per method. Used by pre-flight (does it exist?)
# and by completion (is the result newer than its inputs?).
PATH_COLUMNS_COMMON = ("Config file", "Lake config", "Focal subcatchments",
                       "Replicate file", "Replication file", "TP weights")
PATH_COLUMNS_RESERVOIR = ("Input MCDF", "Input database", "Inflow",
                          "ELS file", "SQ file")

# Columns Bryan never reads. Shown in the table, flagged as informational, and
# copied into run copies untouched. Common in real lists: Design, Catchment,
# Rev ID, Run ID suffix, Basename, Model ID, Comment, EBF on, D50 on, Group.
KNOWN_UNUSED = ("Comment", "Design", "Catchment", "Rev ID", "Run ID suffix",
                "Model ID", "sims_config (not used)", "Group", "Source row")


@dataclass(frozen=True)
class ColumnRequirement:
    """One column Bryan will read for a given row."""

    name: str
    when: str            # 'always' | 'when_running'
    is_path: bool = False
    aliases: tuple[str, ...] = field(default_factory=tuple)

    def present_in(self, columns) -> str | None:
        """The name this requirement is satisfied by, or None."""
        for candidate in (self.name, *self.aliases):
            if candidate in columns:
                return candidate
        return None


def normalise_method(value) -> str:
    """A Method cell as one of METHODS, or '' when unrecognisable.

    Lenient on purpose, so the UI can still classify and display a row whose
    Method reads 'Ensemble'. Whether Bryan would *accept* it is a different
    question - see ``method_is_exact``.
    """
    text = cell_text(value).lower()
    return text if text in METHODS else ""


def method_is_exact(value) -> bool:
    """Whether Bryan would dispatch on this Method cell.

    Main.py:88-95 compares the raw cell against the lowercase literals and
    raises 'Modelling method X not recognised' otherwise, so 'Ensemble' with a
    capital E fails the whole simulation. Case is not forgiven.
    """
    return cell_text(value) in METHODS


def runs_models(row) -> bool:
    """Whether Bryan will execute a model for this row.

    False means an analysis-only row: it needs FEWER inputs but it needs its
    prior RESULTS to exist - see core/completion.py.
    """
    return cell_text(row.get("Run models")).lower() in RUN_MODES_THAT_RUN


def analyses_volumes(row) -> bool:
    """Whether the row switches on the flood volume analysis.

    It matters for reservoir routing: lib/ReservoirRouting.py measures the
    volumes off the inflow hydrographs, so a row that sets this needs its
    'Inflow' file even when Run models is no.
    """
    return cell_text(row.get("Analyse volumes")).lower() not in VOLUME_MODES_THAT_SKIP


def included(row) -> bool:
    """Main.py:64 - an exact, case-sensitive compare, so 'Yes' silently skips."""
    return row.get("Include") == "yes"


def requirements_for(method: str, running: bool,
                     volumes: bool = False) -> list[ColumnRequirement]:
    """Every column Bryan will read for a row of this method and run mode.

    When ``running`` is False the 'when_running' columns are left out entirely:
    ``Simulator.__init__`` returns at its do_runs guard before reading them, and
    a reservoir-routing row that only re-analyses never opens its .sq or its
    inflows. Reporting them as required would flag good analysis-only rows.

    ``volumes`` is the one exception to that: a reservoir-routing row analysing
    the inflow volumes opens the hydrographs whatever Run models says, because
    that is where the volumes are measured. Pass ``analyses_volumes(row)``.
    """
    out: dict[str, ColumnRequirement] = {}

    def add(name, when, is_path=False, aliases=()):
        if when == "when_running" and not running:
            return
        # 'always' wins if a column appears under both.
        existing = out.get(name)
        if existing and existing.when == "always":
            return
        out[name] = ColumnRequirement(name, when, is_path, tuple(aliases))

    for name in UNIVERSAL_ALWAYS:
        add(name, "always", name in PATH_COLUMNS_COMMON)

    if method == RESERVOIR_ROUTING:
        for name in RESERVOIR_ALWAYS:
            add(name, "always")
        for name in RESERVOIR_WHEN_RUNNING:
            add(name, "when_running", name in PATH_COLUMNS_RESERVOIR)
        if volumes:
            # ReservoirRoutingSimulator._ensure_inflows_loaded: the volume
            # analysis reads the hydrographs even for an analysis-only row.
            add("Inflow", "always", True)
        if running:
            add("Input MCDF", "when_running", True, ("Input database",))
    else:
        for name in SIMULATOR_ALWAYS:
            add(name, "always", name in PATH_COLUMNS_COMMON)
        for name in SIMULATOR_WHEN_RUNNING:
            aliases = REPLICATE_FILE_ALIASES[1:] if name == "Replicate file" else ()
            add(name, "when_running", name in PATH_COLUMNS_COMMON, aliases)

    return [out[name] for name in sorted(out)]


def path_columns_for(method: str) -> tuple[str, ...]:
    """Columns holding a path for this method, whatever the run mode."""
    if method == RESERVOIR_ROUTING:
        return PATH_COLUMNS_RESERVOIR + ("Config file", "Lake config")
    return PATH_COLUMNS_COMMON


def climate_requirement(row) -> tuple[str, tuple[str, ...]]:
    """Which climate columns this row needs.

    Returns ('gwl'|'ssp', names). GWL wins if the column exists at all - even
    blank, which Bryan reads as 0 degrees.
    """
    if "GWL" in row.index:
        return "gwl", CLIMATE_GWL
    return "ssp", CLIMATE_SSP


def unused_columns(columns) -> list[str]:
    """Columns present in the sheet that no method of Bryan reads."""
    read = set(UNIVERSAL_ALWAYS) | set(SIMULATOR_ALWAYS) | set(SIMULATOR_WHEN_RUNNING)
    read |= set(RESERVOIR_ALWAYS) | set(RESERVOIR_WHEN_RUNNING) | set(RESERVOIR_OPTIONAL)
    read |= set(PATH_COLUMNS_COMMON) | set(PATH_COLUMNS_RESERVOIR)
    read |= set(CLIMATE_GWL) | set(CLIMATE_SSP) | set(REPLICATE_FILE_ALIASES)
    read |= {"Analyse volumes", "Analyse sub-bursts", "Pre-burst method",
             "Baseflow", "IL percentile", "Preburst percentile"}
    return [str(name) for name in columns if str(name) not in read]


__all__ = [
    "MONTE_CARLO", "ENSEMBLE", "RESERVOIR_ROUTING", "METHODS",
    "ColumnRequirement", "normalise_method", "runs_models", "included",
    "analyses_volumes",
    "requirements_for", "path_columns_for", "climate_requirement",
    "unused_columns", "is_blank", "cell_text",
    "INPUT_DATABASE_ALIASES", "REPLICATE_FILE_ALIASES", "KNOWN_UNUSED",
]
