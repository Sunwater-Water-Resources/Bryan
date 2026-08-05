"""Splitting a selection into chunks that can safely run at the same time.

The hard constraint comes from ``Simulator.initialise_model``, which derives
the URBS working sub-folder from the row's ``Output file``:

    run_id = os.path.split(self.outputfile)[1]      # Simulator.initialise_model
    ...
    if os.path.exists(self.output_folder): shutil.rmtree(self.output_folder)
    if os.path.exists(self.storms_folder): shutil.rmtree(self.storms_folder)
                                                    # UrbsModel.__init__

Two simulations sharing an ``Output file`` delete each other's working
directories mid-run. So rows sharing one are an indivisible atom, and the
packer physically cannot produce a plan that splits them.

Parallelism is off by default. See the note in ui/README.md: reservoir routing
takes seconds, and Monte Carlo already drives thousands of *serial* URBS runs,
so N processes mostly multiply peak memory and storm-file count rather than
saving wall clock.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import pandas as pd

from .columns import ENSEMBLE, MONTE_CARLO, RESERVOIR_ROUTING, normalise_method, runs_models
from .grouping import group_key
from .outputs import collision_key, working_folder_key
from .paths import cell_text

BLOCK = "block"
WARN = "warn"
INFO = "info"

# Rough relative cost of one simulation-hour of storm duration, by method.
# Monte Carlo runs thousands of model invocations; reservoir routing runs none.
METHOD_WEIGHT = {MONTE_CARLO: 60.0, ENSEMBLE: 6.0, RESERVOIR_ROUTING: 0.3}
DEFAULT_WEIGHT = 6.0
NOT_RUNNING_FACTOR = 0.05     # 'Run models' = no: analysis only
DEFAULT_DURATION_HOURS = 24.0


@dataclass(frozen=True)
class Chunk:
    index: int
    rows: tuple
    estimated_minutes: float

    @property
    def size(self) -> int:
        return len(self.rows)


@dataclass(frozen=True)
class Hazard:
    severity: str
    code: str
    message: str
    rows: tuple = ()
    fix_hint: str = ""

    @property
    def blocks(self) -> bool:
        return self.severity == BLOCK


@dataclass(frozen=True)
class RunPlan:
    chunks: tuple
    hazards: tuple
    log_file_overrides: dict = field(default_factory=dict)
    log_renames: tuple = ()

    @property
    def blocked(self) -> bool:
        return any(hazard.blocks for hazard in self.hazards)

    @property
    def rows(self) -> list:
        return [index for chunk in self.chunks for index in chunk.rows]

    @property
    def estimated_minutes(self) -> float:
        """Wall clock, i.e. the slowest chunk - they run at the same time."""
        return max((chunk.estimated_minutes for chunk in self.chunks), default=0.0)


def estimate_minutes(row, history=None) -> float:
    """How long this row is likely to take.

    A heuristic, overridden by measurement wherever the run logs have seen this
    ``Output file`` before - that is what makes the per-chunk estimate worth
    showing rather than decorative.
    """
    name = cell_text(row.get("Output file"))
    if history and name in history:
        return history[name]

    method = normalise_method(row.get("Method"))
    weight = METHOD_WEIGHT.get(method, DEFAULT_WEIGHT)
    try:
        duration = float(row.get("Duration"))
        if pd.isna(duration) or duration <= 0:
            duration = DEFAULT_DURATION_HOURS
    except (TypeError, ValueError):
        duration = DEFAULT_DURATION_HOURS
    minutes = duration * weight / 60.0
    if not runs_models(row):
        minutes *= NOT_RUNNING_FACTOR
    return max(minutes, 0.05)


def observed_minutes(config, extra_logs=()) -> dict:
    """Measured run times per ``Output file``, from the run logs.

    lib/RunLog.py writes 'Start time' and 'End time' as '%Y-%m-%d %H:%M', so
    the resolution is a minute - fine for balancing chunks.
    """
    fmt = "%Y-%m-%d %H:%M"
    out: dict = {}
    candidates = [config.master_run_log, *extra_logs]
    if config.run_root.is_dir():
        candidates.extend(sorted(config.run_root.glob("*/chunk_*_log.csv")))

    for path in candidates:
        try:
            with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as stream:
                for entry in csv.DictReader(stream):
                    name = (entry.get("Simulation") or "").strip()
                    if not name:
                        continue
                    try:
                        start = datetime.strptime(entry["Start time"].strip(), fmt)
                        end = datetime.strptime(entry["End time"].strip(), fmt)
                    except (KeyError, ValueError, AttributeError):
                        continue
                    minutes = (end - start).total_seconds() / 60.0
                    if minutes >= 0:
                        out[name] = max(minutes, 0.05)
        except (OSError, csv.Error):
            continue
    return out


def _atoms(frame, rows, keep_groups_together):
    """Indivisible units, in the order they appear in the sheet."""
    buckets: dict = {}
    order: list = []
    for index in rows:
        row = frame.loc[index]
        if keep_groups_together:
            key = ("group", group_key(row).key)
        else:
            key = ("output", working_folder_key(row))
        # An Output file is indivisible whatever the grouping says, so when
        # groups are off the key IS the output name; when they are on, a group
        # always contains whole output names anyway.
        if key not in buckets:
            buckets[key] = []
            order.append(key)
        buckets[key].append(index)
    return [tuple(buckets[key]) for key in order]


def plan_run(sims, selected_rows, n_chunks=1, *, keep_groups_together=True,
             history=None) -> RunPlan:
    """Split ``selected_rows`` into ``n_chunks`` that can run concurrently."""
    frame = sims.frame if hasattr(sims, "frame") else sims
    rows = [index for index in selected_rows if index in frame.index]
    if not rows:
        return RunPlan(chunks=(), hazards=(Hazard(
            BLOCK, "empty-selection", "No simulations are selected."),))

    n_chunks = max(1, int(n_chunks))
    atoms = _atoms(frame, rows, keep_groups_together)
    costs = {index: estimate_minutes(frame.loc[index], history) for index in rows}

    # Longest-processing-time-first: sort atoms by cost, drop each into the
    # lightest bin. Simple, and good enough for a handful of bins.
    weighted = sorted(
        ((sum(costs[index] for index in atom), atom) for atom in atoms),
        key=lambda item: -item[0],
    )
    bins = [[0.0, []] for _ in range(min(n_chunks, len(weighted)))]
    for cost, atom in weighted:
        target = min(bins, key=lambda entry: entry[0])
        target[0] += cost
        target[1].extend(atom)

    # Keep sheet order within a chunk so 'Simulation N of M' maps predictably.
    order = {index: position for position, index in enumerate(rows)}
    chunks = tuple(
        Chunk(number, tuple(sorted(members, key=order.__getitem__)), round(cost, 2))
        for number, (cost, members) in enumerate(
            ((entry[0], entry[1]) for entry in bins if entry[1]), start=1)
    )

    overrides, renames = _global_log_dedupe(frame, chunks)
    hazards = analyse_hazards(frame, chunks, n_chunks, overrides)
    return RunPlan(chunks=chunks, hazards=tuple(hazards),
                   log_file_overrides=overrides, log_renames=tuple(renames))


def _global_log_dedupe(frame, chunks):
    """Run Bryan's own dedupe once over every chunk, in plan order.

    Returns ({frame index: log path}, [(old, new), ...]).
    """
    ordered = [index for chunk in chunks for index in chunk.rows]
    if not ordered or "Log file" not in frame.columns:
        return {}, []

    from .bryan import resolve_duplicate_logs

    subset = frame.loc[ordered]
    before = subset["Log file"].tolist()
    resolved, _ = resolve_duplicate_logs(subset, capture=True)
    after = resolved["Log file"].tolist()

    overrides = {index: value for index, value in zip(ordered, after)}
    renames = [(old, new) for old, new in zip(before, after)
               if str(old) != str(new)]
    return overrides, renames


def analyse_hazards(frame, chunks, requested_chunks, log_overrides=None) -> list:
    """Everything about this plan worth stopping or warning over."""
    hazards: list = []
    chunk_of = {index: chunk.index
                for chunk in chunks for index in chunk.rows}
    rows = list(chunk_of)

    # 1. The hard one: an Output file split across chunks.
    by_output: dict = {}
    for index in rows:
        by_output.setdefault(working_folder_key(frame.loc[index]), set()).add(chunk_of[index])
    for name, chunk_ids in sorted(by_output.items()):
        if name and len(chunk_ids) > 1:
            hazards.append(Hazard(
                BLOCK, "output-split",
                f"The output name {name!r} is spread across {len(chunk_ids)} "
                f"chunks. Those simulations share a URBS working folder and "
                f"delete it on entry, so they would destroy each other.",
                rows=tuple(i for i in rows if working_folder_key(frame.loc[i]) == name),
                fix_hint="This is a bug in the planner if you see it - rows "
                         "sharing an output name are meant to be indivisible.",
            ))

    # 2. Rows that overwrite each other even run one after another.
    by_collision: dict = {}
    for index in rows:
        by_collision.setdefault(collision_key(frame.loc[index]), []).append(index)
    for key, members in sorted(by_collision.items(), key=lambda item: -len(item[1])):
        if len(members) > 1 and key[0]:
            hazards.append(Hazard(
                BLOCK, "output-collision",
                f"{len(members)} selected rows write to the same output name "
                f"{key[0]!r}. Only the last would survive.",
                rows=tuple(members),
                fix_hint="Give them distinct 'Output file' values, or run them "
                         "one at a time.",
            ))

    # 3. Per-simulation logs that still collide after the global dedupe.
    if log_overrides:
        seen: dict = {}
        for index, value in log_overrides.items():
            text = cell_text(value)
            if text:
                seen.setdefault(text.lower(), []).append(index)
        clashing = [members for members in seen.values() if len(members) > 1]
        if clashing:
            hazards.append(Hazard(
                BLOCK, "log-collision",
                f"{sum(len(m) for m in clashing)} rows would still share a log "
                f"file after deduplication. Each opens it with 'w', so all but "
                f"one run's log is lost.",
                rows=tuple(index for members in clashing for index in members),
            ))

    if len(chunks) > 1:
        cpus = os.cpu_count() or 1
        if len(chunks) > max(1, cpus - 1):
            hazards.append(Hazard(
                WARN, "parallel-oversubscribed",
                f"{len(chunks)} chunks on a machine with {cpus} CPUs. Each "
                f"Bryan process drives its models serially, so a chunk is "
                f"roughly one core.",
            ))

        monte_carlo = sum(
            1 for chunk in chunks
            if any(normalise_method(frame.loc[i].get("Method")) == MONTE_CARLO
                   for i in chunk.rows)
        )
        if monte_carlo > 1:
            hazards.append(Hazard(
                WARN, "memory",
                f"{monte_carlo} chunks run Monte Carlo simulations. Each "
                f"process holds a full realisation database and its "
                f"hydrographs, so peak memory scales with the chunk count.",
            ))

        mop_off = tuple(index for index in rows
                        if cell_text(frame.loc[index].get("Mop up files")).lower() == "no")
        if mop_off:
            hazards.append(Hazard(
                WARN, "mop-up-off",
                f"{len(mop_off)} row(s) have 'Mop up files' = no. Running "
                f"several at once multiplies the storm and result files left "
                f"on disk.",
                rows=mop_off,
            ))

        hazards.append(Hazard(
            INFO, "per-chunk-logs",
            f"Each chunk writes its own run log (chunk_NN_log.csv) because "
            f"lib/RunLog.py derives the name from the simulation list. The UI "
            f"reads all {len(chunks)} of them together.",
        ))

    hazards.extend(stringify_hazards(frame, chunks))

    if requested_chunks > len(chunks) and chunks:
        hazards.append(Hazard(
            INFO, "fewer-chunks",
            f"Asked for {requested_chunks} chunks but the selection only "
            f"splits into {len(chunks)} indivisible unit(s).",
        ))
    return hazards


def stringify_hazards(frame, chunks) -> list:
    """Warn when splitting the list would change ``str(Duration)``.

    An xlsx cell cannot carry the int/float distinction: openpyxl writes 120.0
    as "120", and pandas reads that back as int. So a column is float64 only
    while *some* value in it is fractional - and a chunk that leaves that row
    out silently becomes int64.

    That matters because Bryan stringifies the duration:

        URBSmodel.duration_string(120.0) -> '120.0'   (storm file .120.0h)
        URBSmodel.duration_string(120)   -> '120'     (storm file .120h)

    and ``get_simulation_period`` looks the result up in the URBS config's JSON
    keys, falling back to twice the storm duration when it misses (:287) with
    only a printed note.

    This is live in CLD_RFSL_mc_sims_01.xlsx: it has a 4.5 hour duration, so
    the whole column is float64, and any selection excluding that row would
    behave differently from running the master list. Reservoir-routing rows are
    unaffected - they take their durations from the stored results - so this is
    a warning, not a block.
    """
    if "Duration" not in getattr(frame, "columns", ()):
        return []
    column = frame["Duration"]
    if not pd.api.types.is_float_dtype(column):
        return []

    hazards = []
    for chunk in chunks:
        values = column.loc[list(chunk.rows)].dropna()
        if values.empty or not all(float(v).is_integer() for v in values):
            continue          # stays float, so str() is unchanged

        affected = tuple(
            index for index in chunk.rows
            if not pd.isna(column.loc[index])
        )
        if not affected:
            continue
        methods = {normalise_method(frame.loc[i].get("Method")) for i in affected}
        if methods <= {RESERVOIR_ROUTING}:
            continue          # Duration is not used by reservoir routing

        sample = column.loc[affected[0]]
        hazards.append(Hazard(
            WARN, "duration-dtype",
            f"chunk {chunk.index}: the master list's Duration column is "
            f"fractional somewhere (so pandas reads it as float), but this "
            f"chunk holds whole numbers only. Bryan would see "
            f"{str(int(sample))!r} where the full list gives "
            f"{str(float(sample))!r}, which changes the storm filenames and "
            f"the URBS simulation-period lookup.",
            rows=affected,
            fix_hint="Include the fractional-duration row in the same "
                     "selection, or run the master list, if you need the run "
                     "to match an earlier one exactly.",
        ))
    return hazards


def describe(plan: RunPlan, frame) -> str:
    """A short human summary for the confirm dialog."""
    lines = []
    for chunk in plan.chunks:
        names = [cell_text(frame.loc[i].get("Output file")) or f"row {i + 2}"
                 for i in chunk.rows]
        preview = ", ".join(names[:3]) + (" ..." if len(names) > 3 else "")
        lines.append(f"chunk {chunk.index}: {chunk.size} simulation(s), "
                     f"~{chunk.estimated_minutes:.0f} min  [{preview}]")
    return "\n".join(lines)
