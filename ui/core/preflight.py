"""Checks that run before anything is launched.

Every issue here is something that would otherwise surface as a KeyError or a
FileNotFoundError deep inside a run - after an hour, or overnight. Bryan does
catch them per row and carries on (Main.py:97-108), so a bad list does not stop
the batch; it just quietly produces nothing. ``RunLog.status_for`` exists
precisely because a missing input file is the usual culprit.

All of it is stat calls and column lookups, so it costs milliseconds.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import completion as completion_module
from .columns import (REPLICATE_FILE_ALIASES, RESERVOIR_ROUTING,
                      climate_requirement, method_is_exact, normalise_method,
                      requirements_for, runs_models)
from .outputs import collision_key
from .paths import cell_text, is_blank, resolve_value

BLOCK = "block"
WARN = "warn"
INFO = "info"


@dataclass(frozen=True)
class Issue:
    severity: str
    code: str
    message: str
    rows: tuple = ()
    fix_hint: str = ""

    @property
    def blocks(self) -> bool:
        return self.severity == BLOCK


def check(sims, config, selected_rows, *, completions=None) -> list[Issue]:
    """Everything wrong with this selection. Empty means good to go."""
    frame = sims.frame
    selected = [index for index in selected_rows if index in frame.index]
    issues: list[Issue] = []

    if not selected:
        return [Issue(BLOCK, "empty-selection", "No simulations are selected.")]

    issues.extend(_uncached_formulas(sims, selected))
    issues.extend(_methods(frame, selected))
    issues.extend(_missing_columns(frame, selected))
    issues.extend(_blank_required(frame, selected))
    issues.extend(_missing_inputs(frame, config, selected))
    issues.extend(_collisions(frame, selected))
    issues.extend(_include_case(frame))
    issues.extend(_replicate_alias(frame))
    issues.extend(_needs_prior_results(frame, config, selected, completions))
    return issues


def _uncached_formulas(sims, selected) -> list[Issue]:
    audit = sims.audit
    if not audit.is_damaged:
        return []
    hit = tuple(index for index in selected if audit.affects(index))
    if not hit:
        return [Issue(
            WARN, "uncached-formulas-elsewhere",
            f"{sims.path.name}: {audit.describe()} None of the selected rows "
            f"are affected, but the workbook needs fixing.",
            fix_hint="Open in Excel, press Ctrl+Alt+F9, save.",
        )]
    return [Issue(
        BLOCK, "uncached-formulas",
        f"{sims.path.name}: {audit.describe()} "
        f"{len(hit)} selected row(s) would run with blank values.",
        rows=hit,
        fix_hint="Open the workbook in Excel, press Ctrl+Alt+F9 to "
                 "recalculate, save, then reload here.",
    )]


def _methods(frame, selected) -> list[Issue]:
    """Main.py:88-95 compares the Method cell exactly, so case matters."""
    bad = tuple(index for index in selected
                if not method_is_exact(frame.loc[index].get("Method")))
    if not bad:
        return []

    values = sorted({repr(frame.loc[index].get("Method")) for index in bad})
    recoverable = all(normalise_method(frame.loc[index].get("Method"))
                      for index in bad)
    detail = (
        "only the capitalisation is wrong, but Bryan does not lower the value "
        "before comparing it"
        if recoverable else "Bryan does not recognise this"
    )
    return [Issue(
        BLOCK, "bad-method",
        f"{len(bad)} row(s) have a Method Bryan would reject: "
        f"{', '.join(values)} - {detail}.",
        rows=bad,
        fix_hint="Use exactly 'monte carlo', 'ensemble' or "
                 "'reservoir routing', all lower case. Anything else raises "
                 "'Modelling method X not recognised' and loses the "
                 "simulation.",
    )]


def _missing_columns(frame, selected) -> list[Issue]:
    columns = set(frame.columns)
    missing: dict = {}
    for index in selected:
        row = frame.loc[index]
        method = normalise_method(row.get("Method"))
        if not method:
            continue
        for requirement in requirements_for(method, runs_models(row)):
            if requirement.present_in(columns) is None:
                missing.setdefault(requirement.name, []).append(index)

        kind, names = climate_requirement(row)
        if kind == "ssp" and method != RESERVOIR_ROUTING and runs_models(row):
            for name in names:
                if name not in columns:
                    missing.setdefault(name, []).append(index)

    return [
        Issue(BLOCK, "missing-column",
              f"The sims list has no {name!r} column, which Bryan reads for "
              f"{len(rows)} of the selected row(s).",
              rows=tuple(rows),
              fix_hint="A missing column is a KeyError inside the simulator, "
                       "which loses the whole simulation.")
        for name, rows in sorted(missing.items())
    ]


def _blank_required(frame, selected) -> list[Issue]:
    columns = set(frame.columns)
    blank: dict = {}
    for index in selected:
        row = frame.loc[index]
        method = normalise_method(row.get("Method"))
        if not method:
            continue
        for requirement in requirements_for(method, runs_models(row)):
            name = requirement.present_in(columns)
            # A blank Config file is legitimate for reservoir routing with
            # ensemble input - sim_list.md says to leave it blank.
            if name is None or (name == "Config file" and method == RESERVOIR_ROUTING):
                continue
            if is_blank(row.get(name)):
                blank.setdefault(name, []).append(index)

    return [
        Issue(BLOCK, "blank-required",
              f"{len(rows)} selected row(s) leave {name!r} blank, and Bryan "
              f"reads it for them.",
              rows=tuple(rows))
        for name, rows in sorted(blank.items())
    ]


def _missing_inputs(frame, config, selected) -> list[Issue]:
    issues = []
    for index in selected:
        row = frame.loc[index]
        method = normalise_method(row.get("Method"))
        if not method:
            continue
        missing = []
        for requirement in requirements_for(method, runs_models(row)):
            if not requirement.is_path:
                continue
            name = requirement.present_in(frame.columns)
            if name is None or is_blank(row.get(name)):
                continue
            resolved = resolve_value(config.project_folder, row.get(name))
            if resolved is not None and not resolved.exists():
                missing.append(f"{name}: {resolved}")
        if missing:
            label = cell_text(row.get("Output file")) or f"row {index + 2}"
            issues.append(Issue(
                BLOCK, "missing-input",
                f"{label}: {len(missing)} input file(s) not found.\n  "
                + "\n  ".join(missing),
                rows=(index,),
                fix_hint="Paths in a sims list resolve against the project "
                         "folder, which is where the batch file cd's to.",
            ))
    return issues


def _collisions(frame, selected) -> list[Issue]:
    """Rows that would write over each other - unsafe even run sequentially.

    In CLD_RFSL_mc_sims_01.xlsx twenty-four rows share one (Output file,
    Output suffix, Results folder), differing only by Duration. Only one row of
    that list is Include = yes, which is how it avoids the problem today.
    """
    buckets: dict = {}
    for index in selected:
        buckets.setdefault(collision_key(frame.loc[index]), []).append(index)

    issues = []
    for key, rows in sorted(buckets.items(), key=lambda item: -len(item[1])):
        if len(rows) < 2 or not key[0]:
            continue
        durations = ", ".join(
            cell_text(frame.loc[index].get("Duration")) or "?" for index in rows
        )
        issues.append(Issue(
            BLOCK, "output-collision",
            f"{len(rows)} selected rows share the output name {key[0]!r} "
            f"(durations: {durations}). They all write the same results files, "
            f"so only the last would survive - and if they run at the same "
            f"time they delete each other's working folders.",
            rows=tuple(rows),
            fix_hint="Give the rows distinct 'Output file' values - the name "
                     "needs a duration term - or run them one at a time.",
        ))
    return issues


def _include_case(frame) -> list[Issue]:
    """Main.py:64 compares Include without lowering, so 'Yes' silently skips."""
    if "Include" not in frame.columns:
        return []
    odd = frame["Include"].dropna().astype(str)
    wrong = sorted({value for value in odd
                    if value.strip().lower() == "yes" and value != "yes"})
    if not wrong:
        return []
    return [Issue(
        WARN, "include-case",
        f"The Include column contains {', '.join(repr(v) for v in wrong)}. "
        f"Bryan compares it exactly against 'yes' (Main.py:64), so those rows "
        f"are skipped when the list is run from a batch file.",
        fix_hint="The UI writes 'yes' into its own run copies, so this does "
                 "not affect runs launched from here.",
    )]


def _replicate_alias(frame) -> list[Issue]:
    correct, documented = REPLICATE_FILE_ALIASES
    if documented in frame.columns and correct not in frame.columns:
        return [Issue(
            BLOCK, "replicate-file-alias",
            f"The sims list has a {documented!r} column, which is what "
            f"Manual/SubDocs/sim_list.md documents - but Simulator.__init__ "
            f"reads {correct!r}. Every monte carlo row with 'Run models' set "
            f"would raise a KeyError.",
            fix_hint=f"Rename the column to {correct!r}. The manual is wrong.",
        )]
    return []


def _needs_prior_results(frame, config, selected, completions) -> list[Issue]:
    """The Run models = no inversion - see core/completion.py."""
    if completions is None:
        completions = completion_module.assess_frame(frame, config, rows=selected)
    rows = tuple(index for index in selected
                 if completions.get(index)
                 and completions[index].state == completion_module.NEEDS_PRIOR)
    if not rows:
        return []
    return [Issue(
        BLOCK, "needs-prior-results",
        f"{len(rows)} selected row(s) have 'Run models' = no, so they only "
        f"re-analyse existing results - but no results were found for them.",
        rows=rows,
        fix_hint="Run the row with 'Run models' = yes first, or point it at a "
                 "results folder that has them.",
    )]


def blocking(issues) -> list[Issue]:
    return [issue for issue in issues if issue.blocks]


def blocked_rows(issues) -> set:
    out: set = set()
    for issue in issues:
        if issue.blocks:
            out.update(issue.rows)
    return out
