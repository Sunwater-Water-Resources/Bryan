"""Checks that run before anything is launched.

Every issue here is something that would otherwise surface as a KeyError or a
FileNotFoundError deep inside a run - after an hour, or overnight. Bryan does
catch them per row and carries on (Main.py:97-108), so a bad list does not stop
the batch; it just quietly produces nothing. ``RunLog.status_for`` exists
precisely because a missing input file is the usual culprit.

All of it is stat calls and column lookups, so it costs milliseconds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import completion as completion_module
from .columns import (EXCLUSION_KEYS, REPLICATE_FILE_ALIASES, REPLICATE_KEYS,
                      RESERVOIR_ROUTING, analyses_volumes, climate_requirement,
                      listed_keys, method_is_exact, normalise_method,
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
    issues.extend(_replication(frame, selected))
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
        for requirement in requirements_for(method, runs_models(row),
                                            analyses_volumes(row)):
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
        for requirement in requirements_for(method, runs_models(row),
                                            analyses_volumes(row)):
            name = requirement.present_in(columns)
            # Replicates, Replicate file, Exclusions and a reservoir routing
            # Config file are all legitimately blank - see BLANK_IS_NO_SETTING.
            if name is None or requirement.blank_ok:
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
        for requirement in requirements_for(method, runs_models(row),
                                            analyses_volumes(row)):
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


def _replication(frame, selected) -> list[Issue]:
    """The optional replication and exclusion settings, when they are used.

    Blank is the normal case and means "sample everything afresh", so nothing
    here fires on an empty cell. What does fire: naming a replicate without the
    file to read it from (``pd.read_csv(nan)``), and a key Bryan does not
    recognise - ``set_replicates``/``set_exclusions`` skip those in silence, so
    the run goes ahead sampling, or applying, exactly what it was told not to.
    """
    # None when neither spelling is there - _missing_columns has that already.
    replicate_column = next((candidate for candidate in REPLICATE_FILE_ALIASES
                             if candidate in frame.columns), None)

    no_file, unknown = [], {}
    for index in selected:
        row = frame.loc[index]
        method = normalise_method(row.get("Method"))
        if not method or method == RESERVOIR_ROUTING or not runs_models(row):
            continue
        for column, known in (("Replicates", REPLICATE_KEYS),
                              ("Exclusions", EXCLUSION_KEYS)):
            if column not in frame.columns:
                continue
            keys = listed_keys(row.get(column))
            for key in keys:
                if key not in known:
                    unknown.setdefault((column, key), []).append(index)
            if (column == "Replicates" and replicate_column is not None
                    and any(key in REPLICATE_KEYS for key in keys)
                    and is_blank(row.get(replicate_column))):
                no_file.append(index)

    issues = []
    if no_file:
        issues.append(Issue(
            BLOCK, "replicates-without-file",
            f"{len(no_file)} selected row(s) ask for replicated sampling but "
            f"leave {replicate_column!r} blank.",
            rows=tuple(no_file),
            fix_hint="Point it at the mcdf of the run being replicated, or "
                     "clear 'Replicates' to sample afresh. Simulator.__init__ "
                     "reads the file as soon as a replicate key is "
                     "recognised.",
        ))
    for (column, key), rows in sorted(unknown.items()):
        issues.append(Issue(
            WARN, "unknown-key",
            f"{len(rows)} selected row(s) list {key!r} under {column!r}, which "
            f"Bryan does not recognise. It is skipped without a message, so "
            f"the run looks normal and does not do it.",
            rows=tuple(rows),
            fix_hint=f"Keys for {column!r}: "
                     + ", ".join(REPLICATE_KEYS if column == "Replicates"
                                 else EXCLUSION_KEYS),
        ))
    return issues


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


# ---------------------------------------------------------------------------
# Running what can be run
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Triage:
    """The selection split into what can run and what cannot.

    A ninety-two row list with four bad input paths is the usual case, and
    hunting those four down by hand to deselect them is the whole friction this
    removes. It is not a safety feature: Bryan already catches a bad row and
    carries on (Main.py:97-108). It just means the run does not have to be
    started twice.
    """

    runnable: tuple = ()      # rows that pass every check, in sheet order
    skipped: dict = field(default_factory=dict)   # {row: the Issue that dropped it}
    remaining: tuple = ()     # blocking issues no row can be dropped to clear
    issues: tuple = ()        # every issue for `runnable`, warnings included

    @property
    def can_run(self) -> bool:
        return bool(self.runnable) and not self.remaining

    @property
    def is_whole_selection(self) -> bool:
        return not self.skipped

    def reasons(self) -> list:
        """[(row, code, message)] in row order - what to show the user."""
        return [(index, self.skipped[index].code, self.skipped[index].message)
                for index in sorted(self.skipped)]

    def codes(self) -> list:
        """The distinct reasons rows were dropped, most common first."""
        counts: dict = {}
        for issue in self.skipped.values():
            counts[issue.code] = counts.get(issue.code, 0) + 1
        return [code for code, _ in sorted(counts.items(),
                                           key=lambda item: (-item[1], item[0]))]


def triage(sims, config, selected_rows, *, completions=None, plan_for=None,
           max_passes=5) -> Triage:
    """Drop the rows that block, keep the rest.

    ``plan_for(rows) -> RunPlan`` folds the planner's own blocking hazards in,
    so a subset offered as runnable really is - a collision that only shows up
    once the rows are split into chunks would otherwise stop the run after the
    user had already been told it was fine. Left out, only the pre-flight
    checks are considered.

    Dropping is iterative because an issue can name several rows and clearing
    one can clear another: a collision goes away when its members go. It is
    also deliberately blunt - a collision drops *every* row that shares the
    output name, not the arbitrary all-but-one that would keep the most rows.
    Which of six identically-named rows the user wants is not the UI's to
    guess.
    """
    frame = sims.frame if hasattr(sims, "frame") else sims
    keep = [index for index in selected_rows if index in frame.index]
    skipped: dict = {}

    for _ in range(max_passes):
        issues = list(check(sims, config, keep, completions=completions))
        if keep and plan_for is not None:
            issues += list(plan_for(keep).hazards)

        blockers = [issue for issue in issues if issue.blocks]
        if not blockers:
            return Triage(runnable=tuple(keep), skipped=skipped,
                          issues=tuple(issues))

        drop: dict = {}
        for issue in blockers:
            for index in issue.rows:
                if index in keep:
                    drop.setdefault(index, issue)

        # A blocker naming no row of the selection - a missing column, the
        # 'Replicate file' spelling - cannot be skipped past by dropping rows.
        unclearable = tuple(issue for issue in blockers
                            if not any(index in drop for index in issue.rows))
        if unclearable or not drop:
            return Triage(skipped=skipped, remaining=unclearable or tuple(blockers),
                          issues=tuple(issues))

        skipped.update(drop)
        keep = [index for index in keep if index not in drop]
        if not keep:
            return Triage(skipped=skipped, issues=tuple(issues),
                          remaining=(Issue(BLOCK, "nothing-left",
                                           "Every selected row has a problem, "
                                           "so there is nothing to run."),))

    return Triage(skipped=skipped,
                  remaining=(Issue(BLOCK, "triage-unstable",
                                   f"Still blocked after dropping "
                                   f"{len(skipped)} row(s). Fix the problems "
                                   f"listed rather than skipping them."),))
