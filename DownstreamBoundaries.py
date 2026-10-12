"""Put a dam failure model's outflow into the regional model's inflow files.

    python DownstreamBoundaries.py --selection sims_mc\\results\\GWL1p3_representative_events.json ^
        --model runs\\Regional_E001\\urbs_config_rev2.json ^
        --dam-outflow ..\\dam_failure\\Simulation_List.xlsx --out runs\\Regional_E001\\boundaries

``DownstreamStorms.py --run-urbs`` runs the regional model for each chosen event, and URBS
writes ``<run>_i.csv`` - the flow at every print location marked for it - as the upstream
boundaries of a hydraulic model. One of those locations is the dam, and its column is the
regional model's own routing of the dam's outflow.

Where the dam is to fail, or simply where the study routes the dam in a dam failure model
(Judyth) so that the no-failure case and the failures come from one model, that column has
to come from there instead. This writes a copy of each event's inflow file with the dam's
column replaced, and leaves URBS's own files alone.

**The dam failure model's simulation list says which run belongs to which event.** It is
read for the rows whose ``Inflow_id`` ends in the event's realisation (``sim08029``), whose
``Bryan_group`` - where the list carries that column - is the group the selection was saved
for, and whose ``Gate_scenario`` is the one asked for (blank, by default: the gates
working). Exactly one row has to be left. ``--name-contains`` narrows it where a list holds
the same event more than once.

**The dam's run has to cover the regional model's.** The regional model runs on after the
upstream one stops, for the flood to pass down the valley, and the dam is still releasing
over that time. A dam run that stops early is refused rather than padded: extend it.

Pandas only - no scipy - but run with Bryan's interpreter like the rest.
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import RepresentativeEvents as events                    # noqa: E402
from lib.ConfigPaths import resolve                               # noqa: E402

DURATION = re.compile(r'_(\d+(?:p\d+)?)h_')
GWL = re.compile(r'GWL(\d+(?:p\d+)?)')
SELECTION_SUFFIX = '_' + events.SELECTION_FILE


def from_name(pattern, name, override, what):
    if override is not None:
        return override
    found = pattern.search(os.path.basename(str(name)))
    if not found:
        raise SystemExit(f'cannot read the {what} from "{name}" - pass it explicitly')
    return float(found.group(1).replace('p', '.'))


def run_name(picked, duration, gwl, suffix=''):
    """The regional run for an event - the storm file's name without its extension."""
    return f'{int(picked):06d}_{duration:g}h_GWL{gwl:g}'.replace('.', 'p') + suffix


def group_of(selection):
    """The group a selection was saved for, from its file name; '' where it carries none."""
    name = os.path.basename(str(selection))
    return name[:-len(SELECTION_SUFFIX)] if name.endswith(SELECTION_SUFFIX) else ''


def _blank(value):
    return value is None or (isinstance(value, float) and np.isnan(value)) or str(value).strip() == ''


def dam_run(sims, picked, group='', gate_scenario='', name_contains=''):
    """The simulation list row holding the dam's run for one event."""
    tail = 'sim{:05d}'.format(int(picked))
    rows = sims[sims['Inflow_id'].astype(str).str.endswith(tail)]
    if group and 'Bryan_group' in rows.columns:
        rows = rows[rows['Bryan_group'].astype(str) == group]
    if 'Gate_scenario' in rows.columns:
        wanted = str(gate_scenario or '').strip()
        rows = rows[[(_blank(v) and not wanted) or str(v).strip() == wanted
                     for v in rows['Gate_scenario']]]
    if name_contains:
        rows = rows[rows['Simulation_name'].astype(str).str.contains(name_contains, regex=False)]
    if len(rows) != 1:
        found = ', '.join(rows['Simulation_name'].astype(str)) or 'none'
        raise SystemExit(
            f'{len(rows)} simulations in the dam failure list match realisation {picked}'
            + (f' in group {group}' if group else '')
            + (f' with gate scenario "{gate_scenario}"' if gate_scenario else ' with the gates working')
            + f': {found}. Exactly one is needed - use --name-contains, or --gate-scenario.')
    return rows.iloc[0]


def replace_dam_column(inflow, column, dam_times, dam_flows, what=''):
    """A regional inflow table with one column replaced by the dam model's outflow."""
    if column not in inflow.columns:
        raise SystemExit(
            f'{what} has no "{column}" column. URBS writes a print location to the inflow file '
            'only where its PRINT statement in the vec carries an asterisk after the name.')
    times = inflow.iloc[:, 0].to_numpy(dtype=float)
    if times[-1] > float(dam_times[-1]) + 1e-6:
        raise SystemExit(
            f'the dam run for {what} stops at {float(dam_times[-1]):g} h and the regional run '
            f'goes on to {times[-1]:g} h. Extend the dam run to the regional period.')
    out = inflow.copy()
    out[column] = np.round(np.interp(times, dam_times, dam_flows), 2)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selection', required=True,
                    help='the <group>_representative_events.json the launcher saved')
    ap.add_argument('--model', required=True, help='the regional model URBS config')
    ap.add_argument('--dam-outflow', required=True,
                    help="the dam failure model's simulation list (.xlsx)")
    ap.add_argument('--out', required=True, help='the folder the inflow files are written to')
    ap.add_argument('--gate-scenario', default='',
                    help='the Gate_scenario of the dam runs to use; blank is the gates working')
    ap.add_argument('--name-contains', default='',
                    help='only dam simulations whose name contains this')
    ap.add_argument('--suffix', default='', help='the suffix the storm files were written with')
    ap.add_argument('--duration', type=float, help='override the duration read from the name')
    ap.add_argument('--gwl', type=float, help='override the warming level read from the name')
    ap.add_argument('--group', help='the group, where the selection file name does not carry it')
    ap.add_argument('--outflow-column', default='Total Outflow',
                    help="the outflow column of the dam model's result files")
    ap.add_argument('--time-column', default='Time')
    args = ap.parse_args(argv)

    with open(args.model) as handle:
        model = json.load(handle)
    model_folder = resolve(os.path.dirname(os.path.abspath(args.model)), model['model_folder'])
    column = model.get('max_keys', {}).get('outflow')
    if not column:
        raise SystemExit('the model config names no max_keys "outflow", so there is no dam '
                         'outflow location to replace')

    targets, _ = events.read_selection(args.selection)
    chosen = [t for t in targets if t.picked is not None]
    if not chosen:
        raise SystemExit('nothing to do - no loading has an event picked')
    group = args.group if args.group is not None else group_of(args.selection)

    sims = pd.read_excel(args.dam_outflow)
    missing = [c for c in ('Simulation_name', 'Inflow_id', 'Output_name') if c not in sims.columns]
    if missing:
        raise SystemExit(f'{args.dam_outflow} is not a dam failure simulation list: it has no '
                         + ', '.join(missing))
    dam_folder = os.path.dirname(os.path.abspath(args.dam_outflow))

    os.makedirs(args.out, exist_ok=True)
    tag = f'_{args.gate_scenario}' if args.gate_scenario else ''
    records = []
    print(f'{len(chosen)} events from {os.path.basename(args.selection)}; '
          f'the dam outflow location is {column}')
    for target in chosen:
        duration = from_name(DURATION, target.database, args.duration, 'duration')
        gwl = from_name(GWL, target.database, args.gwl, 'warming level')
        run = run_name(target.picked, duration, gwl, args.suffix)
        source = os.path.join(model_folder, f'{run}_i.csv')
        if not os.path.isfile(source):
            raise SystemExit(f'no regional inflow file for realisation {target.picked}: {source}. '
                             'Run the regional model first (DownstreamStorms.py --run-urbs), '
                             'with the same suffix.')
        row = dam_run(sims, target.picked, group, args.gate_scenario, args.name_contains)
        result = resolve(dam_folder, row['Output_name'])
        if not os.path.isfile(result):
            raise SystemExit(f'{row["Simulation_name"]} has not been run: no {result}')
        dam = pd.read_csv(result)
        for needed in (args.time_column, args.outflow_column):
            if needed not in dam.columns:
                raise SystemExit(f'{result} has no "{needed}" column')
        inflow = pd.read_csv(source)
        urbs_peak = float(inflow[column].max()) if column in inflow.columns else float('nan')
        replaced = replace_dam_column(inflow, column, dam[args.time_column].to_numpy(dtype=float),
                                      dam[args.outflow_column].to_numpy(dtype=float), what=run)
        name = f'{run}{tag}_i.csv'
        replaced.to_csv(os.path.join(args.out, name), index=False, float_format='%.3f')
        records.append(dict(
            regional_run=run, loading=f'{target.kind} {target.value:g}', realisation=int(target.picked),
            dam_simulation=str(row['Simulation_name']), gate_scenario=args.gate_scenario,
            inflow_file=name, subarea_runoff_file=os.path.join(model_folder, f'{run}_r.csv'),
            dam_outflow_column=column, run_ends_h=float(inflow.iloc[-1, 0]),
            peak_dam_outflow=float(dam[args.outflow_column].max()),
            peak_dam_outflow_urbs=urbs_peak,
            time_of_peak_dam_outflow_h=float(dam.loc[dam[args.outflow_column].idxmax(), args.time_column])))
        print(f'  wrote {name}  <- {row["Simulation_name"]}')

    stem = group or os.path.splitext(os.path.basename(args.selection))[0]
    index = os.path.join(args.out, f'{stem}{tag}_boundaries.csv')
    pd.DataFrame(records).to_csv(index, index=False)
    print(f'\n{len(records)} inflow files in {args.out}; what went into them: {index}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
