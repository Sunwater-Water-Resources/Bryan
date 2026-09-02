"""Extract and plot the representative events chosen in the run launcher.

The launcher's Events page picks one realisation per design loading and saves
the list as ``<group>_representative_events.json``. This turns that list into
the things a study needs: the hydrographs of each chosen event, the rainfall
that produced it, and a plot of the three together.

    python util/RepresentativeEvents.py --config sims_config.json ^
        --selection sims_mc\\results\\GWL1p3_representative_events.json

Run with **Bryan's** interpreter: rebuilding the hyetograph drives
``lib/StormGenerator.py``, so scipy and matplotlib are needed. The selection
file itself is written by the launcher, which has neither - that split is the
point of ``lib/RepresentativeEvents.py``, which both sides share so the event
that gets plotted is the event that was chosen.

**The hyetograph is rebuilt, not read.** An mcdf records what was sampled, not
the rainfall series, so ``lib/EventStorm.py`` replays the storm generation for
that one realisation and checks the result against the depths the run itself
recorded. A rebuild that does not match is reported on the plot and in the
workbook rather than passed off as the storm that was modelled.

Time is measured **from the start of the main burst**, so the pre-burst runs at
negative times and events with different pre-burst durations line up with each
other. The stored hydrographs start at the beginning of the storm file, which
is the start of the *pre-burst*, so they are shifted by the rebuilt pre-burst
duration. Where there is no rebuild to take it from, the shift falls back to
the simulation period in the model config, which Bryan lengthens by exactly
that amount (``Simulator.run_models``: ``simulation_period += preburst_duration``).

``util/GetRepresentativeEvents.py`` is the older script that both chooses and
extracts, driven by a workbook of hard-coded paths. It still works; this is the
half of it that the launcher hands off to.
"""

import argparse
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib import EventStorm                                       # noqa: E402
from lib import RepresentativeEvents as events                   # noqa: E402

SERIES_LABELS = {'inflows': 'Inflow', 'outflows': 'Outflow', 'levels': 'Lake level'}

FLOW_COLOURS = {'inflows': '#4E79A7', 'outflows': '#E15759'}
RAIN_COLOUR = '#76B7B2'
LEVEL_COLOUR = '#59A14F'


# -- the project -------------------------------------------------------------

def load_config(path):
    """A sims_config.json, resolved exactly as Main.py resolves it."""
    with open(path) as handle:
        config = json.load(handle)

    folder = os.path.dirname(os.path.abspath(path))
    project_folder = os.path.dirname(path)
    if str(config.get('project_folder', 'default')).lower() != 'default':
        project_folder = config['project_folder']
    project_folder = os.path.abspath(project_folder)

    filepaths = {key: os.path.normpath(os.path.join(folder, value))
                 for key, value in config.get('filepaths', {}).items()}

    sims_list = os.path.join(project_folder, config['simulation_list'])
    return {'project_folder': project_folder, 'filepaths': filepaths,
            'sims_list': sims_list, 'raw': config}


def read_sims_list(path):
    frame = pd.read_excel(path, sheet_name=0)
    frame.columns = [str(name).strip() for name in frame.columns]
    return frame


def resolve(project_folder, value):
    """A sims-list path against the project folder, on either separator."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip()
    if not text:
        return None
    text = text.replace('\\', os.sep).replace('/', os.sep)
    return text if os.path.isabs(text) else os.path.join(project_folder, text)


def find_row(frame, output_file):
    """The sims-list row that produced an event, by its Output file."""
    if 'Output file' not in frame.columns or not output_file:
        return None
    wanted = os.path.basename(str(output_file).replace('\\', '/').rstrip('/'))
    for index in frame.index:
        value = frame.loc[index].get('Output file')
        if pd.isna(value):
            continue
        if os.path.basename(str(value).replace('\\', '/')) == wanted:
            return frame.loc[index]
    return None


def cell(row, name, default=None):
    if row is None or name not in row.index:
        return default
    value = row[name]
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    return value


# -- the stored series -------------------------------------------------------

def read_series(path, column):
    """One simulation's column out of a stored hydrograph file.

    Ragged files are ordinary here - an ensemble file pads every duration out
    to the longest - so the trailing NaN are dropped rather than filled.
    """
    if not path or not os.path.isfile(path):
        return None
    frame = pd.read_csv(path, index_col=0)
    if column not in frame.columns:
        return None
    series = pd.to_numeric(frame[column], errors='coerce').dropna()
    series.index = pd.to_numeric(series.index, errors='coerce')
    return series[series.index.notna()]


def simulation_period(model_config, duration):
    """UrbsModel.get_simulation_period, without instantiating the model.

    Deliberately not through ``UrbsModel``: its constructor rmtree's the run's
    working folder, which is not a thing a post-processing script may do.
    """
    if not model_config or not os.path.isfile(model_config):
        return None
    try:
        with open(model_config) as handle:
            data = json.load(handle)
        periods = data.get('simulation_periods', {})
    except (OSError, ValueError):
        return None
    key = str(duration)
    if key in periods:
        return float(periods[key])
    for candidate in (f'{float(duration):g}', str(int(float(duration)))):
        if candidate in periods:
            return float(periods[candidate])
    return float(duration) * 2


def burst_shift(series, hyeto, period):
    """How far to slide a stored series so t = 0 is the start of the main burst.

    Two independent ways of getting the same number, and they are reported
    against each other: the rebuilt pre-burst duration, and how much longer the
    run was than its nominal simulation period.
    """
    from_storm = hyeto.preburst_hours if hyeto is not None else None
    from_period = None
    if period is not None and series is not None and len(series):
        from_period = float(series.index.max()) - float(period)
        if from_period < 0:
            from_period = None

    if from_storm is None:
        return (from_period or 0.0), ('simulation period' if from_period else 'none')
    if from_period is not None and abs(from_period - from_storm) > max(0.5, 0.02 * from_storm):
        print(f'  NOTE: the rebuilt pre-burst is {from_storm:.2f} h but the run is '
              f'{from_period:.2f} h longer than its simulation period. Using the rebuild.')
    return from_storm, 'rebuilt pre-burst'


# -- one event ---------------------------------------------------------------

def collect(target, frame, config, storms, rebuild=True):
    """Everything about one chosen event: its series, its rainfall, its notes."""
    row = find_row(frame, target.output_file)
    project = config['project_folder']
    sim_id = int(target.picked)
    column = events.sim_label(sim_id)
    notes = []

    database = resolve(project, target.database)
    if database is None and target.output_file:
        database = resolve(project, f'{target.output_file}__mcdf.csv')
    sim = None
    if database and os.path.isfile(database):
        mcdf = events.load_mcdf(database)
        if sim_id in mcdf.index:
            sim = mcdf.loc[sim_id]
        else:
            notes.append(f'simulation {sim_id} is not in {os.path.basename(database)}')
    else:
        notes.append('the results database was not found, so the event has no metrics')

    method = str(cell(row, 'Method', 'monte carlo')).strip().lower()
    hydrographs = events.hydrograph_paths(
        cell(row, 'Output file', target.output_file),
        model=str(cell(row, 'Model type', 'urbs')),
        hydrographs_folder=(resolve(project, cell(row, 'Hydrographs folder'))
                            if method == 'reservoir routing' else None),
        suffix=str(cell(row, 'Output suffix', '') or ''),
        inflow=resolve(project, cell(row, 'Inflow')),
    )

    series = {}
    for kind, path in hydrographs.items():
        found = read_series(resolve(project, path), column)
        if found is not None:
            series[kind] = found
        else:
            notes.append(f'no {kind} for {column} in {os.path.basename(str(path))}')

    hyeto = None
    if rebuild and sim is not None:
        try:
            hyeto = rebuild_hyetograph(row, sim, config, storms)
        except Exception as error:                    # noqa: BLE001 - one event, not the run
            notes.append(f'the hyetograph could not be rebuilt ({error})')
    if hyeto is not None and not hyeto.trustworthy:
        notes.extend(f'hyetograph does not match the run - {problem}'
                     for problem in hyeto.problems)

    period = simulation_period(config['filepaths'].get('model_config'),
                              cell(row, 'Duration', 0))
    # Any one of them dates the run - but `a or b` on a Series raises, so pick
    # the first that is actually there rather than leaning on truthiness.
    reference = next((series[kind] for kind in ('inflows', 'levels', 'outflows')
                      if kind in series), None)
    shift, source = burst_shift(reference, hyeto, period)
    series = {kind: values.set_axis(values.index - shift)
              for kind, values in series.items()}

    return {'target': target, 'row': row, 'sim': sim, 'sim_id': sim_id,
            'column': column, 'series': series, 'hyetograph': hyeto,
            'shift': shift, 'shift_source': source, 'notes': notes}


def rebuild_hyetograph(row, sim, config, storms):
    """The storm behind one realisation, reusing the run's rainfall data."""
    duration = float(cell(row, 'Duration', 0) or 0)
    if not duration:
        raise ValueError('the row has no Duration')

    key = (str(cell(row, 'Output file', '')), duration)
    if key not in storms:
        context = EventStorm.StormContext(
            storm_config=config['filepaths'].get('storm_config'),
            focal_subcatchments=resolve(config['project_folder'],
                                        cell(row, 'Focal subcatchments')),
            duration=duration,
            climate_config=config['filepaths'].get('climate_config'),
            gwl=cell(row, 'GWL'), year=cell(row, 'Year'), ssp=cell(row, 'SSP'),
            initial_loss=cell(row, 'IL'), continuing_loss=cell(row, 'CL'),
            exclusions=cell(row, 'Exclusions', ''),
            preburst_method=cell(row, 'Pre-burst method'),
            **_aep_range(resolve(config['project_folder'], cell(row, 'Config file'))),
        )
        storms[key] = (context,) + EventStorm.build_storm(context)
    context, storm, climate = storms[key]
    return EventStorm.hyetograph(storm, climate, context, sim)


def _aep_range(method_config):
    """The scheme's AEP bounds, which decide whether extreme rainfall is set up."""
    try:
        with open(method_config) as handle:
            data = json.load(handle)
    except (OSError, TypeError, ValueError):
        return {}
    scheme = data.get('scheme_config', data)
    bounds = {}
    for key in ('lower_aep', 'upper_aep'):
        if key in scheme:
            bounds[key] = float(scheme[key])
    return bounds


# -- output ------------------------------------------------------------------

def safe_name(text):
    return ''.join(c if c.isalnum() or c in '-_.' else '_' for c in str(text)).strip('_')


def plot_event(event, path, dpi=150):
    """Draw one event and write it out."""
    figure = figure_for(event)
    figure.savefig(path, dpi=dpi)
    plt.close(figure)


def figure_for(event):
    """Hyetograph, flows and lake level on one time axis.

    Separate from writing it so the figure can be inspected - the rainfall
    axis has to come out *reversed*, and getting that wrong looks like a
    perfectly ordinary plot in every test that only checks a file was written.
    """
    target = event['target']
    hyeto = event['hyetograph']
    series = event['series']

    figure, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True,
                                gridspec_kw={'height_ratios': [1.0, 1.6, 1.4]})

    # Rainfall, drawn downwards from the top - the convention, and it keeps the
    # storm out of the way of the hydrograph rising underneath it.
    rain = axes[0]
    if hyeto is not None and len(hyeto.depths):
        rain.bar(hyeto.depths.index, hyeto.depths.to_numpy(),
                 width=hyeto.timestep, align='edge', color=RAIN_COLOUR,
                 edgecolor='none')
        rain.set_ylabel(f'Rainfall\n(mm per {hyeto.timestep:g} h)')
        # Reversed: zero at the top, the storm hanging down from it. Set as a
        # pair with the larger bound first - `set_ylim(top=0)` on its own
        # collapses the range to nothing, since the other bound is still zero.
        deepest = float(hyeto.depths.max())
        rain.set_ylim(deepest * 1.15 if deepest > 0 else 1.0, 0.0)
        if hyeto.preburst_hours:
            rain.axvline(0.0, color='k', alpha=0.4, lw=0.8)
            rain.annotate('main burst', xy=(0, 0), xytext=(4, 4),
                          textcoords='offset points', fontsize=8, alpha=0.7)
    else:
        rain.set_ylabel('Rainfall')
        rain.text(0.5, 0.5, 'no hyetograph', transform=rain.transAxes,
                  ha='center', va='center', fontsize=9, alpha=0.5)
    rain.grid(alpha=0.25)

    flows = axes[1]
    for kind in ('inflows', 'outflows'):
        if kind in series:
            flows.plot(series[kind].index, series[kind].to_numpy(),
                       color=FLOW_COLOURS[kind], label=SERIES_LABELS[kind], lw=1.6)
    flows.set_ylabel('Flow (m³/s)')
    flows.grid(alpha=0.25)
    if flows.get_legend_handles_labels()[0]:
        flows.legend(loc='upper right', frameon=False)

    level = axes[2]
    if 'levels' in series:
        level.plot(series['levels'].index, series['levels'].to_numpy(),
                   color=LEVEL_COLOUR, lw=1.6)
    full_supply = cell(event['row'], 'FSL')
    if full_supply is not None:
        try:
            level.axhline(float(full_supply), color='k', ls='--', lw=0.8, alpha=0.5)
            level.annotate('FSL', xy=(0.01, float(full_supply)),
                           xycoords=('axes fraction', 'data'), fontsize=8,
                           va='bottom', alpha=0.7)
        except (TypeError, ValueError):
            pass
    level.set_ylabel('Lake level (m AHD)')
    level.set_xlabel('Time from the start of the main burst (hours)'
                     if event['shift'] else 'Time (hours)')
    level.grid(alpha=0.25)

    figure.suptitle(_title(event), fontsize=11)
    if event['notes']:
        figure.text(0.01, 0.005, '\n'.join(event['notes'][:3]), fontsize=7,
                    color='#B03A2E', va='bottom')
    figure.tight_layout(rect=(0, 0.03 if event['notes'] else 0, 1, 0.97))
    return figure


def _title(event):
    target = event['target']
    sim = event['sim']
    bits = [f'{target.label} {target.result_type}', f'sim {event["sim_id"]}']
    if sim is not None:
        rain = sim.get('rain_aep')
        if rain is not None and not pd.isna(rain):
            bits.append(f'rain 1 in {float(rain):,.0f}')
        column = f'{target.result_type}_aep'
        if column in sim.index and not pd.isna(sim[column]) and float(sim[column]) > 0:
            bits.append(f'{target.result_type} 1 in {1 / float(sim[column]):,.0f}')
    if event['target'].source:
        bits.append(event['target'].source)
    return '  |  '.join(bits)


def summary_frame(collected):
    rows = []
    for event in collected:
        target = event['target']
        sim = event['sim']
        hyeto = event['hyetograph']
        row = {'loading': target.label, 'result': target.result_type,
               'source': target.source, 'output file': target.output_file,
               'sim': event['sim_id'], 'hydrograph': event['column']}
        for name in ('rain_aep', 'mean_rain_mm', 'preburst_mm', 'ADV',
                     'inflow', 'level', 'outflow', 'embedded_bursts'):
            row[name] = None if sim is None else sim.get(name)
        if hyeto is not None:
            row['rebuilt burst mm'] = round(hyeto.burst_mm, 2)
            row['rebuilt preburst mm'] = round(hyeto.preburst_mm, 2)
            row['preburst hours'] = hyeto.preburst_hours
            row['hyetograph checks'] = ('matches the run' if hyeto.trustworthy
                                        else '; '.join(hyeto.problems))
        row['time shift (h)'] = round(event['shift'], 3)
        row['notes'] = '; '.join(event['notes'])
        rows.append(row)
    return pd.DataFrame(rows)


def write_workbook(collected, path):
    """One sheet per series, one column per event, plus what was done."""
    with pd.ExcelWriter(path) as writer:
        summary_frame(collected).to_excel(writer, sheet_name='events', index=False)
        for kind in events.HYDROGRAPH_KINDS:
            columns = {}
            for event in collected:
                if kind in event['series']:
                    columns[_column_name(event)] = event['series'][kind]
            if columns:
                pd.concat(columns, axis=1).sort_index().to_excel(writer, sheet_name=kind)
        rain = {_column_name(event): event['hyetograph'].depths
                for event in collected if event['hyetograph'] is not None}
        if rain:
            pd.concat(rain, axis=1).sort_index().to_excel(writer, sheet_name='hyetographs')


def _column_name(event):
    return f'{safe_name(event["target"].label)}_sim{event["sim_id"]:05d}'


# -- entry point -------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Extract and plot the representative events chosen in the '
                    'run launcher.')
    parser.add_argument('--config', required=True,
                        help='the sims_config.json of the project the events came from')
    parser.add_argument('--selection', required=True,
                        help='the <group>_representative_events.json the launcher saved')
    parser.add_argument('--out', default=None,
                        help='where to write (default: beside the selection file)')
    parser.add_argument('--name', default=None,
                        help='basename for the outputs (default: the selection file stem)')
    parser.add_argument('--no-hyetograph', action='store_true',
                        help='skip rebuilding the rainfall - much faster, and needed '
                             'where the storm data is not to hand')
    parser.add_argument('--no-plots', action='store_true', help='write the workbook only')
    parser.add_argument('--dpi', type=int, default=150)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    config = load_config(args.config)
    frame = read_sims_list(config['sims_list'])
    targets, settings = events.read_selection(args.selection)

    chosen = [target for target in targets if target.picked is not None]
    if not chosen:
        print('No events have been chosen in', args.selection)
        print('Open the Events page in the launcher, pick an event for each '
              'loading, and press Save.')
        return 1

    folder = args.out or os.path.dirname(os.path.abspath(args.selection))
    os.makedirs(folder, exist_ok=True)
    stem = args.name or os.path.splitext(os.path.basename(args.selection))[0]

    print(f'{len(chosen)} chosen event(s) from {os.path.basename(args.selection)}')
    if settings:
        print('  settings:', ', '.join(f'{k}={v}' for k, v in settings.items()))

    storms = {}
    collected = []
    for target in chosen:
        print(f'\n{target.label} ({target.result_type}) - simulation {target.picked}')
        event = collect(target, frame, config, storms,
                        rebuild=not args.no_hyetograph)
        for note in event['notes']:
            print('  NOTE:', note)
        found = ', '.join(sorted(event['series'])) or 'nothing'
        print(f'  series: {found}; time shifted by {event["shift"]:.2f} h '
              f'({event["shift_source"]})')
        collected.append(event)

        if not args.no_plots:
            path = os.path.join(folder,
                                f'{stem}_{safe_name(target.label)}_sim{event["sim_id"]:05d}.png')
            plot_event(event, path, dpi=args.dpi)
            print('  wrote', path)

    workbook = os.path.join(folder, f'{stem}_events.xlsx')
    write_workbook(collected, workbook)
    print('\nwrote', workbook)
    return 0


if __name__ == '__main__':
    sys.exit(main())
