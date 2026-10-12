"""Write regional storm files for the representative events chosen in the launcher.

    python DownstreamStorms.py --config downstream_storm_config_CLD_01.json ^
        --selection sims_mc\\results\\GWL1p3_representative_events.json ^
        --model runs\\Regional_E001\\urbs_config_rev2.json

The launcher's Events page picks one realisation per design loading and saves
the list as ``<group>_representative_events.json``. This turns that list into a
storm file per event for the downstream regional model - the same handover
``util/RepresentativeEvents.py`` makes for plotting, pointed at a different
consumer.

**The event list is the input, not a workbook of its own.** That is the whole
difference from ``DownstreamStormGenerator.py``, which this replaces: re-picking
events on the Events page and re-running this is the entire repeat loop, and
nothing else has to be edited for a new set.

Run with **Bryan's** interpreter: assembling the depths drives
``lib/StormGenerator.py``, so scipy is needed. The launcher has neither scipy nor
matplotlib and so runs this as a subprocess, the same way it runs ``Main.py``.
"""
import argparse
import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import RepresentativeEvents as events                    # noqa: E402
from lib.ConfigPaths import resolve                               # noqa: E402
from lib.DownstreamStorms import DownstreamRainfall, DownstreamStormWriter   # noqa: E402

#: Durations and warming levels are not on a saved target, so they are read out
#: of the database's own name - the naming the runs already use, e.g.
#: ``CLD_mc_24h_E009_..._GWL1p3_RGN__mcdf.parquet``. Both can be overridden.
DURATION = re.compile(r'_(\d+(?:p\d+)?)h_')
GWL = re.compile(r'GWL(\d+(?:p\d+)?)')


def number(text):
    return float(text.replace('p', '.'))


def from_name(pattern, name, override, what):
    if override is not None:
        return override
    found = pattern.search(os.path.basename(str(name)))
    if not found:
        raise SystemExit(f'cannot read the {what} from "{name}" - pass it explicitly')
    return number(found.group(1))


def storm_name(target, duration, gwl, suffix):
    # 'p' for the decimal point in the stem, the way the runs already name
    # durations and warming levels; the extension is URBS's, and is the
    # duration as written.
    stem = f'{int(target.picked):06d}_{duration:g}h_GWL{gwl:g}'.replace('.', 'p')
    return f'{stem}{suffix}.{duration:g}'


#: What the upstream run records against a realisation whose pattern it changed.
FILTERED = 'Embedded burst filtered'


class PatternRebuild:
    """Each event's temporal pattern as the upstream run applied it.

    The database records what was sampled - which pattern, at what rainfall - and not the
    series that came of it. Where the run filtered an embedded burst the two differ: the
    same burst depth, spread to keep every sub-duration inside its own IFD depth. The
    regional storm has to carry the pattern the event had, or the dam catchment nested in
    the regional model gets a sharper storm than the one that was chosen, and a larger
    flood. ``util/RepresentativeEvents.py`` already replays the storm generation for one
    realisation and checks the result against the depths the run recorded; this uses the
    same code, so the pattern written here is the one it plots.
    """

    def __init__(self, sims_config):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'util'))
        import RepresentativeEvents as extract
        self.extract = extract
        self.config = extract.load_config(sims_config)
        self.frame = extract.read_sims_list(self.config['sims_list'])
        self.storms = {}

    def pattern(self, target, event):
        """Percent of the burst per timestep, or None where it cannot be rebuilt."""
        row = self.extract.find_row(self.frame, target.output_file)
        source, note = self.extract.storm_row(self.frame, row)
        if source is None:
            print(f'  NOTE: the pattern was not rebuilt - {note}')
            return None
        hyeto = self.extract.rebuild_hyetograph(source, event, self.config, self.storms)
        if not hyeto.trustworthy:
            raise SystemExit(f'the storm rebuilt for realisation {target.picked} does not match '
                             f'what the run recorded: ' + '; '.join(hyeto.problems))
        burst = hyeto.depths[hyeto.depths.index >= 0]
        return 100.0 * burst / burst.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, help='the downstream storm config')
    ap.add_argument('--selection', required=True,
                    help='the <group>_representative_events.json the launcher saved')
    ap.add_argument('--model', required=True, help='the regional model URBS config')
    ap.add_argument('--storm-config', help='the upstream storm config (temporal patterns)')
    ap.add_argument('--climate-config')
    ap.add_argument('--duration', type=float, help='override the duration read from the name')
    ap.add_argument('--gwl', type=float, help='override the warming level read from the name')
    ap.add_argument('--rain-lag', type=float, default=0.0, help='hours to shift the rain')
    ap.add_argument('--suffix', default='')
    ap.add_argument('--dry-run', action='store_true',
                    help='report what would be written, and write nothing')
    ap.add_argument('--run-urbs', action='store_true',
                    help='run URBS on each storm as it is written. Without this the batch '
                         'file is still written and can be run by hand or in a batch later')
    ap.add_argument('--sims-config',
                    help="the upstream run's sims config - the one the Events page had open. "
                         "With it each event's temporal pattern is rebuilt as the run applied "
                         "it, embedded burst filtering included. Needed for any event whose "
                         "embedded bursts were filtered")
    ap.add_argument('--project',
                    help='the project folder the selection was saved in. Only needed when '
                         'its databases are not found from the selection file itself')
    ap.add_argument('--lake-config',
                    help="the main dam's lake config, used only where a realisation has no "
                         "ADV column of its own")
    args = ap.parse_args()
    # URBS is run from the model folder, by a batch file that changes into it first, so a
    # model config given relative to where this was typed points nowhere once it has.
    args.model = os.path.abspath(args.model)

    with open(args.config) as handle:
        cfg = json.load(handle)
    # Through resolve, not os.path.join: these configs are written on Windows
    # and record their paths with backslashes.
    folder = os.path.dirname(os.path.abspath(args.config))
    storm_config = args.storm_config or resolve(folder, cfg['file_paths']['storm_config'])
    climate = args.climate_config or cfg['file_paths'].get('climate_config')
    if climate:
        climate = climate if os.path.isabs(climate) else resolve(folder, climate)

    targets, settings = events.read_selection(args.selection)
    chosen = [t for t in targets if t.picked is not None]
    print(f'{len(targets)} loadings in {os.path.basename(args.selection)}, '
          f'{len(chosen)} with an event picked')
    if not chosen:
        raise SystemExit('nothing to do - no loading has an event picked')

    if args.dry_run and args.run_urbs:
        raise SystemExit('--dry-run writes nothing, so there is no storm for --run-urbs '
                         'to run. Choose one.')

    rainfall = DownstreamRainfall(args.config)
    lake_config = args.lake_config or cfg['file_paths'].get('lake_config')
    if lake_config and not os.path.isabs(lake_config):
        lake_config = resolve(folder, lake_config)
    writer = None if args.dry_run else DownstreamStormWriter(
        rainfall, storm_config, args.model, climate, main_lake_config=lake_config)

    rebuild = PatternRebuild(args.sims_config) if args.sims_config else None

    records = []
    for target in chosen:
        duration = from_name(DURATION, target.database, args.duration, 'duration')
        gwl = from_name(GWL, target.database, args.gwl, 'warming level')
        # Saved relative to the project folder, which is not always where this is run.
        database = events.locate_database(target.database, args.selection, args.project)
        if not os.path.isfile(database):
            raise SystemExit(f'the database for realisation {target.picked} is not on disk: '
                             f'"{target.database}". It is saved relative to the project '
                             f'folder - pass that as --project')
        event = DownstreamStormWriter.read_event(database, int(target.picked))
        name = storm_name(target, duration, gwl, args.suffix)
        print(f'\n{target.kind} {target.value:g} ({target.result_type}) -> realisation '
              f'{target.picked}, {duration:g} h, GWL {gwl:g}')
        if args.dry_run:
            records.append(dict(filename=name, duration=duration, gwl=gwl,
                                driver_aep=float(event['rain_aep']),
                                storm_method=str(event['storm_method'])))
            print(f'  would write {name}')
            continue
        pattern = rebuild.pattern(target, event) if rebuild else None
        if pattern is None and FILTERED in str(event.get('embedded_bursts', '')):
            raise SystemExit(
                f'realisation {target.picked} had its embedded bursts filtered in the run that '
                f'produced it, and the pattern read from the database is the one from before '
                f'the filter. Pass --sims-config (the upstream sims config) so the filtered '
                f'pattern is rebuilt; without it the regional storm carries the same depth in '
                f'a sharper burst than the event had.')
        records.append(writer.write(event, duration, name, gwl=gwl,
                                    rain_lag_hours=args.rain_lag, pattern=pattern))
        print(f'  wrote {name}')
        # The batch file carries every dam's antecedent level beside the URBS command
        # line, so it is always written: it is the record of what the run was, whether
        # URBS is run now or by hand afterwards.
        result_name = os.path.splitext(name)[0]
        writer.run(event, name, result_name, execute=args.run_urbs)

    table = pd.DataFrame(records)
    out = os.path.splitext(args.selection)[0] + '_downstream_storms.csv'
    table.to_csv(out, index=False)
    print(f'\n{len(records)} storm files; what went into them: {out}')
    if not args.dry_run:
        flagged = table[(table.duration_dip != '') |
                        (~table.embedded_bursts.isin(['', 'No embedded bursts']))]
        if len(flagged):
            print(f'{len(flagged)} need a look - see duration_dip and embedded_bursts')


if __name__ == '__main__':
    main()
