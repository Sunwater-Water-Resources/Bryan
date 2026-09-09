"""Write regional storm files for the representative events chosen in the launcher.

    python DownstreamStorms.py --config downstream_storm_config_CLD_01.json ^
        --selection sims_mc\\results\\GWL1p3_representative_events.json ^
        --model runs\\E010\\Downstream\\Regional\\rev_2\\urbs_config_rev2.json

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
    ap.add_argument('--dry-run', action='store_true', help='report without writing')
    args = ap.parse_args()

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

    rainfall = DownstreamRainfall(args.config)
    writer = None if args.dry_run else DownstreamStormWriter(
        rainfall, storm_config, args.model, climate)

    records = []
    for target in chosen:
        duration = from_name(DURATION, target.database, args.duration, 'duration')
        gwl = from_name(GWL, target.database, args.gwl, 'warming level')
        event = DownstreamStormWriter.read_event(target.database, int(target.picked))
        name = storm_name(target, duration, gwl, args.suffix)
        print(f'\n{target.kind} {target.value:g} ({target.result_type}) -> realisation '
              f'{target.picked}, {duration:g} h, GWL {gwl:g}')
        if args.dry_run:
            records.append(dict(filename=name, duration=duration, gwl=gwl,
                                driver_aep=float(event['rain_aep']),
                                storm_method=str(event['storm_method'])))
            print(f'  would write {name}')
            continue
        records.append(writer.write(event, duration, name, gwl=gwl,
                                    rain_lag_hours=args.rain_lag))
        print(f'  wrote {name}')

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
