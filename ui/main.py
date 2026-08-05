"""Entry point: ``python ui/main.py [sims_config.json] [--port N]``.

Bryan has no packaging, so the repo root goes on sys.path before anything
imports ``lib.*``. Run this with the UI's own interpreter - it does not need
scipy or matplotlib, and it must not be the frozen bryan29 environment.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

UI_ROOT = Path(__file__).resolve().parent
BRYAN_ROOT = UI_ROOT.parent
for path in (str(UI_ROOT), str(BRYAN_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="bryan-ui",
        description="Choose and run Bryan simulations from a sims list")
    parser.add_argument("config", nargs="?", default=None,
                        help="sims_config.json to open on startup")
    parser.add_argument("--host", default=None)
    parser.add_argument("--port", type=int, default=8081)
    parser.add_argument("--native", action="store_true",
                        help="a desktop window instead of a browser tab")
    parser.add_argument("--no-show", action="store_true",
                        help="do not open a browser tab")
    parser.add_argument("--bryan-python", default=None,
                        help="interpreter to run Bryan with (saved)")
    parser.add_argument("--bryan-main", default=None,
                        help="path to Bryan's Main.py (saved)")
    parser.add_argument("--max-parallel", type=int, default=None,
                        help="how many Bryan processes to run at once "
                             "(default 1 - see ui/README.md)")
    args = parser.parse_args()

    if args.native:
        try:
            import webview  # noqa: F401  - imported only to check it is there
        except ImportError:
            parser.error("--native needs pywebview: pip install 'nicegui[native]'")

    from nicegui import app, ui

    import pages
    from state import STATE

    if args.bryan_python:
        STATE.settings.bryan_python = args.bryan_python
    if args.bryan_main:
        STATE.settings.bryan_main = args.bryan_main
    if args.max_parallel is not None:
        STATE.settings.max_parallel = max(1, args.max_parallel)
    if any((args.bryan_python, args.bryan_main, args.max_parallel is not None)):
        STATE.apply_settings()

    pages.register_all()

    if args.config:
        path = Path(args.config).expanduser()
        if not path.is_file():
            parser.error(f"config file not found: {path}")
        try:
            STATE.open_project(path)
        except Exception as exc:      # noqa: BLE001 - report, do not crash the UI
            print(f"could not open {path}: {exc}", file=sys.stderr)

    # Never leave orphaned Bryan or URBS processes behind when the UI closes.
    app.on_shutdown(STATE.manager.shutdown)

    ui.run(
        title="Bryan",
        host=args.host,
        port=args.port,
        native=args.native,
        window_size=(1500, 950) if args.native else None,
        show=not args.no_show and not args.native,
        reload=False,
        dark=None,          # follow the OS theme
    )


if __name__ in {"__main__", "__mp_main__"}:
    main()
