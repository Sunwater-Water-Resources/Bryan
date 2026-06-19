"""
Render the Bryan manual's Markdown files to browsable HTML.

The manual is authored as plain Markdown (``*.md``) but every cross-document link
points at ``*.md.html`` (e.g. ``SubDocs/what.md.html``). That, plus the use of
Markdeep-only syntax in the source (``++underline++`` inserts, ``>`` call-outs),
means the docs are designed to be rendered with Markdeep
(https://casual-effects.com/markdeep/).

This script walks the Manual folder, and for every ``foo.md`` it writes a sibling
``foo.md.html`` Markdeep file (the original Markdown plus a Markdeep footer). Because
the output keeps the ``.md.html`` name, all of the existing relative links between
documents resolve correctly - just open ``Manual.md.html`` in a browser.

Markdeep renders client-side from a single ``markdeep.min.js``. The generated pages
reference a *local* copy first (vendored next to this script) and fall back to the
public CDN if it is missing. On first run the script tries to download the local copy
so the manual also works offline; pass ``--no-download`` to skip that.

Usage (from anywhere)::

    python Manual/render_manual.py            # render everything
    python Manual/render_manual.py --clean    # delete generated *.md.html first
    python Manual/render_manual.py --no-download

Only the Python standard library is required.
"""

import argparse
import os
import sys
import urllib.request

# Folder this script lives in - the manual root.
MANUAL_DIR = os.path.dirname(os.path.abspath(__file__))

# Local (vendored) Markdeep file and the CDN fallback.
MARKDEEP_FILENAME = 'markdeep.min.js'
MARKDEEP_LOCAL = os.path.join(MANUAL_DIR, MARKDEEP_FILENAME)
MARKDEEP_CDN = 'https://morgan3d.github.io/markdeep/latest/markdeep.min.js'

# Markdeep rendering options - tweak the look here. See casual-effects.com/markdeep.
MARKDEEP_OPTIONS = "window.markdeepOptions = {tocStyle: 'medium'};"


def ensure_markdeep(allow_download=True):
    """Make sure a local markdeep.min.js exists; try to fetch it if missing.

    Returns True if the local copy is present after the call, else False (in which
    case generated pages still work online via the CDN fallback).
    """
    if os.path.isfile(MARKDEEP_LOCAL):
        return True
    if not allow_download:
        print(f'  (no local {MARKDEEP_FILENAME}; pages will use the online CDN)')
        return False
    print(f'  Local {MARKDEEP_FILENAME} not found - attempting to download it...')
    try:
        with urllib.request.urlopen(MARKDEEP_CDN, timeout=15) as response:
            data = response.read()
        with open(MARKDEEP_LOCAL, 'wb') as f:
            f.write(data)
        print(f'  Saved {MARKDEEP_LOCAL} ({len(data) // 1024} KB) - manual will work offline.')
        return True
    except Exception as exc:  # noqa: BLE001 - any failure just means "online only"
        print(f'  Could not download Markdeep ({exc}).')
        print('  Pages will still render when online via the CDN. To enable offline')
        print(f'  viewing, place a copy of {MARKDEEP_FILENAME} next to this script.')
        return False


def find_markdown_files(root):
    """Yield every *.md file under root (excluding hidden folders)."""
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith('.')]
        for name in sorted(filenames):
            if name.lower().endswith('.md'):
                yield os.path.join(dirpath, name)


def build_footer(output_path, have_local):
    """Build the Markdeep footer, pointing at the local markdeep.min.js relatively."""
    rel_js = os.path.relpath(MARKDEEP_LOCAL, os.path.dirname(output_path))
    rel_js = rel_js.replace(os.sep, '/')
    local_script = ''
    if have_local:
        local_script = f'<script src="{rel_js}" charset="utf-8"></script>'
    return (
        '\n\n<!-- Markdeep footer (added by render_manual.py) -->\n'
        '<style class="fallback">body{visibility:hidden;white-space:pre;'
        'font-family:monospace}</style>\n'
        f'<script>{MARKDEEP_OPTIONS}</script>\n'
        f'{local_script}'
        f'<script src="{MARKDEEP_CDN}" charset="utf-8"></script>\n'
        '<script>window.alreadyProcessedMarkdeep||'
        '(document.body.style.visibility="visible")</script>\n'
    )


def render_file(md_path, have_local):
    """Render a single Markdown file to a sibling .md.html Markdeep file."""
    output_path = md_path + '.html'  # foo.md -> foo.md.html
    with open(md_path, 'r', encoding='utf-8') as f:
        markdown = f.read()
    html = '<meta charset="utf-8">\n' + markdown + build_footer(output_path, have_local)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    return output_path


def clean(root):
    """Delete previously generated *.md.html files."""
    removed = 0
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith('.')]
        for name in filenames:
            if name.lower().endswith('.md.html'):
                os.remove(os.path.join(dirpath, name))
                removed += 1
    print(f'Removed {removed} generated *.md.html file(s).')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--clean', action='store_true',
                        help='delete generated *.md.html files instead of rendering')
    parser.add_argument('--no-download', action='store_true',
                        help='do not try to download markdeep.min.js for offline use')
    args = parser.parse_args()

    if args.clean:
        clean(MANUAL_DIR)
        return

    print(f'Rendering Markdown in: {MANUAL_DIR}')
    have_local = ensure_markdeep(allow_download=not args.no_download)

    count = 0
    for md_path in find_markdown_files(MANUAL_DIR):
        output_path = render_file(md_path, have_local)
        print(f'  {os.path.relpath(md_path, MANUAL_DIR)} -> '
              f'{os.path.relpath(output_path, MANUAL_DIR)}')
        count += 1

    if count == 0:
        print('No Markdown files found.')
        sys.exit(1)

    landing = os.path.join(MANUAL_DIR, 'Manual.md.html')
    print(f'\nRendered {count} file(s). Open the manual at:\n  {landing}')


if __name__ == '__main__':
    main()
