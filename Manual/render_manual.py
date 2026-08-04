"""
Render the Bryan manual's Markdown files to browsable HTML.

The manual is authored as plain Markdown (``*.md``) with ordinary cross-document
links (e.g. ``[Getting started](SubDocs/getting_started.md)``) so the docs browse
correctly on GitHub. The use of Markdeep-only syntax in the source
(``++underline++`` inserts, ``>`` call-outs) means the docs are designed to also be
rendered with Markdeep (https://casual-effects.com/markdeep/).

This script walks the Manual folder, and for every ``foo.md`` it writes a sibling
``foo.md.html`` Markdeep file (the original Markdown plus a Markdeep footer). Along
the way it rewrites relative links that point at another ``*.md`` file in the manual
to point at that file's rendered ``*.md.html`` sibling instead, so the generated
pages link to each other correctly - just open ``Manual.md.html`` in a browser.

**Markdeep is not the same dialect the pages are drafted in.** The manual is written
in GitHub-flavoured Markdown - that is what Dillinger, where the pages are drafted,
and GitHub both render. Markdeep accepts less, and where it disagrees it does not
complain, it quietly produces something wrong. So the source is left in the dialect
it is authored in and this script translates on the way out:

* ```inline code``` (any run of backticks) becomes `inline code`. Markdeep only
  understands single backticks, and the leftovers swallow table cell separators.
* A backslash-escaped pipe in a table cell becomes ``&#124;``. Markdeep splits the
  cell on an escaped pipe anyway.
* A blank line is inserted between a table and any prose butted against it, which
  Markdeep would otherwise pull into the table as a row.

Every one of those faults shows up the same way - a row that disagrees with the rest
of its table about how many cells it holds - so after translating, ``check_tables``
checks exactly that and the script exits non-zero if anything is wrong. Add new pages
in whatever Markdown your editor produces; if it renders something Markdeep cannot
take, this will say so rather than publishing a mangled table.

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
import re
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

# Matches Markdown links: [text](target) - target is captured on its own so it can
# be rewritten without disturbing the link text.
MD_LINK_RE = re.compile(r'(\[[^\]]*\]\()([^)]+)(\))')
# Absolute URLs (http:, mailto:, etc.) - these should be left untouched.
EXTERNAL_SCHEME_RE = re.compile(r'^[a-zA-Z][a-zA-Z0-9+.-]*:')

# An inline code span written with two or more backticks, e.g. ```Include```.
# The delimiters must be the same length at both ends and the content must not
# itself contain a backtick.
MULTI_BACKTICK_RE = re.compile(r'(?<!`)(`{2,})([^`\n]+?)\1(?!`)')
# The start (or end) of a fenced code block, which must be left exactly as it is.
# Only the backtick run and an optional info string may be on the line: a line that
# merely *begins* with an inline span - "```Vf``` and ```Vc``` are..." - opens no
# fence, and treating it as one would skip the rest of that paragraph.
FENCE_RE = re.compile(r'^\s*(`{3,}|~{3,})[^`~]*$')
# A table's header separator row: pipes, dashes, colons and spaces only, holding at
# least one dash and at least one pipe. The leading pipe is optional, and requiring
# a pipe is what keeps a "---" horizontal rule from being read as a table.
TABLE_SEPARATOR_RE = re.compile(r'^\s*\|?[\s:|-]*-[\s:|-]*$')


def _is_table_separator(line):
    return bool(TABLE_SEPARATOR_RE.match(line)) and '|' in line and '-' in line


def _rewrite_link_target(target):
    """Point a relative *.md link at its rendered *.md.html sibling."""
    if EXTERNAL_SCHEME_RE.match(target):
        return target
    path, sep, fragment = target.partition('#')
    if path.lower().endswith('.md'):
        path += '.html'
    return path + sep + fragment


def rewrite_markdown_links(markdown):
    """Rewrite relative *.md links in Markdown source to the *.md.html files this
    script generates, so the rendered pages link to each other correctly."""
    return MD_LINK_RE.sub(
        lambda m: m.group(1) + _rewrite_link_target(m.group(2)) + m.group(3),
        markdown,
    )


def collapse_inline_code_fences(markdown):
    """Rewrite ```inline code``` to `inline code` for Markdeep.

    The manual writes inline code with three backticks throughout. CommonMark
    allows that - a code span may be delimited by any number of backticks - so it
    renders correctly on GitHub. **Markdeep only understands single backticks.**
    It opens a code span on the first backtick of the run and closes it on the
    second, so ```Include``` comes out as a mangled span with two stray backticks
    trailing it.

    In ordinary prose that is merely ugly. In a table it is destructive: the stray
    backtick opens a *second* code span that runs on to the next backtick further
    down the line, swallowing the ``|`` cell separators on the way. Markdeep then
    sees one enormous cell instead of a row of them, and the whole table collapses
    into a single row.

    Fenced code blocks are left alone - only spans that open and close on the one
    line are touched.

    Returns (markdown, converted, skipped): counts for reporting, where `skipped`
    is spans holding a backtick, which is the one case a single-backtick span
    cannot express.
    """
    lines = markdown.split('\n')
    in_fence = False
    converted = 0
    skipped = 0
    for i, line in enumerate(lines):
        if FENCE_RE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        line, count = MULTI_BACKTICK_RE.subn(lambda m: '`' + m.group(2) + '`', line)
        converted += count
        # Anything left with a run of backticks holds a backtick of its own, so it
        # cannot be written as a single-backtick span. Rare, but say so rather than
        # rendering it wrongly and quietly.
        skipped += len(re.findall(r'(?<!`)`{2,}(?!`)', line))
        lines[i] = line
    return '\n'.join(lines), converted, skipped


def _mark_table_lines(lines):
    """Flag which lines belong to a table, anchored on the ``|---|---|`` separator.

    Not simply "starts with a pipe": the leading pipe is optional in Markdown and
    the manual does leave it off in places. Anchoring on the separator row and
    growing outwards over the lines that hold pipes finds the whole block however
    it is punctuated, and - importantly - keeps a header row joined to its
    separator, which a narrower test would split apart.
    """
    in_fence = False
    code = [False] * len(lines)
    for i, line in enumerate(lines):
        if FENCE_RE.match(line):
            in_fence = not in_fence
            code[i] = True
            continue
        code[i] = in_fence

    is_table = [False] * len(lines)
    for i, line in enumerate(lines):
        if code[i] or not _is_table_separator(line):
            continue
        is_table[i] = True
        for step in (-1, 1):  # grow up, then down, over the rows either side
            j = i + step
            while 0 <= j < len(lines) and not code[j] and '|' in lines[j] and lines[j].strip():
                is_table[j] = True
                j += step
    return is_table


def escape_table_pipes(markdown):
    """Turn an escaped pipe inside a table cell into an HTML entity.

    ``\\|`` is how Markdown escapes a literal pipe in a table cell, and GitHub
    honours it. Markdeep does not - it splits the cell there anyway, giving that
    row one column more than the rest of the table. ``&#124;`` renders as a pipe
    in both, and is not a separator in either.

    Returns (markdown, escaped).
    """
    lines = markdown.split('\n')
    is_table = _mark_table_lines(lines)
    escaped = 0
    for i, line in enumerate(lines):
        if not is_table[i] or '\\|' not in line:
            continue
        escaped += line.count('\\|')
        lines[i] = line.replace('\\|', '&#124;')
    return '\n'.join(lines), escaped


def separate_tables(markdown):
    """Put a blank line between a table and any prose butted up against it.

    Markdeep needs a table to be a block of its own. Where a line of prose sits
    directly above or below the rows with no blank line between - which is how
    most of the manual is written, because GitHub does not care - Markdeep pulls
    that prose into the table as an extra row.

    The manual's own table captions get away with it only by accident: they end
    with a trailing space, which is a Markdown line break, and that is enough to
    start the table cleanly. A caption without one is swallowed.

    A Markdeep caption below a table ("[Caption]") is left attached, as the blank
    line would detach it.

    Returns (markdown, inserted).
    """
    lines = markdown.split('\n')
    is_table = _mark_table_lines(lines)

    out = []
    inserted = 0
    for i, line in enumerate(lines):
        starts_table = is_table[i] and not (i and is_table[i - 1])
        ends_table = is_table[i] and not (i + 1 < len(lines) and is_table[i + 1])
        if starts_table and out and out[-1].strip():
            out.append('')
            inserted += 1
        out.append(line)
        if ends_table and i + 1 < len(lines):
            following = lines[i + 1].strip()
            if following and not following.startswith('['):
                out.append('')
                inserted += 1
    return '\n'.join(out), inserted


def _count_cells(line):
    """How many cells Markdeep will split a table row into.

    The outer pipes only delimit the row, so they are dropped before splitting -
    the count is the number of gaps between the inner ones. This mirrors what
    Markdeep does, deliberately including its blind spots: a pipe inside a code
    span still splits a cell, which is why escape_table_pipes() runs first.
    """
    row = line.strip()
    if row.startswith('|'):
        row = row[1:]
    if row.endswith('|'):
        row = row[:-1]
    return len(row.split('|'))


def check_tables(markdown, path):
    """Check every table renders as a table. Returns a list of problem strings.

    Run on the *transformed* Markdown - what Markdeep is actually handed - so it
    sees what will be published rather than what was written. A table whose rows
    disagree on how many cells they hold has not been understood: something has
    swallowed a separator or introduced one. Every rendering fault found in this
    manual so far showed up exactly that way, which is what makes one cheap check
    worth having.
    """
    lines = markdown.split('\n')
    is_table = _mark_table_lines(lines)
    problems = []
    block = []
    for i, line in enumerate(lines):
        if is_table[i]:
            block.append((i, line))
            continue
        if block:
            problems += _check_table_block(block, path)
            block = []
    if block:
        problems += _check_table_block(block, path)
    return problems


def _check_table_block(block, path):
    problems = []
    widths = [(i, line, _count_cells(line)) for i, line in block]
    counts = {width for _, _, width in widths}

    if len(counts) > 1:
        commonest = max(counts, key=lambda w: sum(1 for _, _, x in widths if x == w))
        for i, line, width in widths:
            if width != commonest:
                problems.append(
                    f'{path}:{i + 1}: row has {width} cell(s), the rest of the table '
                    f'has {commonest}\n      {line.strip()[:110]}')

    # The separator belongs on the second row. Anywhere else means the rows above
    # it are prose that has been pulled in, or the header has been split off.
    separators = [i for i, line in block if _is_table_separator(line)]
    if not separators:
        problems.append(f'{path}:{block[0][0] + 1}: table block with no |---|---| separator row')
    elif len(block) > 1 and separators[0] != block[1][0]:
        problems.append(
            f'{path}:{separators[0] + 1}: the |---|---| separator is not the second row '
            f'of its table - the rows above it may be prose that has been pulled in')
    return problems


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
    markdown = rewrite_markdown_links(markdown)
    markdown, converted, skipped = collapse_inline_code_fences(markdown)
    markdown, escaped = escape_table_pipes(markdown)
    markdown, separated = separate_tables(markdown)
    problems = check_tables(markdown, os.path.relpath(md_path, MANUAL_DIR))
    html = '<meta charset="utf-8">\n' + markdown + build_footer(output_path, have_local)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    return output_path, converted, skipped, separated, escaped, problems


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
    total_converted = 0
    total_skipped = 0
    total_separated = 0
    table_problems = []
    for md_path in find_markdown_files(MANUAL_DIR):
        output_path, converted, skipped, separated, escaped, problems = render_file(md_path, have_local)
        fixes = []
        if converted:
            fixes.append(f'{converted} inline code span(s)')
        if separated:
            fixes.append(f'{separated} blank line(s) around tables')
        if escaped:
            fixes.append(f'{escaped} escaped pipe(s) in table cells')
        note = f'  [{", ".join(fixes)}]' if fixes else ''
        if skipped:
            note += f'  WARNING: {skipped} span(s) hold a backtick and were left as they are'
        print(f'  {os.path.relpath(md_path, MANUAL_DIR)} -> '
              f'{os.path.relpath(output_path, MANUAL_DIR)}{note}')
        count += 1
        total_converted += converted
        total_skipped += skipped
        total_separated += separated
        table_problems += problems

    if count == 0:
        print('No Markdown files found.')
        sys.exit(1)

    landing = os.path.join(MANUAL_DIR, 'Manual.md.html')
    print(f'\nRendered {count} file(s), converting {total_converted} multi-backtick '
          f'inline code span(s) that Markdeep does not understand and inserting '
          f'{total_separated} blank line(s) to keep tables clear of the surrounding text.')
    if total_skipped:
        print(f'{total_skipped} span(s) contain a backtick and could not be converted - '
              f'check how they render.')
    if table_problems:
        print(f'\n{"!" * 78}')
        print(f'{len(table_problems)} table problem(s) - these will not render as tables:')
        for problem in table_problems:
            print(f'  {problem}')
        print('A row that disagrees with the rest of its table on how many cells it holds')
        print('has not been understood by Markdeep. Usual causes: a pipe inside a cell that')
        print('needs escaping, or prose sitting against the table with no blank line.')
        print(f'{"!" * 78}')
        sys.exit(1)

    print('All tables check out.')
    print(f'Open the manual at:\n  {landing}')


if __name__ == '__main__':
    main()
