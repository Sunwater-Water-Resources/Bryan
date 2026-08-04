"""
lib/ConsoleTitle.py

Puts the running simulation in the console window's title bar, so a long batch
can be followed from the taskbar without reading the scroll.

Works on both platforms Bryan gets run on, by two different mechanisms:

  Windows - the kernel32 SetConsoleTitleW API. This is what cmd.exe's own
            "title" command calls, and it works in cmd.exe, PowerShell and
            Windows Terminal alike. It does not go through stdout, so it still
            works when a batch file redirects the output to a file.
  Linux /  - the xterm OSC escape sequence "ESC ] 0 ; <text> BEL", understood by
  macOS     essentially every terminal emulator. This *does* go through stdout,
            so it is only written when stdout is a terminal - otherwise the
            escape codes would end up as rubbish in a redirected log file.

There is no read-back: Windows could be asked with GetConsoleTitleW, but the
terminals that understand the OSC sequence have no portable way to report the
current title, so the title is set to a closing message at the end of a batch
rather than restored to whatever it was.

Nothing here is allowed to break a run: every call is best-effort and swallows
its errors. A title bar is a convenience, not a result.

The escape sequence is written to sys.__stdout__ rather than sys.stdout on
purpose. The simulators replace sys.stdout with a Logger that tees to their log
file, and the escape codes would be written into those logs as well.
"""
import os
import sys

_IS_WINDOWS = os.name == 'nt'
_MAX_LENGTH = 200  # titles are truncated by the terminal well before this


def _clean(text):
    """Strip out anything that would break the escape sequence or the API call."""
    text = ' '.join(str(text).split())
    text = ''.join(ch for ch in text if ch.isprintable())
    if len(text) > _MAX_LENGTH:
        text = text[:_MAX_LENGTH - 3] + '...'
    return text


def set_console_title(text):
    """Set the console window title. Returns True if it was set.

    Never raises - a failure here means the run carries on with whatever title
    the window already had.
    """
    title = _clean(text)
    if not title:
        return False
    try:
        if _IS_WINDOWS:
            import ctypes
            return bool(ctypes.windll.kernel32.SetConsoleTitleW(title))
        stream = sys.__stdout__
        if stream is None or not stream.isatty():
            # Output is redirected to a file or a pipe - writing the escape
            # sequence there would corrupt it.
            return False
        stream.write(f'\033]0;{title}\007')
        stream.flush()
        return True
    except Exception:
        return False
