"""
File removal that survives a file being briefly held open by something else.

A Monte Carlo run writes and then deletes tens of thousands of small files. On Windows a
file that was written moments ago can still be held by an antivirus scanner, the search
indexer, or a drive sync agent, and ``os.remove`` then raises

    PermissionError: [WinError 32] The process cannot access the file because it is
    being used by another process

That lock is transient - it clears in well under a second - but a bare ``os.remove`` turns
it into a failed simulation, and for the mop-up that means losing a run over a file being
deleted purely for tidiness, after its results have already been read.

Note this is not about the hydrologic model still writing. URBS and RORB are both launched
synchronously and waited on, and a process that has exited holds no handles; ``run_storm``
also blocks on ``wait_for_results`` until the output files can be opened for append. By the
time anything here is called, the model is long gone.
"""
import os
import shutil
import time

ATTEMPTS = 5
INITIAL_DELAY = 0.2


def remove_file(path, attempts=ATTEMPTS, delay=INITIAL_DELAY):
    """Delete a file, retrying with backoff while it is locked.

    Returns ``(removed, reason)``. ``removed`` is True if the file is gone - including
    when it was never there. ``reason`` is the last error message when it is not.

    The caller decides what a failure means. For mop-up it is a warning; nothing depends
    on the file having gone.
    """
    if not os.path.exists(path):
        return True, None
    for attempt in range(attempts):
        try:
            os.remove(path)
            return True, None
        except FileNotFoundError:
            return True, None                       # something else removed it - fine
        except OSError as error:
            if attempt == attempts - 1:
                return False, str(error)
            time.sleep(delay * 2 ** attempt)
    return False, 'retries exhausted'


def remove_tree(path, attempts=ATTEMPTS, delay=INITIAL_DELAY):
    """Delete a directory tree, retrying with backoff, then raising.

    Unlike the mop-up, the callers of this one cannot carry on without it: they are
    clearing a working folder before a run, and a stale folder would leave results from a
    previous simulation in place. So the last failure is raised rather than reported.
    """
    if not os.path.exists(path):
        return
    for attempt in range(attempts):
        try:
            shutil.rmtree(path)
            return
        except FileNotFoundError:
            return
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(delay * 2 ** attempt)


class MopWarnings:
    """Counts and reports files the mop-up could not delete, without drowning the log.

    A simulation that cannot delete anything - a read-only folder, say - would otherwise
    print one line per realisation. Report the first few in full, then say so once and
    keep counting.
    """
    LIMIT = 20

    def __init__(self):
        self.count = 0

    def report(self, path, reason):
        self.count += 1
        if self.count < self.LIMIT:
            print(f'WARNING: could not delete {os.path.basename(path)} - {reason}')
        elif self.count == self.LIMIT:
            print(f'WARNING: could not delete {os.path.basename(path)} - {reason}')
            print(f'WARNING: {self.LIMIT} files could not be deleted during mop-up. '
                  f'Further mop-up warnings are suppressed for this simulation.')
            print('         The results are unaffected - these files are deleted only to '
                  'keep the folder tidy.')
