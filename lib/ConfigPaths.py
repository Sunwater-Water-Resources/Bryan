"""Config paths that resolve on Windows and on Linux.

Bryan's configs are written on Windows and record their relative paths that
way, with backslash separators: ``"rare_ifds": "ifd\\_ifd_files_03.json"``.
Windows treats both separators as separators, so nothing about that is visible
there. Linux does not: a backslash is an ordinary character in a filename, so
``os.path.join('storm_data', 'ifd\\_ifd_files_03.json')`` names a single file
called ``ifd\\_ifd_files_03.json`` that has never existed, rather than a file
inside a directory that has. The failure is a missing file several frames
later, which reads as a missing input rather than as a path that was never
going to resolve.

``resolve`` is the one place that difference is absorbed - every config-relative
path in Bryan goes through it, so a config authored on either platform runs on
both. It only touches separators. **Capitalisation is left alone on purpose**:
Windows resolves ``ifd`` to a directory named ``IFD`` and Linux does not, but
guessing which of several similarly-named directories was meant would turn a
config typo into a silent substitution. A config that names a directory it does
not have is wrong on Windows too - it just runs there.

Absolute paths are returned unchanged by ``os.path.join`` already, which is what
keeps ``model_exe`` (an absolute Windows path to ``urbs32.exe``) intact for the
platform that can run it.
"""
import os


def resolve(folder, relpath):
    """``relpath``, read from a config in ``folder``, as a path for this OS."""
    return os.path.normpath(os.path.join(folder, str(relpath).replace('\\', os.sep)))


def resolve_all(folder, filepaths):
    """``resolve`` every value of a config's filepaths mapping, in place."""
    for key, relpath in filepaths.items():
        filepaths[key] = resolve(folder, relpath)
    return filepaths


def regional_results_folder(config_file, config_data):
    """Where a regional model's event results are: its ``results_folder``, else its model folder.

    URBS writes its results into the folder it is run from. ``DownstreamStorms.py`` runs
    the regional model from the folder this returns, so anything reading the results back
    finds them through this too and the two cannot disagree. A config with no
    ``results_folder`` keeps the results beside the model, as they were before the key
    was read.
    """
    folder = os.path.dirname(os.path.abspath(config_file))
    return resolve(folder, config_data.get('results_folder') or config_data['model_folder'])
