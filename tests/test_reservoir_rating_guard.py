"""The reservoir routing method refuses a rating it cannot route.

The storage-indication solve inverts psi(S) = 2000*S/dt + O(S) with np.interp,
which needs psi to increase with storage.  A flat run of outflow keeps psi
increasing (the storage term does), so it routes correctly and is accepted.  A
discharge that FALLS as storage rises can make psi fall, and np.interp does not
raise on that -- it returns wrong storages.  So the .sq reader refuses it, and a
storage column that does not rise.

The case that prompted this: Callide rating 10.3 switched from ogee to orifice
control at 219.82 m while the orifice already passed less, so the spillway flow
fell 56 m3/s there; its tables hid the fall only by the spacing of their rows.
"""

import pytest

from lib.ReservoirRouting import _read_sq

HEADER = "TEST DAM\n* note\n* note\n* note\n{n} PAIRS:\n"


def write_sq(tmp_path, pairs):
    path = tmp_path / "rating.sq"
    body = "".join(f"{s} {q}\n" for s, q in pairs)
    path.write_text(HEADER.format(n=len(pairs)) + body)
    return path


def test_an_increasing_rating_is_read(tmp_path):
    s, o = _read_sq(write_sq(tmp_path, [(0, 0.0), (100, 50.0), (250, 300.0)]), fsv=1000.0)
    assert list(s) == [1000.0, 1100.0, 1250.0]
    assert list(o) == [0.0, 50.0, 300.0]


def test_a_flat_run_is_accepted_because_it_still_routes(tmp_path):
    # gates held shut above full supply: zero outflow over the first storages
    s, o = _read_sq(write_sq(tmp_path, [(0, 0.0), (1296, 0.0), (1633, 27.32)]), fsv=0.0)
    assert list(o) == [0.0, 0.0, 27.32]


def test_a_falling_discharge_is_refused(tmp_path):
    path = write_sq(tmp_path, [(0, 0.0), (57500, 7440.2), (57602, 7384.2), (60000, 8000.0)])
    with pytest.raises(ValueError, match="falls from 7440.2 to 7384.2"):
        _read_sq(path, fsv=0.0)


def test_a_storage_that_does_not_rise_is_refused(tmp_path):
    path = write_sq(tmp_path, [(0, 0.0), (100, 10.0), (100, 20.0)])
    with pytest.raises(ValueError, match="storage"):
        _read_sq(path, fsv=0.0)
