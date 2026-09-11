"""UrbsModel makes its storms folder whichever way it is constructed.

A simulation passes ``sub_folder`` (the sims-list 'Output file'), and that branch has
always created the folder it names. A caller that passes none - DownstreamStorms.py
builds the model that way, because it is writing storm files for a regional model
rather than running a simulation - fell through to a branch that only resolved the
path. Every writer after it assumed the folder was there, so a project whose storms
folder had never been created failed at the first storm file rather than at
construction, with no indication that a missing directory was the problem.

RorbModel has always made it in both branches. This is URBS catching up.
"""
import json
import shutil
from pathlib import Path

import pytest

pytest.importorskip('scipy')

from lib.URBSmodel import UrbsModel

EXAMPLE = Path(__file__).resolve().parents[1] / 'example_project'


@pytest.fixture
def model_folder(tmp_path):
    """A copy of the example URBS model, with no storms folder."""
    if not (EXAMPLE / 'model' / 'urbs_config.json').is_file():
        pytest.skip('example project has not been generated')
    shutil.copytree(EXAMPLE / 'model' / 'urbs', tmp_path / 'urbs')
    shutil.copy(EXAMPLE / 'model' / 'urbs_config.json', tmp_path / 'urbs_config.json')
    shutil.rmtree(tmp_path / 'urbs' / 'storms', ignore_errors=True)
    return tmp_path


def test_the_storms_folder_is_made_without_a_sub_folder(model_folder):
    model = UrbsModel(str(model_folder / 'urbs_config.json'), dam_location='JUNIPER')
    assert Path(model.storms_folder).is_dir()


def test_it_is_made_with_a_sub_folder_too(model_folder):
    model = UrbsModel(str(model_folder / 'urbs_config.json'),
                      dam_location='JUNIPER', sub_folder='a_run')
    assert Path(model.storms_folder).is_dir()
    assert Path(model.storms_folder).name == 'a_run'


def test_an_existing_folder_and_its_contents_survive(model_folder):
    """exist_ok, not a rebuild: the no-sub_folder branch must not clear the folder the
    way the sub_folder branch deliberately does."""
    storms = model_folder / 'urbs' / 'storms'
    storms.mkdir(parents=True)
    keep = storms / 'previous.024'
    keep.write_text('a storm file from an earlier run')
    UrbsModel(str(model_folder / 'urbs_config.json'), dam_location='JUNIPER')
    assert keep.is_file()


def test_the_sub_folder_branch_still_clears_its_own(model_folder):
    """The other half of the contract, so the fix above cannot be copied onto it."""
    stale = model_folder / 'urbs' / 'storms' / 'a_run' / 'stale.024'
    stale.parent.mkdir(parents=True)
    stale.write_text('left over from the last run')
    UrbsModel(str(model_folder / 'urbs_config.json'),
              dam_location='JUNIPER', sub_folder='a_run')
    assert not stale.exists()
