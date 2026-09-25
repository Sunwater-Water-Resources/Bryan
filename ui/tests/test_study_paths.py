"""Paths typed or pasted into the UI: a run added to a study, a config opened.

The report that started this: Add run refused a sims_config.json that was
there. A path copied from a file's Properties > Security tab carries an
invisible U+202A in front, and a relative path is read from the study file's
folder - neither of which the message used to say.
"""

from __future__ import annotations

import json

import pytest

from core import study as studies
from core.config import load_sims_config
from core.paths import clean_path_text

LTR_EMBEDDING = "\u202a"       # what the Security tab's "Object name" puts in front


@pytest.fixture
def config(tmp_path):
    model = tmp_path / "model"
    model.mkdir()
    (model / "sims.xlsx").write_bytes(b"")
    path = model / "KRO_mc_sims_01.json"
    path.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return path


@pytest.fixture
def study(tmp_path):
    return studies.new_study(tmp_path / "report" / "bryan_study.json", "Kroombit")


@pytest.mark.parametrize("text", [
    "{path}", '"{path}"', " {path} ", LTR_EMBEDDING + "{path}",
    LTR_EMBEDDING + '"{path}"', "{path}\u200e", "\ufeff{path}"])
def test_a_pasted_path_is_added_however_it_was_copied(study, config, text):
    entry = study.add_run("E018", text.format(path=config))
    assert studies.resolve(study.folder, entry["sims_config"]) == config.resolve()


def test_clean_path_text_keeps_everything_that_is_part_of_the_path():
    assert clean_path_text(LTR_EMBEDDING + r'"C:\Models\KRO 2025\a.json" ') == \
        r"C:\Models\KRO 2025\a.json"
    assert clean_path_text(r"\\server\share\a.json") == r"\\server\share\a.json"


def test_a_relative_path_says_where_it_was_looked_for(study, config):
    with pytest.raises(studies.StudyError) as caught:
        study.add_run("E018", config.name)
    message = str(caught.value)
    assert str(study.folder / config.name) in message
    assert "read from the study file's folder" in message


def test_a_folder_is_named_as_one_with_the_configs_in_it(study, config):
    with pytest.raises(studies.StudyError, match="is a folder") as caught:
        study.add_run("E018", str(config.parent))
    assert config.name in str(caught.value)


def test_an_absolute_path_that_is_not_there_says_so_plainly(study, config):
    missing = config.with_name("nope.json")
    with pytest.raises(studies.StudyError) as caught:
        study.add_run("E018", str(missing))
    assert str(missing) in str(caught.value)
    assert "study file's folder" not in str(caught.value)


def test_the_project_page_opens_a_pasted_config_too(config):
    assert load_sims_config(LTR_EMBEDDING + str(config)).config_path == config.resolve()
