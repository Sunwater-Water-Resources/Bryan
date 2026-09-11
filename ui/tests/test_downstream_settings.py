"""The Downstream page's two config paths are remembered per project.

They are the only thing on that page worth keeping: a downstream storm config and a
regional model config are long absolute paths, typed by hand, and unchanged from one
visit to the next. The duration and warming level are deliberately not kept - blank
means "read it from the database name", and a stale override is worse than retyping.

Keyed on the **resolved path** of the sims config, never its name. Two projects
routinely carry the same config filename at different paths, one per model revision,
and keying on the name would hand one project's regional model to the other.
"""
from pathlib import Path

import pytest

from settings import UiSettings


@pytest.fixture
def settings(tmp_path, monkeypatch):
    monkeypatch.setattr('settings.SETTINGS_PATH', tmp_path / '.bryan_ui.json')
    return UiSettings()


def test_an_unknown_project_comes_back_blank(settings, tmp_path):
    assert settings.downstream_for(tmp_path / 'sims.json') == {"config": "", "model": ""}


def test_what_was_remembered_comes_back(settings, tmp_path):
    cfg = tmp_path / 'sims.json'
    settings.remember_downstream(cfg, config='D:/s/downstream.json', model='D:/m/urbs.json')
    assert settings.downstream_for(cfg) == {"config": 'D:/s/downstream.json',
                                            "model": 'D:/m/urbs.json'}


def test_it_survives_a_save_and_reload(settings, tmp_path):
    cfg = tmp_path / 'sims.json'
    settings.remember_downstream(cfg, config='a.json', model='b.json')
    assert UiSettings.load().downstream_for(cfg) == {"config": 'a.json', "model": 'b.json'}


def test_same_filename_at_two_paths_stays_apart(settings, tmp_path):
    """E010 and E011 name their sims config the same thing. They must not share."""
    e010 = tmp_path / 'E010' / 'sims.json'
    e011 = tmp_path / 'E011' / 'sims.json'
    e010.parent.mkdir(parents=True)
    e011.parent.mkdir(parents=True)
    settings.remember_downstream(e010, config='ds_E010.json', model='urbs_E010.json')
    settings.remember_downstream(e011, config='ds_E011.json', model='urbs_E011.json')
    assert settings.downstream_for(e010)["model"] == 'urbs_E010.json'
    assert settings.downstream_for(e011)["model"] == 'urbs_E011.json'


def test_the_key_is_resolved_so_the_same_project_matches_either_way(settings, tmp_path):
    cfg = tmp_path / 'sims.json'
    cfg.write_text('{}')
    settings.remember_downstream(cfg, config='a.json', model='b.json')
    roundabout = tmp_path / 'sub' / '..' / 'sims.json'
    (tmp_path / 'sub').mkdir()
    assert settings.downstream_for(roundabout)["config"] == 'a.json'


def test_clearing_both_fields_forgets_the_project(settings, tmp_path):
    cfg = tmp_path / 'sims.json'
    settings.remember_downstream(cfg, config='a.json', model='b.json')
    settings.remember_downstream(cfg, config='', model='')
    assert settings.downstream_for(cfg) == {"config": "", "model": ""}
    assert settings.downstream_configs == {}


def test_clearing_one_field_does_not_restore_the_old_value(settings, tmp_path):
    """A merge-on-write would make a cleared field come back on the next keystroke."""
    cfg = tmp_path / 'sims.json'
    settings.remember_downstream(cfg, config='a.json', model='b.json')
    settings.remember_downstream(cfg, config='', model='b.json')
    assert settings.downstream_for(cfg)["config"] == ''


def test_settings_written_before_the_field_existed_still_load(settings, tmp_path):
    """An older ~/.bryan_ui.json has no downstream_configs key at all."""
    import json
    from core.paths import atomic_write_json
    atomic_write_json(Path(settings.__class__.__module__ and tmp_path / '.bryan_ui.json'),
                      {"bryan_python": "/usr/bin/python3", "recent_configs": []})
    loaded = UiSettings.load()
    assert loaded.downstream_configs == {}
    assert loaded.downstream_for(tmp_path / 'sims.json') == {"config": "", "model": ""}
