

def test_a_study_file_is_told_apart_from_a_sims_config(tmp_path):
    """The launcher opens either from its command line; a batch file beside the study
    naming bryan_study.json failed as a sims config with no 'simulation_list'."""
    from core import study as studies
    created = studies.new_study(tmp_path / "bryan_study.json", "Callide")
    created.save()
    sims = tmp_path / "sims_config.json"
    sims.write_text('{"simulation_list": "sims.xlsx", "filepaths": {}}')
    assert studies.is_study_file(tmp_path / "bryan_study.json")
    assert not studies.is_study_file(sims)
    assert not studies.is_study_file(tmp_path / "missing.json")
