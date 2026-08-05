"""The on-disk record of a run: ``_ui_runs/<run_id>/run.json``.

Disk is authoritative and memory is a cache, which is what lets a run survive
the browser tab closing, the UI being restarted, or a second UI being opened on
another port. This mirrors the registry in ``swe2d_ui/runner.py:202-232``, with
one addition: a process **create time** is stored beside the pid.

That matters. A pid on its own is not proof of identity - on a long-lived
Windows box pids get reused, and reattaching to a stranger's process would let
the UI report someone else's work as this run, or kill it.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .paths import atomic_write_json, read_json

RUN_STATE_FILENAME = "run.json"

QUEUED = "queued"
RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"
CANCELLED = "cancelled"
DIED = "died"

LIVE_STATES = (QUEUED, RUNNING)
FINAL_STATES = (COMPLETED, FAILED, CANCELLED, DIED)


@dataclass
class ChunkRecord:
    index: int
    config: str
    sims_list: str
    console_log: str
    run_log: str
    rows: list = field(default_factory=list)
    names: list = field(default_factory=list)
    log_paths: list = field(default_factory=list)
    status: str = QUEUED
    pid: int | None = None
    create_time: float | None = None
    started: float | None = None
    ended: float | None = None
    returncode: int | None = None
    note: str = ""

    @property
    def is_live(self) -> bool:
        return self.status in LIVE_STATES

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ChunkRecord":
        known = {key: data.get(key) for key in cls.__dataclass_fields__}
        known["rows"] = list(known.get("rows") or [])
        known["names"] = list(known.get("names") or [])
        known["log_paths"] = list(known.get("log_paths") or [])
        return cls(**known)


@dataclass
class RunRecord:
    run_id: str
    folder: str
    project_folder: str
    source_config: str
    source_sims_list: str
    created: float = field(default_factory=time.time)
    chunks: list = field(default_factory=list)
    max_parallel: int = 1
    stop_requested: bool = False
    note: str = ""

    @property
    def path(self) -> Path:
        return Path(self.folder) / RUN_STATE_FILENAME

    @property
    def is_live(self) -> bool:
        return any(chunk.is_live for chunk in self.chunks)

    @property
    def status(self) -> str:
        states = {chunk.status for chunk in self.chunks}
        if states & {RUNNING}:
            return RUNNING
        if states & {QUEUED}:
            return QUEUED
        if states & {CANCELLED}:
            return CANCELLED
        if states & {FAILED, DIED}:
            return FAILED
        return COMPLETED if states else QUEUED

    def chunk(self, index: int) -> ChunkRecord | None:
        for chunk in self.chunks:
            if chunk.index == index:
                return chunk
        return None

    def to_dict(self) -> dict:
        data = asdict(self)
        data["chunks"] = [chunk.to_dict() for chunk in self.chunks]
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "RunRecord":
        chunks = [ChunkRecord.from_dict(entry) for entry in data.get("chunks", [])]
        known = {key: data.get(key) for key in cls.__dataclass_fields__
                 if key != "chunks"}
        known["created"] = known.get("created") or time.time()
        return cls(chunks=chunks, **known)

    def save(self) -> None:
        atomic_write_json(self.path, self.to_dict())


def load(folder) -> RunRecord | None:
    data = read_json(Path(folder) / RUN_STATE_FILENAME)
    if not isinstance(data, dict) or "run_id" not in data:
        return None
    try:
        return RunRecord.from_dict(data)
    except (TypeError, ValueError):
        return None


def discover(run_root) -> list:
    """Every run under ``_ui_runs/``, newest first."""
    root = Path(run_root)
    if not root.is_dir():
        return []
    runs = []
    for folder in root.iterdir():
        if not folder.is_dir():
            continue
        record = load(folder)
        if record is not None:
            runs.append(record)
    return sorted(runs, key=lambda record: record.created, reverse=True)


def from_run_folder(run_folder, config, max_parallel=1, names=None,
                    log_paths=None) -> RunRecord:
    """Build the initial record for a freshly written run folder."""
    names = names or {}
    log_paths = log_paths or {}
    chunks = [
        ChunkRecord(
            index=chunk.index,
            config=str(chunk.config),
            sims_list=str(chunk.sims_list),
            console_log=str(chunk.console_log),
            run_log=str(chunk.run_log),
            rows=list(chunk.rows),
            names=[names.get(row, "") for row in chunk.rows],
            log_paths=[str(log_paths.get(row, "") or "") for row in chunk.rows],
        )
        for chunk in run_folder.chunks
    ]
    return RunRecord(
        run_id=run_folder.run_id,
        folder=str(run_folder.folder),
        project_folder=str(config.project_folder),
        source_config=str(config.config_path),
        source_sims_list=str(config.sims_list_path),
        chunks=chunks,
        max_parallel=max(1, int(max_parallel)),
    )


def prune(run_root, older_than_days: float = 30.0, keep_live=True) -> list:
    """Delete finished run folders older than N days. Returns what went.

    ``_ui_runs/`` grows a folder per launch, so without this it accumulates
    indefinitely.
    """
    import shutil

    cutoff = time.time() - older_than_days * 86400
    removed = []
    for record in discover(run_root):
        if record.created > cutoff:
            continue
        if keep_live and record.is_live:
            continue
        try:
            shutil.rmtree(record.folder)
            removed.append(record.run_id)
        except OSError:
            continue
    return removed
