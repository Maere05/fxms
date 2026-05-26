from __future__ import annotations

import csv
import datetime as dt
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


CASE_NAMES = ("all", "naca", "ahmed", "skijumper")


@dataclass
class RunSpec:
    case: str = "skijumper"
    memory_mb: int = 512
    gpu_ids: Sequence[str] = field(default_factory=lambda: ("0",))
    label: str = "debug"
    location: str = "local"
    notes: str = ""

    def validate(self) -> None:
        if self.case not in CASE_NAMES:
            raise ValueError(f"Unknown case {self.case!r}; expected one of {CASE_NAMES}")
        if self.memory_mb <= 0:
            raise ValueError("memory_mb must be greater than 0")


@dataclass
class RemoteConfig:
    host: str
    user: str
    repo_path: str
    ssh_port: int = 22
    executable: str = "bin/FluidX3D"

    @property
    def target(self) -> str:
        return f"{self.user}@{self.host}"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def local_executable(root: Path | None = None) -> Path:
    root = repo_root() if root is None else Path(root)
    exe = root / "bin" / "FluidX3D.exe"
    return exe if exe.exists() else root / "bin" / "FluidX3D"


def export_dir(root: Path | None = None) -> Path:
    root = repo_root() if root is None else Path(root)
    return root / "bin" / "export" / "force_validation"


def run_log_path(root: Path | None = None) -> Path:
    root = repo_root() if root is None else Path(root)
    return root / "py" / "run_log.csv"


def command_for_local_run(spec: RunSpec, root: Path | None = None) -> list[str]:
    spec.validate()
    exe = local_executable(root)
    return [str(exe), spec.case, str(spec.memory_mb), *[str(gpu) for gpu in spec.gpu_ids]]


def command_for_remote_run(spec: RunSpec, remote: RemoteConfig) -> str:
    spec.validate()
    args = " ".join(shlex.quote(str(x)) for x in [spec.case, spec.memory_mb, *spec.gpu_ids])
    return f"cd {shlex.quote(remote.repo_path)} && {shlex.quote(remote.executable)} {args}"


def run_local(spec: RunSpec, root: Path | None = None, dry_run: bool = True) -> subprocess.CompletedProcess[str] | list[str]:
    cmd = command_for_local_run(spec, root)
    if dry_run:
        return cmd
    result = subprocess.run(cmd, cwd=repo_root() if root is None else Path(root), text=True, capture_output=True)
    append_run_log(spec, status="ok" if result.returncode == 0 else "failed", command=" ".join(cmd), root=root)
    return result


def run_remote(spec: RunSpec, remote: RemoteConfig, dry_run: bool = True) -> subprocess.CompletedProcess[str] | list[str]:
    remote_cmd = command_for_remote_run(spec, remote)
    cmd = ["ssh", "-p", str(remote.ssh_port), remote.target, remote_cmd]
    if dry_run:
        return cmd
    result = subprocess.run(cmd, text=True, capture_output=True)
    append_run_log(spec, status="ok" if result.returncode == 0 else "failed", command=" ".join(cmd))
    return result


def copy_remote_csvs(remote: RemoteConfig, destination: Path | None = None, dry_run: bool = True) -> subprocess.CompletedProcess[str] | list[str]:
    destination = export_dir() if destination is None else Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    source = f"{remote.target}:{remote.repo_path}/bin/export/force_validation/*.csv"
    cmd = ["scp", "-P", str(remote.ssh_port), source, str(destination)]
    if dry_run:
        return cmd
    return subprocess.run(cmd, text=True, capture_output=True)


def append_run_log(spec: RunSpec, status: str, command: str, root: Path | None = None) -> None:
    path = run_log_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["timestamp", "label", "location", "case", "memory_mb", "gpu_ids", "status", "notes", "command"])
        if not exists:
            writer.writeheader()
        writer.writerow(
            {
                "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                "label": spec.label,
                "location": spec.location,
                "case": spec.case,
                "memory_mb": spec.memory_mb,
                "gpu_ids": " ".join(str(g) for g in spec.gpu_ids),
                "status": status,
                "notes": spec.notes,
                "command": command,
            }
        )


def csv_files(root: Path | None = None) -> list[Path]:
    path = export_dir(root)
    return sorted(path.glob("*.csv")) if path.exists() else []


def load_force_csvs(root: Path | None = None):
    import pandas as pd

    frames = []
    for path in csv_files(root):
        frame = pd.read_csv(path)
        frame["source_file"] = path.name
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_run_log(root: Path | None = None):
    import pandas as pd

    path = run_log_path(root)
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def summarize_runs(df):
    if df.empty:
        return df
    columns = ["Fx_N", "Fy_drag_N", "Fz_lift_N", "Cd", "Cl", "Cs"]
    return (
        df.groupby(["case", "memory_mb", "source_file"], dropna=False)[columns]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
    )


def latest_samples(df):
    if df.empty:
        return df
    return df.sort_values("step").groupby(["case", "memory_mb", "source_file"], as_index=False).tail(1)
