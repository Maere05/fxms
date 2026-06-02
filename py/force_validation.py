from __future__ import annotations

import csv
import datetime as dt
import math
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


CASE_NAMES = ("all", "naca", "ahmed", "skijumper")
RESEARCH_MEMORY_TIERS_MB = (20_000, 30_000, 40_000)
TARGET_BANDS = {
    "ahmed": {"Cd": (0.285, 0.298)},
    "naca": {"Cd": (0.006, 0.008)},
    "skijumper": {"Fy_drag_N": (40.0, 45.0), "Fz_lift_N": (28.0, 32.0)},
}


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


def command_for_local_build(root: Path | None = None) -> list[str]:
    root = repo_root() if root is None else Path(root)
    if os.name == "nt":
        msbuild = shutil.which("MSBuild.exe") or "MSBuild.exe"
        return [msbuild, str(root / "FluidX3D.sln"), "/p:Configuration=Release", "/p:Platform=x64"]
    return ["make", "Linux"]


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


def build_local(root: Path | None = None, dry_run: bool = True) -> subprocess.CompletedProcess[str] | list[str]:
    root = repo_root() if root is None else Path(root)
    cmd = command_for_local_build(root)
    if dry_run:
        return cmd
    return subprocess.run(cmd, cwd=root, text=True, capture_output=True)


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


def target_band_checks(df):
    import pandas as pd

    if df.empty:
        return pd.DataFrame()
    rows = []
    latest = latest_samples(df)
    for _, row in latest.iterrows():
        for metric, (target_min, target_max) in TARGET_BANDS.get(row["case"], {}).items():
            value = row.get(metric)
            rows.append(
                {
                    "case": row["case"],
                    "memory_mb": row["memory_mb"],
                    "source_file": row["source_file"],
                    "metric": metric,
                    "value": value,
                    "target_min": target_min,
                    "target_max": target_max,
                    "in_band": target_min <= value <= target_max if pd.notna(value) else False,
                }
            )
    return pd.DataFrame(rows)


def grid_independence_checks(df, low_memory_mb: int = 20_000, high_memory_mb: int = 40_000, tolerance: float = 0.05):
    import pandas as pd

    if df.empty:
        return pd.DataFrame()
    rows = []
    metrics = ("Cd", "Cl")
    means = df.groupby(["case", "memory_mb"], dropna=False)[list(metrics)].mean(numeric_only=True).reset_index()
    for case in sorted(means["case"].dropna().unique()):
        low = means[(means["case"] == case) & (means["memory_mb"] == low_memory_mb)]
        high = means[(means["case"] == case) & (means["memory_mb"] == high_memory_mb)]
        if low.empty or high.empty:
            continue
        for metric in metrics:
            low_value = float(low.iloc[0][metric])
            high_value = float(high.iloc[0][metric])
            denom = max(abs(high_value), 1.0e-12)
            rel_delta = abs(low_value - high_value) / denom
            rows.append(
                {
                    "case": case,
                    "metric": metric,
                    "low_memory_mb": low_memory_mb,
                    "high_memory_mb": high_memory_mb,
                    "low_value": low_value,
                    "high_value": high_value,
                    "relative_delta": rel_delta,
                    "tolerance": tolerance,
                    "passes": rel_delta <= tolerance,
                }
            )
    return pd.DataFrame(rows)


def divergence_checks(df, growth_ratio_limit: float = 10.0):
    import pandas as pd

    if df.empty:
        return pd.DataFrame()
    metrics = ["Fx_N", "Fy_drag_N", "Fz_lift_N", "Cd", "Cl", "Cs"]
    rows = []
    for (case, memory_mb, source_file), run_df in df.sort_values("step").groupby(["case", "memory_mb", "source_file"], dropna=False):
        for metric in metrics:
            values = [float(v) for v in run_df[metric].dropna()] if metric in run_df else []
            finite = all(math.isfinite(v) for v in values)
            growth_ratio = 0.0
            if len(values) >= 2:
                first = max(abs(values[0]), 1.0e-12)
                last = abs(values[-1])
                growth_ratio = last / first
            rows.append(
                {
                    "case": case,
                    "memory_mb": memory_mb,
                    "source_file": source_file,
                    "metric": metric,
                    "finite": finite,
                    "growth_ratio": growth_ratio,
                    "growth_ratio_limit": growth_ratio_limit,
                    "passes": finite and growth_ratio <= growth_ratio_limit,
                }
            )
    return pd.DataFrame(rows)
