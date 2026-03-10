import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional


def find_value(obj: Any, key: str) -> Optional[Any]:
    """Recursively find the first occurrence of `key` in a JSON-like object."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for value in obj.values():
            result = find_value(value, key)
            if result is not None:
                return result
    elif isinstance(obj, list):
        for item in obj:
            result = find_value(item, key)
            if result is not None:
                return result
    return None


def resolve_path(raw_path: str, app_dir: Path, conf_dir: Path) -> Path:
    """
    Resolve a path from conf.json.

    Tries, in order:
    1. Absolute path as-is
    2. Relative to app/
    3. Relative to conf/

    If none exist yet, returns the app-relative path.
    """
    path = Path(raw_path).expanduser()

    if path.is_absolute():
        return path

    candidates = [
        app_dir / path,
        conf_dir / path,
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return (app_dir / path).resolve()


def main() -> int:
    app_dir = Path(__file__).resolve().parent
    conf_dir = app_dir / "conf"
    config_path = conf_dir / "conf.json"

    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        config = json.load(f)

    freq_file_raw = find_value(config, "freq_file")
    pops_count_file_raw = find_value(config, "pops_count_file")
    graph_files_path_raw = find_value(config, "graph_files_path")

    missing_keys = [
        key
        for key, value in {
            "freq_file": freq_file_raw,
            "pops_count_file": pops_count_file_raw,
            "graph_files_path": graph_files_path_raw,
        }.items()
        if value is None
    ]
    if missing_keys:
        raise KeyError(f"Missing keys in conf.json: {', '.join(missing_keys)}")

    if not all(isinstance(v, str) for v in [freq_file_raw, pops_count_file_raw, graph_files_path_raw]):
        raise TypeError("freq_file, pops_count_file, and graph_files_path must all be strings")

    freq_file = resolve_path(freq_file_raw, app_dir, conf_dir)
    pops_count_file = resolve_path(pops_count_file_raw, app_dir, conf_dir)
    graph_files_path = resolve_path(graph_files_path_raw, app_dir, conf_dir)

    have_needed_inputs = (
        freq_file.is_file()
        and pops_count_file.is_file()
        and (graph_files_path.is_dir() and any(graph_files_path.iterdir()))
    )

    if not have_needed_inputs:
        subprocess.run([sys.executable, "produce_example_graph_file.py"],  cwd=app_dir, check=True)
    subprocess.run([sys.executable, "app.py"], cwd=app_dir, check=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}: {e.cmd}", file=sys.stderr)
        raise SystemExit(e.returncode)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        raise SystemExit(1)