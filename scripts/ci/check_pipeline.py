"""Syntax-check science source and shell jobs without running processing code."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    locations = [root / part for part in ("src", "tests", "scripts", "cluster", "viz")]
    python_count = shell_count = 0
    for location in locations:
        if not location.is_dir():
            raise ValueError(f"Missing source directory: {location}")
        for path in sorted(location.rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix == ".py":
                compile(path.read_bytes(), str(path), "exec")
                python_count += 1
            elif path.suffix in {".sh", ".sbatch"}:
                subprocess.run(["bash", "-n", str(path)], check=True)
                shell_count += 1
    if not python_count or not shell_count:
        raise ValueError("Expected both Python source and cluster shell jobs")
    print(f"Syntax valid: {python_count} Python files and {shell_count} shell job files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
