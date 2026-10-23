#!/usr/bin/env python3
"""One-command SSH-independent, fail-closed CRIT-v7 fixed-pilot launcher.

Usage (from the pinned standalone Git worktree):
  python3 scripts/crit_v7_launch_detached.py \
      /opt/2026/RLR/relative-logistic-regression

No network access, no daemon, no shell heredoc, no results overwrite.
All child stdin/stdout/stderr disconnected from SSH; one immutable log.
The child itself enforces exact SHA for main/original/archive.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BRANCH = "feature/crit-v7-m10-synthetic-repair-01"
REF_SHA = "4873878f57f806bf5b52144756f8087b17fa540c"
REL_RESULTS = "results/crit_v7_original_rlr_2x3_v1"
RUNNER = "scripts/crit_v7_2x3_guarded_batch.sh"


def checked_git(path: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python3 scripts/crit_v7_launch_detached.py "
              "/absolute/path/to/published/original/RLR", file=sys.stderr)
        return 2
    reference = Path(sys.argv[1]).resolve(strict=True)
    head = checked_git(REPO, "rev-parse", "HEAD")
    expected = checked_git(REPO, "rev-parse", "origin/" + BRANCH)
    if head != expected:
        print("STOP: checkout differs from the verified remote branch",
              file=sys.stderr)
        return 2
    if checked_git(reference, "rev-parse", "HEAD") != REF_SHA:
        print("STOP: original RLR source SHA mismatch", file=sys.stderr)
        return 2
    for path in [
        REPO / REL_RESULTS,
        REPO / "crit_v7_background.log",
        REPO / "pilot.pid",
    ]:
        if path.exists() or path.is_symlink():
            print(f"STOP: existing state detected, refusing overwrite: {path}",
                  file=sys.stderr)
            return 2

    log_path = REPO / "crit_v7_background.log"
    pid_path = REPO / "pilot.pid"
    with log_path.open("xb") as logfile:
        child = subprocess.Popen(
            ["nice", "-n", "10", "bash", RUNNER, str(reference)],
            cwd=REPO,
            stdin=subprocess.DEVNULL,
            stdout=logfile,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )
        pid_path.write_text(str(child.pid) + "\n", encoding="ascii")
    time.sleep(3)
    if child.poll() is not None:
        print("CRIT7_LAUNCH_EXITED_EARLY=" + str(child.returncode))
        print(log_path.read_text(encoding="utf-8", errors="replace")[-5000:])
        return 1
    print("CRIT7_DETACHED_STARTED=PASS")
    print("PID=" + str(child.pid))
    print("CODE_SHA=" + head)
    print("LOG=" + str(log_path))
    print("Closing SSH will not send terminal HUP to this process.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
