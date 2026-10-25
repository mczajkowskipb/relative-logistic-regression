#!/usr/bin/env python3
"""Launch repaired CRIT-v7 as a PERSISTENT user systemd SERVICE, not a session child.

Usage:
  python3 scripts/crit_v7_launch_systemd.py \
    /opt/2026/RLR/relative-logistic-regression

The user must authorize any permission changes for loginctl lingering.
Never overwrites the output directory, existing service, or historical logs.
The service is detached from SSH and owned by user@UID.service.
"""
from __future__ import annotations

import getpass
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REF_SHA = "4873878f57f806bf5b52144756f8087b17fa540c"
MAIN_SHA = "4079a5ac400f1d2cfded523c7000f655470925ce"
BRANCH = "feature/crit-v7-systemd-env-fix3"
SERVICE = "crit7-fix3.service"
OUT = HERE / "results/crit_v7_original_rlr_2x3_v1"
RUNNER = "scripts/crit_v7_2x3_guarded_batch.sh"


def command(*argv: str) -> str:
    r = subprocess.run(argv, text=True, capture_output=True)
    if r.returncode:
        raise RuntimeError(
            f"Command failed rc={r.returncode}: {argv!r}\n"
            f"{r.stdout}\n{r.stderr}"
        )
    return r.stdout.strip()


def git(repo: Path, *args: str) -> str:
    return command("git", "-C", str(repo), *args)


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError(
            "Usage: python3 scripts/crit_v7_launch_systemd.py "
            "/absolute/path/to/frozen/original/RLR"
        )
    ref = Path(sys.argv[1]).resolve(strict=True)
    sha = git(HERE, "rev-parse", "HEAD")
    remote = git(HERE, "rev-parse", "origin/" + BRANCH)
    if sha != remote:
        raise RuntimeError("STOP: fix2 branch HEAD differs from origin")
    if git(HERE, "rev-parse", "origin/main") != MAIN_SHA:
        raise RuntimeError("STOP: original GitHub main SHA unexpectedly changed")
    if git(ref, "rev-parse", "HEAD") != REF_SHA:
        raise RuntimeError("STOP: original RLR paper checkout SHA changed")
    if OUT.exists() or OUT.is_symlink():
        raise RuntimeError("STOP: new result directory already exists: " + str(OUT))
    if not (HERE / RUNNER).is_file():
        raise RuntimeError("STOP: runner absent")

    loaded = subprocess.run(
        ["systemctl", "--user", "show", SERVICE, "--property=LoadState",
         "--value"], text=True, capture_output=True
    )
    if loaded.returncode:
        raise RuntimeError(
            "STOP: user systemd manager unavailable:\n" + loaded.stderr
        )
    if loaded.stdout.strip() != "not-found":
        raise RuntimeError("STOP: systemd service already exists in state "
                           + loaded.stdout.strip())

    user = getpass.getuser()
    linger = command("loginctl", "show-user", user,
                     "--property=Linger", "--value")
    if linger != "yes":
        # This changes only the caller's login manager setting.
        # On systems where polkit does not permit it, this fails closed.
        print("Linger=no. Enabling persistence for user's services "
              "after SSH logout...", flush=True)
        command("loginctl", "enable-linger", user)
        if command("loginctl", "show-user", user,
                   "--property=Linger", "--value") != "yes":
            raise RuntimeError("STOP: lingering not confirmed, will not start")
    print("CRIT7_USER_LINGER=PASS", flush=True)

    runner_command = [
        "systemd-run", "--user",
        "--unit=" + SERVICE.removesuffix(".service"),
        "--description=CRIT7 original RLR 2x3 controlled pilot",
        "--service-type=exec",
        "--working-directory=" + str(HERE),
        "--nice=10",
    ]
    # The NAME-only form is documented, but fix2's environment block
    # was rejected on Ubuntu. Pass explicit NAME=VALUE to rule out
    # ambiguous environment inheritance and validate control characters.
    env_vars = {
        "PATH": os.environ.get("PATH") or "/usr/local/bin:/usr/bin:/bin",
        "HOME": os.environ.get("HOME") or str(Path.home()),
    }
    for var in ("R_LIBS_USER", "R_LIBS_SITE", "LD_LIBRARY_PATH"):
        if var in os.environ:
            env_vars[var] = os.environ[var]
    for name, value in env_vars.items():
        if any(c in value for c in ("\x00", "\r", "\n")):
            raise RuntimeError("STOP: invalid environment value for " + name)
        runner_command.append("--setenv=" + name + "=" + value)
    runner_command.extend(["/usr/bin/bash", RUNNER, str(ref)])
    result = command(*runner_command)
    print(result, flush=True)
    time.sleep(2)
    print(command("systemctl", "--user", "show", SERVICE,
                  "--property=ActiveState", "--property=SubState",
                  "--property=Result", "--property=ExecMainStatus"),
          flush=True)
    print("CRIT7_SYSTEMD_SERVICE_SUBMITTED=PASS", flush=True)
    print("JOURNAL_COMMAND=journalctl --user -u " + SERVICE +
          " -n 80 --no-pager", flush=True)
    print("STATUS_FILE=" + str(OUT / "status.txt"), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError) as e:
        print("CRIT7_SYSTEMD_LAUNCH_STOP: " + str(e), file=sys.stderr)
        sys.exit(2)
