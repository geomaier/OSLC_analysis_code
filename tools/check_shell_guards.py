#!/usr/bin/env python3
"""Check shell syntax and output guards using temporary files, without Slurm jobs."""
from pathlib import Path
import json
import os
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent.parent
SHELL_DIR = ROOT / "upstream_provenance" / "read_processing"
checks = []


def record(name, passed):
    checks.append({"check": name, "passed": bool(passed)})
    if not passed:
        raise AssertionError(name)


for source in sorted(SHELL_DIR.glob("*.sh")):
    result = subprocess.run(["bash", "-n", str(source)], capture_output=True, text=True)
    record(f"Syntax: {source.name}", result.returncode == 0)

launchers = ["launch_check_strandedness.sh", "launch_trim_map_assign.sh"]
for launcher in launchers:
    for existing in ("work", "output"):
        with tempfile.TemporaryDirectory(prefix="oslc_shell_guard_") as td:
            folder = Path(td)
            (folder / existing).mkdir()
            sentinel = folder / existing / "existing_result.txt"
            sentinel.write_text("preserve this result\n")
            samples = folder / "samples.tsv"
            samples.write_text("id\tname\nsample1\tsample1\n")
            config = {"WORKING_DIR": folder / "work", "OUTPUT_DIR": folder / "output",
                      "SAMPLE_INFO": samples}
            (folder / "config.txt").write_text("\n".join(
                f"{key}={shlex.quote(str(value))}" for key, value in config.items()) + "\n")
            bindir = folder / "bin"
            bindir.mkdir()
            blocker = bindir / "sbatch"
            blocker.write_text("#!/bin/sh\nexit 97\n")
            blocker.chmod(0o700)
            env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ.get("PATH", ""))
            result = subprocess.run(["bash", str(SHELL_DIR / launcher)], cwd=folder,
                                    env=env, capture_output=True, text=True, timeout=10)
            record(f"Preserve existing {existing}: {launcher}",
                   result.returncode != 0 and sentinel.is_file()
                   and sentinel.read_text() == "preserve this result\n"
                   and "File exists" in result.stdout + result.stderr)

source = (SHELL_DIR / "trim_map_assign.sh").read_text()
branch = re.search(r"^if \[ \$STRANDEDNESS -eq 1 \]; then\n.*?^fi$",
                   source, flags=re.MULTILINE | re.DOTALL).group(0)
for value, expected in [("0", "NONE"), ("1", "FIRST_READ_TRANSCRIPTION_STRAND"),
                        ("2", "SECOND_READ_TRANSCRIPTION_STRAND")]:
    result = subprocess.run(["bash", "-c", branch + '\nprintf "%s" "$PIC_STRAND"'],
                            env=dict(os.environ, STRANDEDNESS=value),
                            capture_output=True, text=True)
    record(f"Picard strand {value}", result.returncode == 0 and result.stdout == expected)

for launcher in launchers:
    code = (SHELL_DIR / launcher).read_text()
    counter = next(line.strip() for line in code.splitlines() if "RUNNING=$(" in line)
    for listing, expected in [("JOBID NAME", "0"), ("JOBID NAME\n123 sample", "1")]:
        command = ('set -euo pipefail\nJOB_ID=123\n'
                   'squeue() { printf "%s\\n" "$OSLC_TEST_QUEUE"; }\n')
        command += counter + '\nprintf "%s" "$RUNNING"'
        result = subprocess.run(["bash", "-c", command],
                                env=dict(os.environ, OSLC_TEST_QUEUE=listing),
                                capture_output=True, text=True)
        record(f"Queue count {expected}: {launcher}",
               result.returncode == 0 and result.stdout == expected)

with tempfile.TemporaryDirectory(prefix="oslc_shell_header_") as td:
    folder = Path(td)
    (folder / "featureCounts_sample_summaries.txt").write_text("header\n" + "row\n" * 50000)
    output = folder / "summary.txt"
    code = (SHELL_DIR / "launch_check_strandedness.sh").read_text()
    selector = next(line for line in code.splitlines() if "sed -n '1p'" in line)
    result = subprocess.run(["bash", "-o", "pipefail", "-c", selector],
                            env=dict(os.environ, OUTPUT_DIR=str(folder), SUMMARY=str(output)),
                            capture_output=True, text=True)
    record("Header selection consumes input under pipefail",
           result.returncode == 0 and output.read_text() == "header\n")

print(json.dumps({"statistical_execution": False, "checks": checks}, indent=2))
