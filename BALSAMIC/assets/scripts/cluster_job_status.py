#!/usr/bin/env python3

import subprocess
import sys

jobid = sys.argv[1]

try:
    output = subprocess.check_output(
        ["/usr/bin/sacct", "-j", jobid, "--format=State", "--noheader"],
        stderr=subprocess.DEVNULL,
    )
    state = output.decode().split()[0].strip()
except subprocess.CalledProcessError:
    # Job not found, maybe it was purged
    print("unknown")
    sys.exit(0)

# Normalize states for Snakemake
if "COMPLETED" in state:
    print("success")
elif any(
    x in state for x in ["FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY"]
):
    print("failed")
else:
    print("running")
