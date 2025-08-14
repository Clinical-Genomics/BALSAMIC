#!/usr/bin/env python3
import subprocess
import sys
import time

jobid = sys.argv[1]


def run(cmd, timeout=10):
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, timeout=timeout)
        return out.decode().strip()
    except Exception:
        return ""


def normalize_and_print(state):
    s = state.upper()
    if "COMPLETED" in s:
        print("success")
    elif any(
        x in s
        for x in [
            "FAILED",
            "CANCELLED",
            "TIMEOUT",
            "NODE_FAIL",
            "OUT_OF_MEMORY",
            "PREEMPTED",
        ]
    ):
        print("failed")
    else:
        # covers PENDING, RUNNING, SUSPENDED, CONFIGURING, COMPLETING, etc.
        print("running")


# Try a few times because sacct can lag right after submission/completion.
for attempt in range(5):
    # 1) First check scheduler views (good for PENDING/RUNNING)
    # Use squeue (fast) then scontrol (more detail); either might be absent on some nodes, so we guard both.
    sq = run(["/usr/bin/squeue", "-j", jobid, "-h", "-o", "%T"])
    if sq:
        normalize_and_print(sq)
        sys.exit(0)

    sc = run(["/usr/bin/scontrol", "show", "job", "-o", jobid])
    if sc:
        # Parse JobState=<STATE> from single-line -o output
        for token in sc.split():
            if token.startswith("JobState="):
                normalize_and_print(token.split("=", 1)[1])
                sys.exit(0)

    # 2) Then check accounting (good for COMPLETED/FAILED and anything after finish)
    # Query parent job line only to avoid .batch/.extern noise; -X avoids steps unless asked.
    sa = run(["/usr/bin/sacct", "-X", "-j", jobid, "-n", "--format=State"])
    if sa:
        # pick first non-empty token from the first non-empty line
        first_line = next((l for l in sa.splitlines() if l.strip()), "")
        state = first_line.split()[0] if first_line else ""
        if state:
            normalize_and_print(state)
            sys.exit(0)

    # Nothing yet — brief backoff and try again.
    time.sleep(2**attempt * 0.5)

# Still nothing — safest is to report unknown/running so Snakemake keeps waiting.
print("running")
