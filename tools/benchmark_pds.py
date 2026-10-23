#!/usr/bin/env python3
"""Serial, bounded PDS comparison against SCIP/SoPlex; retain every solver log."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]


def match(pattern, text):
    found = re.search(pattern, text, re.MULTILINE)
    return found.group(1) if found else ""


def run(command, log, timeout, env):
    started = time.monotonic()
    timed_out = False
    with log.open("w") as output:
        process = subprocess.Popen(command, cwd=ROOT, env=env,
                                   stdout=output, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
    elapsed = time.monotonic() - started
    return process.returncode, elapsed, timed_out, log.read_text(errors="replace")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--iterations", type=int, default=300000)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/pds-benchmark/tools/test_netlib")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    env = os.environ.copy()
    env.pop("LP_SIMPLEX_PROFILE", None)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", LC_ALL="C")
    sources = sorted([*ROOT.glob("src/**/*.c"), *ROOT.glob("src/**/*.h"),
                      *ROOT.glob("include/**/*.h"), ROOT / "CMakeLists.txt",
                      ROOT / "tools/test_netlib.c"])
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sources}
    metadata = dict(platform=platform.platform(), machine=platform.machine(),
                    head=subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                 cwd=ROOT, text=True).strip(),
                    wall_limit_seconds=args.seconds, iterations=args.iterations,
                    binary=str(args.binary.resolve()), source_sha256=hashes,
                    binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                    thread_environment={k: env[k] for k in (
                        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    fields = ["instance", "solver", "status", "returncode", "wall_seconds",
              "solve_seconds", "iterations", "objective", "primal_residual",
              "dual_residual", "postsolve_primal", "reference_objective", "log"]
    with (out / "results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        handle.flush()
        for size in range(20, 101, 10):
            name = f"pds-{size}"
            model = ROOT / "data/pds/mps" / f"{name}.mps"
            reference = ""
            for solver in ("scip-soplex", "lp_simplex"):
                if solver == "scip-soplex":
                    command = ["scip"]
                    for item in ("set parallel maxnthreads 1", "set lp advanced threads 1",
                                 "set presolving milp threads 1", "set lp initalgorithm d",
                                 "set lp resolvealgorithm d",
                                 f"set limits time {args.seconds}", f"read {model}",
                                 "optimize", "display statistics", "quit"):
                        command.extend(["-c", item])
                else:
                    command = [str(args.binary.resolve()), "--algorithm", "dual-revised",
                               "--iterations", str(args.iterations)]
                    if reference:
                        command.extend(["--expected", reference])
                    command.append(str(model))
                log = out / "logs" / f"{solver}-{name}.log"
                (out / "logs" / f"{solver}-{name}.command.json").write_text(
                    json.dumps(command, indent=2) + "\n")
                # SCIP's internal time limit starts after model reading. A slightly
                # larger external limit also bounds loading/shutdown for that run.
                limit = args.seconds + 10 if solver == "scip-soplex" else args.seconds
                code, elapsed, timed_out, text = run(command, log, limit, env)
                row = dict(instance=name, solver=solver, returncode=code,
                           wall_seconds=f"{elapsed:.6f}", reference_objective=reference,
                           log=str(log.relative_to(out)))
                if solver == "scip-soplex":
                    status = match(r"^SCIP Status\s*:\s*(.+)$", text)
                    row.update(status="TIMEOUT" if timed_out else status or "ERROR",
                               solve_seconds=match(r"^Solving Time \(sec\)\s*:\s*([\d.]+)", text))
                    if "optimal solution found" in status and code == 0 and not timed_out:
                        reference = match(r"^Primal Bound\s*:\s*([-+\d.eE]+)", text)
                        row["objective"] = reference
                else:
                    row.update(status="TIMEOUT" if timed_out else
                               match(r"^result:.*\(([^)]+)\)", text) or "ERROR",
                               solve_seconds=match(r"^time:\s*([\d.]+)", text),
                               iterations=match(r"^iterations:\s*(\d+)", text),
                               objective=match(r"^objective:\s*([-+\d.eE]+)", text),
                               primal_residual=match(r"^residuals:\s*primal=([^,]+)", text),
                               dual_residual=match(r"^residuals:.*dual=(\S+)", text),
                               postsolve_primal=match(r"^postsolve:\s*primal=([^,]+)", text))
                    if not timed_out and row["status"] == "success" and code != 0:
                        row["status"] = "VALIDATION_FAILED"
                writer.writerow(row)
                handle.flush()
                print(f"{name} {solver}: {row['status']}; wall={elapsed:.3f}s; "
                      f"objective={row.get('objective', '')}", flush=True)


if __name__ == "__main__":
    main()
