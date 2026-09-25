"""Level 3 driver: run the Meep reference cases (in parallel, Meep is single-threaded here) then compare.py.

Meep interpreter: $MEEP_PYTHON (default ~/micromamba/envs/mp/bin/python). Missing Meep -> exit 3 (BLOCKED).
Wave 1: vac10_{s,p} (Level 1 geometry), vac8_{s,p} (Level 2 normalization)
Wave 2: med_{s,p}_avg{1,0} (n=1.5 half-space, normalized by vac8_{pol})
Runs are reused if meep_meta.json is newer than meep_ref.py (set FDTD_FRESH=1 to force).
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MEEP_PY = os.environ.get("MEEP_PYTHON", os.path.expanduser("~/micromamba/envs/mp/bin/python"))
RUNS = os.path.join(ROOT, "runs", "meep")
SCRIPT = os.path.join(ROOT, "meep_ref.py")
FRESH = os.environ.get("FDTD_FRESH", "0") == "1"


def needed(tag):
    m = os.path.join(RUNS, tag, "meep_meta.json")
    return FRESH or not os.path.exists(m) or os.path.getmtime(m) < os.path.getmtime(SCRIPT)


def wave(jobs):
    procs = []
    for tag, args in jobs:
        if not needed(tag):
            print(f"  reuse {tag}")
            continue
        os.makedirs(os.path.join(RUNS, tag), exist_ok=True)
        log = open(os.path.join(RUNS, tag, "run.log"), "w")
        cmd = [MEEP_PY, SCRIPT, "--nl", "20", "--out", os.path.join(RUNS, tag)] + args
        print("  start", tag)
        procs.append((tag, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                            env=dict(os.environ, OMP_NUM_THREADS="1")), log))
    bad = []
    for tag, p, log in procs:
        rc = p.wait()
        log.close()
        print(f"  done {tag} (exit {rc})")
        if rc != 0:
            bad.append((tag, rc))
    return bad


def main():
    if not os.path.exists(MEEP_PY) or subprocess.run([MEEP_PY, "-c", "import meep"], capture_output=True).returncode:
        print("BLOCKED: Meep interpreter not available at", MEEP_PY)
        return 3
    only = sys.argv[1:] == ["--runs-only"]
    w1 = [(f"vac10_{p}", ["--pol", p, "--tf_l", "10"]) for p in ("s", "p")] + \
         [(f"vac8_{p}", ["--pol", p, "--tf_l", "8"]) for p in ("s", "p")]
    w2 = [(f"med_{p}_avg{a}", ["--pol", p, "--tf_l", "8", "--medium", "--eps_averaging", str(a),
                               "--norm", os.path.join(RUNS, f"vac8_{p}")]) for p in ("s", "p") for a in (1, 0)]
    for w in (w1, w2):
        bad = wave(w)
        if bad:
            print("Meep runs failed:", bad)
            return 1
    if only:
        return 0
    return subprocess.run([sys.executable, os.path.join(ROOT, "compare.py")], cwd=ROOT).returncode


if __name__ == "__main__":
    sys.exit(main())
