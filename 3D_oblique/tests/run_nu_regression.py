"""Full regression for the nonuniform-grid work (SPEC_nonuniform §19, Stage 10).

    python3 tests/run_nu_regression.py             # stop at the first FAIL
    python3 tests/run_nu_regression.py --fresh     # re-run every simulation (C and Meep)
    python3 tests/run_nu_regression.py --continue  # keep going after a FAIL (diagnosis only; grants no release)
    python3 tests/run_nu_regression.py --skip-legacy

Order: frozen-hash check -> legacy regression (run_all_regression.py) -> lint -> unit tests -> gate 0 -> gate 1 ->
gate 2 -> gate 3 (Meep) -> stage B. Exit 0 = PASS, 3 = BLOCKED (Meep missing only), 1 = FAIL.
"""
import hashlib
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEPS = [
    ("legacy regression (D7)", "tests/run_all_regression.py"),
    ("lint: coordinates from grid files", "tests/nu_lint_coords.py"),
    ("unit tests: grid generator", "tests/test_grid_gen.py"),
    ("theory preflight (--check)", "tests/nu_predict.py"),
    ("gate 0-1 regression parity", "tests/nu_gate0_regression.py"),
    ("gate 0-2 adjointness", "tests/nu_gate0_adjoint.py"),
    ("gate 0-3/0-4 energy, stability", "tests/nu_gate0_energy_stability.py"),
    ("gate 1 vacuum (+ 0-5)", "tests/nu_gate1_vacuum.py"),
    ("gate 2 film / stack", "tests/nu_gate2_film.py"),
    ("gate 3 Meep", "tests/nu_gate3_meep.py"),
    ("stage B gate 0", "tests/nuB_gate0.py"),
    ("stage B B1/B2", "tests/nuB_gate1.py"),
    ("stage B B3 grating", "tests/nuB_gate3_grating.py"),
]


def frozen_hashes():
    """thresholds_nu.json and nu_predictions.json must match the hashes recorded in CLAUDE.md."""
    txt = open(os.path.join(ROOT, "CLAUDE.md"), encoding="utf-8").read()
    ok = True
    for rel in ("tests/thresholds_nu.json", "results/nu_predictions.json"):
        m = re.search(re.escape(rel) + r": ([0-9a-f]{64})", txt)
        h = hashlib.sha256(open(os.path.join(ROOT, rel), "rb").read()).hexdigest()
        good = m is not None and m.group(1) == h
        print(f"  {rel}: {'OK' if good else 'CHANGED WITHOUT A RECORDED DECISION'}")
        ok &= good
    return ok


def main():
    env = dict(os.environ)
    if "--fresh" in sys.argv:
        env["FDTD_FRESH"] = "1"
    subprocess.run(["make", "-s"], cwd=ROOT, check=True)
    subprocess.run(["make", "-s", "ref"], cwd=ROOT, check=True)
    print("===== frozen thresholds / predictions =====")
    summary = [("frozen hashes", "PASS" if frozen_hashes() else "FAIL", 0.0)]
    steps = STEPS[1:] if "--skip-legacy" in sys.argv else STEPS
    for name, script in steps:
        if summary[-1][1] == "FAIL" and "--continue" not in sys.argv:
            break
        t0 = time.time()
        print(f"\n===== {name}: {script} =====", flush=True)
        args = [sys.executable, os.path.join(ROOT, script)]
        if script.endswith("run_all_regression.py") and "--fresh" in sys.argv:
            args.append("--fresh")
        if script.endswith("nu_predict.py"):
            args.append("--check")
        if script.endswith("run_all_regression.py"):
            with open(os.path.join(ROOT, "results", "regression.log"), "w") as fh:
                rc = subprocess.run(args, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
            print(open(os.path.join(ROOT, "results", "regression.log")).read()[-1200:])
        else:
            rc = subprocess.run(args, cwd=ROOT, env=env).returncode
        status = "PASS" if rc == 0 else ("BLOCKED" if rc == 3 else "FAIL")
        summary.append((name, status, time.time() - t0))
    print("\n===== nonuniform regression summary =====")
    for name, status, dt in summary:
        print(f"{name:<40} {status:<8} {dt:8.0f} s")
    states = [s for _, s, _ in summary]
    ran_all = len(summary) == len(steps) + 1
    final = "PASS" if ran_all and all(s == "PASS" for s in states) else ("FAIL" if "FAIL" in states or not ran_all
                                                                         else "BLOCKED")
    if "--continue" in sys.argv and "FAIL" in states:
        print("--continue: steps after the first FAIL were run for diagnosis only; they grant no release")
    print("REGRESSION:", final)
    with open(os.path.join(ROOT, "results", "nu_regression.log"), "a", encoding="utf-8") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + " ".join(f"{n}={s}" for n, s, _ in summary) + f" -> {final}\n")
    return {"PASS": 0, "BLOCKED": 3}.get(final, 1)


if __name__ == "__main__":
    sys.exit(main())
