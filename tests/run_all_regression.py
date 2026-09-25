"""Full regression: every validation script in order; stops at the first failure.

    python3 tests/run_all_regression.py            # reuse cached simulations whose inputs are unchanged
    python3 tests/run_all_regression.py --fresh    # force every simulation (C and Meep) to re-run
    python3 tests/run_all_regression.py --continue # keep going after a failure (diagnosis only)

Exit code 0 only if every step passes. Level 3 exits 3 (BLOCKED) if Meep is unavailable; that is reported
as BLOCKED and the regression result is BLOCKED, never PASS.
"""
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEPS = [
    ("Level 0   exact injection", "tests/level0_exact_injection.py"),
    ("Stage 1   C multistep exactness", "tests/level0b_multistep.py"),
    ("Stage 2   CPML stability", "tests/level_cpml_stability.py"),
    ("Level 1-1 leakage", "tests/level1_leakage.py"),
    ("Level 1-2..8", "tests/level1_all.py"),
    ("Level 2   Fresnel / Snell", "tests/level2_fresnel.py"),
    ("Level 3   Meep comparison", "tests/level3_meep.py"),
]


def main():
    env = dict(os.environ)
    if "--fresh" in sys.argv:
        env["FDTD_FRESH"] = "1"
    subprocess.run(["make", "-s"], cwd=ROOT, check=True)
    summary = []
    for name, script in STEPS:
        t0 = time.time()
        print(f"\n===== {name}: {script} =====", flush=True)
        rc = subprocess.run([sys.executable, os.path.join(ROOT, script)], cwd=ROOT, env=env).returncode
        status = "PASS" if rc == 0 else ("BLOCKED" if rc == 3 else "FAIL")
        summary.append((name, status, time.time() - t0))
        if rc != 0 and "--continue" not in sys.argv:
            break
    print("\n===== regression summary =====")
    for name, status, dt in summary:
        print(f"{name:<34} {status:<8} {dt:7.0f} s")
    if len(summary) < len(STEPS):
        print(f"stopped after '{summary[-1][0]}' ({summary[-1][1]}); later steps not run")
    states = [s for _, s, _ in summary]
    final = "PASS" if all(s == "PASS" for s in states) and len(summary) == len(STEPS) else \
        ("FAIL" if "FAIL" in states else "BLOCKED")
    if "--continue" in sys.argv and "FAIL" in states:
        first = states.index("FAIL")
        print(f"--continue: steps after '{summary[first][0]}' were run for diagnosis only; "
              "per SPEC a level is not entered until the previous one passes, so they grant no release")
    print("REGRESSION:", final)
    return {"PASS": 0, "BLOCKED": 3}.get(final, 1)


if __name__ == "__main__":
    sys.exit(main())
