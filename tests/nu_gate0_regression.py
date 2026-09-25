"""Gate 0-1: the nonuniform code on an equally spaced grid vs the frozen uniform code (SPEC_nonuniform §20.1).

0-1a: fields after 1, 10, 100, 1000, 5000 steps from
      (i)  the new binary with mesh=uniform (legacy layout built internally),
      (ii) the new binary with mesh=file on an equally spaced grid file,
      against ref/fdtd3d_oblique_uniform_ref.c (run once per checkpoint with nsteps=N, dump=1).
      Cases: s-pol vacuum; p-pol with an eps=2.25 half-space (arithmetic-mean interface, CPML, TF/SF, aux line).
      Metric per component and checkpoint: max|F_new - F_ref| / max|F_ref| < 1e-12 (half-grid j = Ny excluded).
0-1b: legacy regression (tests/run_all_regression.py --fresh) with the new binary -- read from its log.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402

STEPS = nc.thresholds()["gate0"]["0-1a"]["steps"]
THR = nc.thresholds()["gate0"]["0-1a"]["max_rel"]
CASES = {
    "s_vac": (dict(pol="s"), "P_uniform_nl20"),
    "p_med": (dict(pol="p", eps2=2.25, y1=60), "P_uniform_nl20_med"),
}


def rel_diff(A, B):
    out = {}
    for c in nc.COMPS:
        a, b = A[c], B[c]
        if c in ("Ey", "Hx", "Hz"):
            a, b = a[:, :-1, :], b[:, :-1, :]
        s = np.abs(a).max()
        out[c] = float(np.abs(a - b).max() / s) if s > 0 else float(np.abs(a - b).max())
    return out


def main():
    nc.build()
    T = nc.Table("gate0_regression")
    base = os.path.join(nc.RUNS, "gate0", "parity")
    worst = {}
    for case, (kw, gname) in CASES.items():
        ref = {}
        for N in STEPS:
            out = os.path.join(base, f"{case}_ref_n{N}")
            meta, _, _ = nc.run(out, binary=nc.BIN_REF, nsteps=N, dump=1, energy_every=0, **kw)
            ref[N] = nc.load_fields(out, "fields_final.bin", meta)
        out_u = os.path.join(base, f"{case}_new_uniform")
        mu, _, _ = nc.run(out_u, nsteps=max(STEPS), dump_at=STEPS, energy_every=0, **kw)
        kwf = {k: v for k, v in kw.items() if k not in ("eps2", "y1")}
        out_f = os.path.join(base, f"{case}_new_file")
        mf, _, _ = nc.run(out_f, grid=gname, nsteps=max(STEPS), dump_at=STEPS, energy_every=0, **kwf)
        for label, out, meta in (("mesh=uniform", out_u, mu), ("mesh=file (equal spacing)", out_f, mf)):
            rows = []
            for N in STEPS:
                F = nc.load_fields(out, f"fields_n{N}.bin", meta)
                d = rel_diff(F, ref[N])
                rows.append(max(d.values()))
                worst[(case, label, N)] = d
            m = max(rows)
            T.row("0-1a", f"{case}, {label}: max over comps of max|ΔF|/max|F_ref| at n = {STEPS}", "0",
                  ", ".join(f"{r:.1e}" for r in rows), f"< {THR:g}", m < THR, grid=f"Δ=λ0/20 uniform, {gname}",
                  note="bit-identical" if m == 0.0 else "")
    # 0-1b: legacy regression log (produced by tests/run_nu_regression.py or run_all_regression.py --fresh)
    log = os.path.join(nc.RESULTS, "regression.log")
    status = None
    if os.path.exists(log):
        txt = open(log, encoding="utf-8", errors="replace").read()
        m = re.findall(r"REGRESSION: (\w+)", txt)
        status = m[-1] if m else None
    T.row("0-1b", "legacy regression (run_all_regression.py --fresh, new binary, D7)", "PASS", str(status),
          "PASS", status == "PASS", note="results/regression.log")
    T.save(extra=dict(detail={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in worst.items()}))
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
