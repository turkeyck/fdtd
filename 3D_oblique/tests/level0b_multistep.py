"""Stage 1: the C engine propagates the exact discrete plane wave for many steps.

The whole grid is filled with the analytic wave (derivation.md §2.4), the y ends are driven
with the exact values (Dirichlet), no PML, no source. After N steps every component on every
node must still equal the analytic value to < 1e-12 relative; boundary (PBC-wrapped) cells are
reported separately. Also cross-checks the C theory numbers (ky, K~, E0, H0) against fdtd_theory.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import fdtd_io  # noqa: E402
from fdtd_theory import COMPS, Setup, analytic  # noqa: E402

TOL = 1e-12
NSTEPS = 200


def check(nl, pol, m, n):
    out = os.path.join(ROOT, "runs", "stage1", f"nl{nl}_{pol}_m{m}_n{n}")
    meta = fdtd_io.run(out, nl=nl, pol=pol, m=m, n=n, npml=0, sf=4, tf=16, inc="0", init="a",
                       nsteps=NSTEPS, dump=1, energy_every=0)
    st = Setup(n_lambda=nl, pol=pol, m=m, n=n)
    e, h = st.amplitudes()
    theory = max(abs(meta["ky"] - st.ky) / st.ky,
                 float(np.max(np.abs(np.array(meta["Ktilde"]) - st.Kt))),
                 float(np.max(np.abs(np.array(meta["E0"]) - e))),
                 float(np.max(np.abs(np.array(meta["H0"]) - h))))
    F = fdtd_io.load_fields(out, meta)
    Nx, Ny, Nz = meta["Nx"], meta["Ny"], meta["Nz"]
    I = np.arange(Nx)[:, None, None]
    J = np.arange(Ny + 1)[None, :, None]
    K = np.arange(Nz)[None, None, :]
    bnd = np.zeros((Nx, Ny + 1, Nz), dtype=bool)
    bnd[0] = bnd[-1] = True
    bnd[:, :, 0] = bnd[:, :, -1] = True
    res = {}
    for grp, comps in (("E", COMPS[:3]), ("H", COMPS[3:])):
        worst = worst_b = worst_i = 0.0
        scale = 0.0
        diffs = []
        for c in comps:
            nstep = NSTEPS if c[0] == "E" else NSTEPS - 1
            ref = analytic(st, c, I, J, K, nstep)
            valid = np.ones_like(bnd)
            if meta["offsets"][c][1] == 0.5:
                valid[:, -1, :] = False
            d = np.abs(F[c] - ref)
            d[~valid] = 0.0
            scale = max(scale, float(np.abs(ref[valid]).max()))
            diffs.append(d)
        for d in diffs:
            worst = max(worst, d.max() / scale)
            worst_b = max(worst_b, d[bnd].max() / scale)
            worst_i = max(worst_i, d[~bnd].max() / scale)
        res[grp] = dict(max_rel=float(worst), boundary=float(worst_b), interior=float(worst_i))
    passed = res["E"]["max_rel"] < TOL and res["H"]["max_rel"] < TOL and theory < 1e-12
    return dict(nl=nl, pol=pol, m=m, n=n, nsteps=NSTEPS, theory_mismatch=theory, **res, passed=bool(passed))


def main():
    fdtd_io.build()
    rows = [check(nl, pol, m, n) for nl in (10, 20) for pol in ("s", "p") for (m, n) in ((1, 1), (-1, 1))]
    print(f"{'nl':>3} {'pol':>3} {'(m,n)':>7} | {'E max':>9} {'E bnd':>9} {'E int':>9} | "
          f"{'H max':>9} {'H bnd':>9} {'H int':>9} | {'C vs py':>9} | PASS")
    for r in rows:
        print(f"{r['nl']:>3} {r['pol']:>3} {str((r['m'], r['n'])):>7} | {r['E']['max_rel']:9.2e} "
              f"{r['E']['boundary']:9.2e} {r['E']['interior']:9.2e} | {r['H']['max_rel']:9.2e} "
              f"{r['H']['boundary']:9.2e} {r['H']['interior']:9.2e} | {r['theory_mismatch']:9.2e} | "
              f"{'PASS' if r['passed'] else 'FAIL'}")
    ok = all(r["passed"] for r in rows)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "stage1_multistep.json"), "w") as fh:
        json.dump(dict(tolerance=TOL, rows=rows, passed=ok), fh, indent=2)
    print("STAGE 1 multistep:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
