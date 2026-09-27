"""Stage B gate 0 (SPEC_nonuniform §20.5 B0-1, B0-2; amendment D16).

B0-1: x (or z) nonuniform + inc=a, or + auxref=1 -> exit 2 with the exact message.
B0-2: gate-0 checks 0-2 (adjointness + C kernel), 0-3 (energy 1e5 steps), 0-4 (stability limit), 0-5 (divergence)
      on grids graded in x AND z (B_ops_xz_graded_pec) / in x (B_vac_k1 with the analytic injection inc=p).
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import analyze_nu as an  # noqa: E402

TH = nc.thresholds()
G0 = TH["gate0"]
T = nc.Table("B_gate0")
MSG = "phasor aux line requires uniform x and z"
GR = "B_ops_xz_graded_pec"


def main():
    nc.build()
    base = os.path.join(nc.RUNS, "stageB", "gate0")
    # ---- B0-1
    for kw, what in ((dict(inc="a"), "inc=a"), (dict(inc="p", auxref=1), "inc=p + auxref=1")):
        _, rc, err = nc.run(os.path.join(base, "abort"), grid="B_vac_k1", nsteps=5, energy_every=0, m=1, n=0,
                            expect_fail=True, **kw)
        ok = rc == 2 and MSG in err
        T.row("B0-1", f"x nonuniform + {what} is refused", "exit 2 + message", f"exit {rc}: {err.strip()[:90]}",
              "exact", ok, grid="B_vac_k1")
    # ---- B0-2 / 0-2 adjointness + kernel
    out = os.path.join(base, "geom")
    nc.run(out, grid=GR, inc="0", init="r", nsteps=2, dump_at=[1, 2], energy_every=0)
    used = nc.load_used(out)
    geo = nc.Geo(nc.load_grid(GR))
    W = an.weights(used)
    for seed in (1, 2):
        F = an.random_state(used, seed)
        CH, CE = an.curl_H(F, used), an.curl_E(F, used)
        lhs = an.inner(F, CH, W, ("Ex", "Ey", "Ez"))
        rhs = an.inner(CE, F, W, ("Hx", "Hy", "Hz"))
        rel = abs(lhs - rhs) / (np.sqrt(an.inner(F, F, W, ("Ex", "Ey", "Ez"))) * np.sqrt(an.inner(CH, CH, W, ("Ex", "Ey", "Ez"))))
        T.row("B0-2 (0-2)", f"adjointness, x and z graded, seed {seed}", "0", f"{rel:.2e}", f"< {G0['0-2']['rel']:g}",
              rel < G0["0-2"]["rel"], grid=GR)
    meta = nc.load_meta(out)
    F1, F2 = nc.load_fields(out, "fields_n1.bin", meta), nc.load_fields(out, "fields_n2.bin", meta)
    eps = an.eps_fields(used)
    CE1 = an.curl_E(F1, used)
    eH = max(np.abs(F2[c] - (F1[c] - meta["dt"] * CE1[c])).max() / np.abs(F1[c]).max() for c in ("Hx", "Hy", "Hz"))
    CH2 = an.curl_H(F2, used)
    eE = max(np.abs(F2[c] - (F1[c] + meta["dt"] * CH2[c] / eps[c])).max() / np.abs(F1[c]).max() for c in ("Ex", "Ey", "Ez"))
    T.row("B0-2 (0-2)", "C kernel vs Python operator, x and z graded", "0", f"H {eH:.1e}, E {eE:.1e}", "< 1e-13",
          max(eH, eE) < 1e-13, grid=GR)
    # ---- 0-3 energy
    out = os.path.join(base, "energy")
    m3, _, _ = nc.run(out, grid=GR, inc="0", init="r", seed=5, nsteps=G0["0-3"]["nsteps"], energy_every=100)
    L = nc.load_log(out)
    Wm = L["W_mod"][~np.isnan(L["W_mod"])]
    drift = float(np.max(np.abs(Wm / Wm[0] - 1)))
    T.row("B0-2 (0-3)", f"energy drift over {G0['0-3']['nsteps']:.0e} steps, x and z graded", "0", f"{drift:.2e}",
          f"< {G0['0-3']['rel_drift']:g}", drift < G0["0-3"]["rel_drift"], grid=nc.grid_info(geo, m3["dt"]))
    # ---- 0-4 stability limit
    out = os.path.join(base, "eig")
    nc.run(out, grid=GR, inc="0", mode="e", eig_tol=G0["0-4"]["power_iter_rq_tol"], seed=3)
    with open(os.path.join(out, "eig.json")) as fh:
        E = json.load(fh)
    dtmax = E["dt_max"]
    dtc = geo.dt(S=1.0 / np.sqrt(3.0))
    T.row("B0-2 (0-4)", "Courant formula Δt ≤ Δt_max (power iteration), x and z graded", f"{dtc:.10g}",
          f"{dtmax:.10g} ({E['iterations']} it)", "Δt_C ≤ Δt_max", dtc <= dtmax and E["last_rel_change"] < 1e-10, grid=GR)
    ns = G0["0-4"]["nsteps"]
    _, rc_s, _ = nc.run(os.path.join(base, "s099"), grid=GR, inc="0", init="r", seed=11, dt=0.99 * dtmax, nsteps=ns,
                        energy_every=100, expect_fail=True)
    _, rc_u, err = nc.run(os.path.join(base, "s102"), grid=GR, inc="0", init="r", seed=11, dt=1.02 * dtmax, nsteps=ns,
                          energy_every=100, expect_fail=True)
    T.row("B0-2 (0-4)", f"0.99 Δt_max stable / 1.02 Δt_max diverges ({ns} steps)", "stable / diverges",
          f"exit {rc_s} / exit {rc_u}", "0 / non-finite", rc_s == 0 and rc_u != 0 and "instability" in err, grid=GR)
    # ---- 0-5 divergence with the analytic injection (x graded), final-step dump
    out = os.path.join(base, "div")
    gv = nc.Geo(nc.load_grid("B_vac_k1"))
    dt = gv.dt_run()
    Np = int(round(1 / dt))
    md, _, _ = nc.run(out, grid="B_vac_k1", inc="p", m=1, n=0, dt=repr(dt), nsteps=40 * Np, dump_at=[40 * Np],
                      energy_every=0, div_every=Np)
    used = nc.load_used(out)
    F = nc.load_fields(out, f"fields_n{40 * Np}.bin", md)
    lo, hi = gv.j0 + 2, gv.Ny - gv.npml_hi - 2
    K = md["omega_tilde"]
    Emax = max(np.abs(F[c]).max() for c in ("Ex", "Ey", "Ez"))
    Hmax = max(np.abs(F[c]).max() for c in ("Hx", "Hy", "Hz"))
    dE = np.abs(an.div_E(F, used)[:, lo:hi + 1, :]).max() / (K * Emax)
    dH = np.abs(an.div_H(F, used)[:, lo:hi + 1, :]).max() / (K * Hmax)
    dEu = np.abs(an.div_E(F, used, uniform_delta=used["hx"][0])[:, lo:hi + 1, :]).max() / (K * Emax)
    T.row("B0-2 (0-5)", "divergence in the TF interior, x graded, inc=p", "0", f"E {dE:.1e}, H {dH:.1e}",
          f"< {G0['0-5']['rel']:g}", max(dE, dH) < G0["0-5"]["rel"], grid="B_vac_k1",
          note=f"control (uniform operator): {dEu:.1e}")
    T.save()
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
