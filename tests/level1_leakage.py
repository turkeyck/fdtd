"""Level 1-1 / Stage 3: TF/SF leakage into the scattered-field region.

Main test (derivation.md §6): incident field from the 1D modal line, erf ramp (default) and the
10-period raised cosine; max |E| in the SF region (PML excluded) over the causal window that ends
before any far-PML echo can reach y0 (signal speed bounded by c), relative to |E0| = 1.

Controls (analytic incident field x ramp, the TF/SF plane fed directly by derivation.md §2.4):
  * 'analytic, discrete ky'   -- shows the ramp itself leaks (derivation.md §6.1)
  * 'analytic, continuous ky' -- steady-state leakage from the ky mismatch. Measured as the DFT at w0
    on an SF plane over a post-ramp window, minus the same quantity of the discrete-ky run (removes the
    shared far-PML echo and ramp remnant by linearity). Must scale as Delta^2 (nl = 20 vs 40).
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import fdtd_io  # noqa: E402

TOL = 1e-10
RUNS = os.path.join(ROOT, "runs", "level1_leak")


def causal_end_step(meta, src_j, margin_periods=1.0):
    """Last step before a signal leaving src_j at t=0 at speed c can reflect off the far PML face and reach y0."""
    D = meta["Delta"]
    y_far = (meta["Ny"] - meta["npml"]) * D
    y0 = meta["j0"] * D
    t_end = (y_far - src_j * D) + (y_far - y0) - margin_periods
    return int(math.floor(t_end / meta["dt"]))


def geom(nl):
    return dict(nl=nl, npml=20, sf=nl, tf=10 * nl)  # SF 1 lambda0, TF 10 lambda0


def main_runs():
    rows = []
    for pol in ("s", "p"):
        for ramp, extra in (("e", {}), ("r", {"ramp_T": 10})):
            nl = 20
            out = os.path.join(RUNS, f"aux_{pol}_{ramp}")
            g = geom(nl)
            # run past the causal window so the echo arrival is visible in the figure
            meta = fdtd_io.run(out, inc="a", pol=pol, ramp=ramp, nsteps=1400, energy_every=0, **g, **extra)
            log = fdtd_io.load_log(out)
            n_end = causal_end_step(meta, meta["aux"]["ja"])
            sel = log["n"] <= n_end
            leak_E = float(log["maxE_SF"][sel].max())
            leak_H = float(log["maxH_SF"][sel].max())  # |H0| = |E0|/eta0 = 1
            rows.append(dict(pol=pol, ramp={"e": "erf(t0=20,tau=4)", "r": "raised-cosine 10T0"}[ramp], nl=nl,
                             Delta=meta["Delta"], dt=meta["dt"], causal_end_step=n_end,
                             leak_E=leak_E, leak_H=leak_H, passed=bool(max(leak_E, leak_H) < TOL), out=out))
    return rows


def control_runs():
    recs = []
    for nl in (20, 40):
        g = geom(nl)
        per = int(round(nl / 0.5))           # steps per period
        d0, d1 = 48 * per, 64 * per          # after the erf ramp (t0 + 5 tau = 40 T0) and its transient
        jplane = g["npml"] + g["sf"] // 2    # SF plane half-way between near PML and y0
        res = {}
        for ky in ("d", "c"):
            out = os.path.join(RUNS, f"analytic_{ky}_nl{nl}")
            meta = fdtd_io.run(out, inc="n", ky=ky, nsteps=d1, dft0=d0, dft1=d1, yplanes=[jplane],
                               energy_every=0, **g)
            F, _ = fdtd_io.load_dft(out, f"y{jplane}", meta)
            log = fdtd_io.load_log(out)
            n_end = causal_end_step(meta, meta["j0"])
            res[ky] = dict(F=F, meta=meta, window_leak=float(log["maxE_SF"][log["n"] <= n_end].max()),
                           log=log, out=out)
        diff = max(float(np.abs(res["c"]["F"][c] - res["d"]["F"][c]).max()) for c in ("Ex", "Ey", "Ez"))
        floor = max(float(np.abs(res["d"]["F"][c]).max()) for c in ("Ex", "Ey", "Ez"))
        m = res["c"]["meta"]
        recs.append(dict(nl=nl, Delta=m["Delta"], k0Delta=2 * math.pi * m["Delta"], dft_window=[d0, d1],
                         sf_plane=jplane, ky_leak_steady=diff, discrete_run_sf_amplitude=floor,
                         window_leak_discrete=res["d"]["window_leak"], window_leak_continuous=res["c"]["window_leak"],
                         ky_cont=m["ky_cont"], ky_disc=m["ky_disc"]))
        recs[-1]["_logs"] = (res["d"]["log"], res["c"]["log"])
    slope = math.log(recs[0]["ky_leak_steady"] / recs[1]["ky_leak_steady"]) / math.log(2.0)
    return recs, slope


def main():
    fdtd_io.build()
    rows = main_runs()
    print(f"{'pol':>3} {'ramp':>20} | {'leak E':>9} {'leak H':>9} | {'window end':>10} | PASS")
    for r in rows:
        print(f"{r['pol']:>3} {r['ramp']:>20} | {r['leak_E']:9.2e} {r['leak_H']:9.2e} | "
              f"{r['causal_end_step']:>10} | {'PASS' if r['passed'] else 'FAIL'}")
    recs, slope = control_runs()
    print("\nControls (s, (1,1)):")
    for c in recs:
        print(f"nl={c['nl']:>3}  steady ky-leak |E|={c['ky_leak_steady']:.3e}  "
              f"(discrete-ky run SF amplitude {c['discrete_run_sf_amplitude']:.2e}; "
              f"causal-window max: disc {c['window_leak_discrete']:.2e}, cont {c['window_leak_continuous']:.2e})")
    slope_ok = abs(slope - 2.0) < 0.1
    print(f"ky-leak scaling exponent (nl 20 -> 40): {slope:.3f}  -> {'PASS' if slope_ok else 'FAIL'} (2 ± 0.1)")
    ok = all(r["passed"] for r in rows) and slope_ok

    make_figure(rows, recs)
    out = dict(tolerance=TOL, main=[{k: v for k, v in r.items() if k != "out"} for r in rows],
               controls=[{k: v for k, v in c.items() if not k.startswith("_")} for c in recs],
               control_scaling_exponent=slope, control_scaling_ok=slope_ok, passed=bool(ok))
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "level1_leakage.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("LEVEL 1-1 leakage:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def make_figure(rows, recs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(ROOT, "figures", "level1"), exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
    ax = axes[0]
    for r in rows:
        if r["pol"] != "s":
            continue
        log = fdtd_io.load_log(r["out"])
        ax.semilogy(log["t"], np.maximum(log["maxE_SF"], 1e-18), lw=1, label=f"modal line, {r['ramp']}")
    ld, lc = recs[0]["_logs"]
    ax.semilogy(ld["t"], ld["maxE_SF"], lw=1, label="analytic×ramp, discrete ky")
    ax.semilogy(lc["t"], lc["maxE_SF"], lw=1, label="analytic×ramp, continuous ky")
    tend = rows[0]["causal_end_step"] * rows[0]["dt"]
    ax.axvline(tend, color="k", ls="--", lw=0.8, label="causal window end (modal line)")
    ax.axhline(TOL, color="r", ls=":", lw=0.8, label="threshold 1e-10")
    ax.set_xlim(0, 64)
    ax.set_xlabel("t  [λ0/c]")
    ax.set_ylabel("max |E| in SF region / |E0|")
    ax.set_title("SF leakage, s-pol, Δ=λ0/20, S=0.5")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    kd = np.array([c["k0Delta"] for c in recs])
    lk = np.array([c["ky_leak_steady"] for c in recs])
    ax.loglog(kd, lk, "o-", label="steady ky-mismatch leak (cont − disc)")
    ax.loglog(kd, lk[0] * (kd / kd[0]) ** 2, "k:", label="slope 2")
    ax.set_xlabel("k0 Δ")
    ax.set_ylabel("|E| on SF plane (DFT at ω0) / |E0|")
    ax.set_title("Continuous-ky control scales as Δ²")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "level1", "L1_1_leakage.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
