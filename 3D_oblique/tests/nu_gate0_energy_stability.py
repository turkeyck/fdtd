"""Gate 0-3 (energy conservation, 1e5 steps, PEC y / PBC x-z, random fields, no source) and
gate 0-4 (stability limit from power iteration in the C solver; 0.99 dt_max stable, 1.02 dt_max diverges;
Courant formula dt <= dt_max). SPEC_nonuniform §20.1. Grids: the gate-1 taper and abrupt grids with PEC walls on
a thin x-z cell (4 x 6 cells of the same spacing; the y operator carries the nonuniformity)."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402

TH = nc.thresholds()["gate0"]
GRIDS = ["L1_taper_r1.1_base20", "L1_abrupt_r4_base20"]


def figure(curves, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for lab, (n, W) in curves["energy"].items():
        axes[0].plot(n, np.abs(W / W[0] - 1) + 1e-18, label=lab)
    axes[0].set_yscale("log")
    axes[0].axhline(TH["0-3"]["rel_drift"], color="#b91c1c", ls="--", lw=0.8, label="threshold 1e-10")
    axes[0].set_xlabel("step n")
    axes[0].set_ylabel("|W_mod(n)/W_mod(0) − 1|")
    axes[0].set_title("0-3: conserved discrete energy (PEC y, PBC x/z)")
    axes[0].legend(fontsize=7)
    for lab, (it, lam) in curves["eig"].items():
        axes[1].semilogy(it, np.abs(lam / lam[-1] - 1) + 1e-18, label=lab)
    axes[1].set_xlabel("iteration")
    axes[1].set_ylabel("|λ_k/λ_final − 1|")
    axes[1].set_title("0-4: power iteration (C solver)")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    nc.build()
    T = nc.Table("gate0_energy_stability")
    curves = {"energy": {}, "eig": {}}
    base = os.path.join(nc.RUNS, "gate0")
    for g in GRIDS:
        thin = g + "_pec_thin"
        geo = nc.Geo(nc.load_grid(thin))
        # ---- 0-3 energy
        out = os.path.join(base, f"energy_{thin}")
        n3 = TH["0-3"]["nsteps"]
        meta, _, _ = nc.run(out, grid=thin, inc="0", init="r", seed=7, nsteps=n3, energy_every=100)
        L = nc.load_log(out)
        W = L["W_mod"][~np.isnan(L["W_mod"])]
        n = L["n"][~np.isnan(L["W_mod"])]
        drift = float(np.max(np.abs(W / W[0] - 1)))
        curves["energy"][g] = (n, W)
        T.row("0-3", f"energy W_mod relative drift over {n3:.0e} steps ({g}, PEC)", "0", f"{drift:.2e}",
              f"< {TH['0-3']['rel_drift']:g}", drift < TH["0-3"]["rel_drift"], grid=nc.grid_info(geo, meta["dt"]))
        Wo = L["W_total"][~np.isnan(L["W_total"])]
        T.row("0-3", f"legacy-form energy ½Σ(εE^n·E^(n+1)+|H^(n+½)|²) drift ({g})", "0",
              f"{np.max(np.abs(Wo / Wo[0] - 1)):.2e}", "INFO", None, grid=thin)
        # ---- 0-4 stability limit
        out = os.path.join(base, f"eig_{thin}")
        nc.run(out, grid=thin, inc="0", mode="e", eig_tol=TH["0-4"]["power_iter_rq_tol"], seed=3)
        with open(os.path.join(out, "eig.json")) as fh:
            E = json.load(fh)
        lg = np.genfromtxt(os.path.join(out, "eig_log.csv"), delimiter=",", names=True)
        curves["eig"][g] = (np.atleast_1d(lg["it"]), np.atleast_1d(lg["rayleigh"]))
        dtmax = E["dt_max"]
        dtc = geo.dt(S=1.0 / np.sqrt(3.0))       # the minimum-spacing Courant bound (S_max)
        T.row("0-4", f"power iteration λ_max, Δt_max = 2/sqrt(λ_max) ({g})", "—",
              f"λ={E['lambda_max']:.10g}, Δt_max={dtmax:.10g} T0 ({E['iterations']} it)", "RQ change < 1e-10",
              E["last_rel_change"] < TH["0-4"]["power_iter_rq_tol"], grid=thin)
        T.row("0-4", f"Courant formula Δt ≤ Δt_max ({g})", f"{dtc:.10g}", f"Δt_max {dtmax:.10g}", "Δt_C ≤ Δt_max",
              dtc <= dtmax, grid=thin, note=f"margin {dtmax / dtc - 1:.2%}")
        dt_base = 0.5 * geo.hy.max()
        T.row("0-4", f"Δt cost of the smallest cell ({g})", "—",
              f"Δt(S=0.5) = {geo.dt():.6g} T0 vs {dt_base:.6g} on the base grid: ×{dt_base / geo.dt():.3f} steps",
              "INFO", None, grid=nc.grid_info(geo))
        ns = TH["0-4"]["nsteps"]
        out_s = os.path.join(base, f"stab099_{thin}")
        ms, rc, err = nc.run(out_s, grid=thin, inc="0", init="r", seed=11, dt=0.99 * dtmax, nsteps=ns,
                             energy_every=100, expect_fail=True)
        ok_s = rc == 0
        dr = np.nan
        if ok_s:
            L = nc.load_log(out_s)
            Wm = L["W_mod"][~np.isnan(L["W_mod"])]
            dr = float(np.max(np.abs(Wm / Wm[0] - 1)))
            ok_s = dr < TH["0-3"]["rel_drift"]
        T.row("0-4", f"0.99 Δt_max: {ns} steps finite, energy drift ({g})", "stable",
              f"exit {rc}, drift {dr:.2e}", "finite & drift < 1e-10", ok_s, grid=thin)
        out_u = os.path.join(base, f"stab102_{thin}")
        mu, rc, err = nc.run(out_u, grid=thin, inc="0", init="r", seed=11, dt=1.02 * dtmax, nsteps=ns,
                             energy_every=100, expect_fail=True)
        blew = rc != 0 and "instability detected" in err
        step = err.split("non-finite field at step")[-1].split()[0] if blew else "—"
        T.row("0-4", f"1.02 Δt_max: diverges within {ns} steps ({g})", "diverges", f"exit {rc}, step {step}",
              "non-finite before the end", blew, grid=thin)
    figure(curves, os.path.join(nc.FIGS, "gate0", "energy_eig.png"))
    T.save()
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
