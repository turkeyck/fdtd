"""Stage 2: CPML stability without a source.

(a) Energy diagnostic check: the same packet in a lossless PEC/PBC cavity (npml = 0) must keep the
    Yee-conserved energy W^{n+1/2} = 1/2 sum(eps E^n.E^{n+1} + |H^{n+1/2}|^2) constant to round-off.
(b) CPML run: a divergence-free packet (Ez(x,y), Ex(y,z); carrier ky0 = 5, Gaussian width 1.2) starts
    in the middle and must leave through the y-CPMLs. Over 5000 steps: finite everywhere (the solver
    aborts otherwise), W_total and W_phys never increase between one-period samples, and W decays
    below 1e-6 of its maximum.
Writes results/stage2_cpml.json and figures/stage2/cpml_energy.png.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import fdtd_io  # noqa: E402

NSTEPS = 5000


def growth_stats(W, per):
    Wp = W[per - 1::per]
    return float((np.diff(Wp) / Wp[:-1]).max()), float((np.diff(W) / W[:-1]).max())


def main():
    fdtd_io.build()
    cav = os.path.join(ROOT, "runs", "stage2", "cavity")
    mc = fdtd_io.run(cav, inc="0", init="b", npml=0, sf=20, tf=240, nsteps=1000, energy_every=1)
    Wc = fdtd_io.load_log(cav)["W_total"]
    cav_drift = float(np.abs(Wc / Wc[0] - 1).max())

    out = os.path.join(ROOT, "runs", "stage2", "blob")
    meta = fdtd_io.run(out, quiet=False, inc="0", init="b", npml=20, sf=20, tf=200, nsteps=NSTEPS,
                       energy_every=1)
    log = fdtd_io.load_log(out)
    t, W, Wph = log["t"], log["W_total"], log["W_phys"]
    per = int(round(1.0 / meta["dt"]))
    gT, gT_step = growth_stats(W, per)
    gP, gP_step = growth_stats(Wph, per)
    decay = float(W[-1] / W.max())
    passed = bool(np.all(np.isfinite(W)) and cav_drift < 1e-12 and gT <= 0.0 and gP <= 0.0 and decay < 1e-6)
    res = dict(nsteps=NSTEPS, grid=[meta["Nx"], meta["Ny"] + 1, meta["Nz"]], Delta=meta["Delta"], dt=meta["dt"],
               pml=meta["pml"], cavity_energy_drift=cav_drift, W_initial=float(W[0]), W_final=float(W[-1]),
               decay_ratio=decay, max_growth_per_period_total=gT, max_growth_per_step_total=gT_step,
               max_growth_per_period_phys=gP, max_growth_per_step_phys=gP_step,
               runtime_s=meta["runtime_s"], passed=passed)
    print(json.dumps(res, indent=2))
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "stage2_cpml.json"), "w") as fh:
        json.dump(res, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(ROOT, "figures", "stage2"), exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    ax.semilogy(t, W / W.max(), lw=1, label="W_total (incl. PML)")
    ax.semilogy(t, Wph / W.max(), lw=1, ls="--", label="W_phys (outside PML)")
    ax.axhline(1e-6, color="k", lw=0.8, ls=":", label="1e-6")
    ax.set_xlabel("t  [λ0/c]")
    ax.set_ylabel("W / max W")
    ax.set_title(f"CPML, no source, Δ=λ0/{meta['nl']}, N_pml={meta['npml']}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.plot(fdtd_io.load_log(cav)["t"], Wc / Wc[0] - 1, lw=1)
    ax.set_xlabel("t  [λ0/c]")
    ax.set_ylabel("W(t)/W(0) - 1")
    ax.set_title("Lossless PEC/PBC cavity: Yee energy conservation")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "stage2", "cpml_energy.png"), dpi=130)
    print("STAGE 2 CPML stability:", "PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
