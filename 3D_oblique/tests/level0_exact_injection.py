"""Level 0: prove the injected discrete plane wave is an exact solution of the PBC Yee update.

Plugs the analytic fields (derivation.md §2) into one Yee step with x/z periodic
index wrap (derivation.md §4) and checks residuals, discrete divergences, and the
continuous-ky controls. Writes results/level0.json and figures/level0/*.png.
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fdtd_theory import (ECOMP, HCOMP, EPS0, MU0, Setup, analytic, ktilde)  # noqa: E402

TOL_RES = 1e-12
TOL_DIV = 1e-12
NY = 8  # y extent of the test block; y is not periodic, only interior j in [1, NY-1] is compared


def fields(setup, nE, nH, continuous):
    I = np.arange(setup.Nx)[:, None, None]
    J = np.arange(NY + 1)[None, :, None]
    K = np.arange(setup.Nz)[None, None, :]
    f = {c: analytic(setup, c, I, J, K, nE, continuous=continuous) for c in ECOMP}
    f.update({c: analytic(setup, c, I, J, K, nH, continuous=continuous) for c in HCOMP})
    return f


# Periodic differences. "fwd" is used where the target sits half a cell *after* the
# integer-indexed source (F[i+1]-F[i], wrap F[N]=F[0]); "bwd" where the target sits
# half a cell *before* (F[i]-F[i-1], wrap F[-1]=F[N-1]).  See derivation.md §4.
def fwd(F, ax):
    return np.roll(F, -1, axis=ax) - F


def bwd(F, ax):
    return F - np.roll(F, 1, axis=ax)


def curl_E(f, D):
    """(curl E) at the H points. y differences are only valid for j <= NY-1."""
    return {
        "Hx": (fwd(f["Ez"], 1) - fwd(f["Ey"], 2)) / D,
        "Hy": (fwd(f["Ex"], 2) - fwd(f["Ez"], 0)) / D,
        "Hz": (fwd(f["Ey"], 0) - fwd(f["Ex"], 1)) / D,
    }


def curl_H(f, D):
    """(curl H) at the E points. y differences are only valid for j >= 1."""
    return {
        "Ex": (bwd(f["Hz"], 1) - bwd(f["Hy"], 2)) / D,
        "Ey": (bwd(f["Hx"], 2) - bwd(f["Hz"], 0)) / D,
        "Ez": (bwd(f["Hy"], 0) - bwd(f["Hx"], 1)) / D,
    }


def div_E(f, D):
    return (bwd(f["Ex"], 0) + bwd(f["Ey"], 1) + bwd(f["Ez"], 2)) / D  # at (i, j, k), valid j >= 1


def div_H(f, D):
    return (fwd(f["Hx"], 0) + fwd(f["Hy"], 1) + fwd(f["Hz"], 2)) / D  # at cell centres, valid j <= NY-1


INTERIOR = (slice(None), slice(1, NY), slice(None))


def boundary_mask(setup):
    m = np.zeros((setup.Nx, NY + 1, setup.Nz), dtype=bool)
    m[0, :, :] = m[-1, :, :] = True
    m[:, :, 0] = m[:, :, -1] = True
    return m[INTERIOR]


def one_step(setup, continuous=None, nstep=7):
    """Residuals of one H half-step and one E half-step, relative to max|field| (spec definition)."""
    D, dt = setup.Delta, setup.dt
    now = fields(setup, nstep, nstep - 1, continuous)      # E^n, H^{n-1/2}
    nxt = fields(setup, nstep + 1, nstep, continuous)      # E^{n+1}, H^{n+1/2}
    cE = curl_E(now, D)
    cH = curl_H(nxt, D)  # uses the analytic H^{n+1/2}, so E and H residuals are independent
    bmask = boundary_mask(setup)
    out = {}
    for group, comps in (("E", ECOMP), ("H", HCOMP)):
        scale = max(np.max(np.abs(nxt[c][INTERIOR])) for c in comps)
        worst, worst_b, worst_i = 0.0, 0.0, 0.0
        per = {}
        for c in comps:
            if group == "H":
                upd = now[c] - dt / MU0 * cE[c]
            else:
                upd = now[c] + dt / EPS0 * cH[c]
            r = np.abs(upd - nxt[c])[INTERIOR]
            per[c] = float(r.max() / scale)
            worst = max(worst, per[c])
            worst_b = max(worst_b, float(r[bmask].max() / scale))
            worst_i = max(worst_i, float(r[~bmask].max() / scale))
        out[group] = {"max_rel": worst, "per_component": per,
                      "boundary_cells": worst_b, "interior_cells": worst_i}
    # Divergences, normalized by |K~| * max|field| (each difference term is O(|K~| |F|)).
    Kn = float(np.linalg.norm(setup.Kt))
    Escale = max(np.max(np.abs(now[c][INTERIOR])) for c in ECOMP)
    Hscale = max(np.max(np.abs(now[c][INTERIOR])) for c in HCOMP)
    out["divE_rel"] = float(np.abs(div_E(now, D)[INTERIOR]).max() / (Kn * Escale))
    out["divH_rel"] = float(np.abs(div_H(now, D)[INTERIOR]).max() / (Kn * Hscale))
    return out


def predicted_control(setup, continuous):
    """Closed-form per-step residual amplitude (relative) for the continuous-ky controls."""
    w, dt = setup.omega, setup.dt
    kc = np.array([setup.kx, setup.ky_cont, setup.kz])
    Kc = ktilde(kc, setup.Delta)
    e, h = setup.amplitudes(1.0, continuous)
    ph = np.exp(-1j * w * dt) - 1.0
    half = np.exp(-1j * w * dt / 2.0)
    rE = ph * e - 1j * dt * half * np.cross(Kc, h) / EPS0
    rH = ph * h + 1j * dt * half * np.cross(Kc, e) / MU0
    return float(np.abs(rE).max() / np.abs(e).max()), float(np.abs(rH).max() / np.abs(h).max())


def main():
    cases = []
    for nl in (10, 20, 40):
        for pol in ("s", "p"):
            for (m, n) in ((1, 1), (-1, 1), (1, -1)):
                cases.append(dict(n_lambda=nl, pol=pol, m=m, n=n))

    results = {"tolerance_residual": TOL_RES, "tolerance_div": TOL_DIV, "cases": [], "controls": []}
    ok = True
    print(f"{'N/lam':>5} {'pol':>3} {'(m,n)':>7} | {'E res':>9} {'H res':>9} {'E res bnd':>9} "
          f"{'H res bnd':>9} | {'divE':>9} {'divH':>9} | PASS")
    for cs in cases:
        st = Setup(**cs)
        r = one_step(st)
        passed = (r["E"]["max_rel"] < TOL_RES and r["H"]["max_rel"] < TOL_RES
                  and r["divE_rel"] < TOL_DIV and r["divH_rel"] < TOL_DIV)
        ok &= passed
        results["cases"].append(dict(cs, Delta=st.Delta, dt=st.dt, ky=st.ky, **r, passed=passed))
        print(f"{cs['n_lambda']:>5} {cs['pol']:>3} {str((cs['m'], cs['n'])):>7} | "
              f"{r['E']['max_rel']:9.2e} {r['H']['max_rel']:9.2e} {r['E']['boundary_cells']:9.2e} "
              f"{r['H']['boundary_cells']:9.2e} | {r['divE_rel']:9.2e} {r['divH_rel']:9.2e} | "
              f"{'PASS' if passed else 'FAIL'}")

    print("\nControls (continuous ky), s and p, (m,n)=(1,1):")
    print(f"{'ctrl':>5} {'pol':>3} {'N/lam':>5} | {'E res':>9} {'pred':>9} | {'H res':>9} {'pred':>9} | "
          f"{'E res/(w dt)':>12} | {'divE':>9}")
    for ctrl in ("ky", "full"):
        for pol in ("s", "p"):
            for nl in (10, 20, 40):
                st = Setup(n_lambda=nl, pol=pol)
                r = one_step(st, continuous=ctrl)
                pE, pH = predicted_control(st, ctrl)
                wdt = st.omega * st.dt
                rec = dict(control=ctrl, pol=pol, n_lambda=nl, Delta=st.Delta, k0Delta=st.k0 * st.Delta,
                           E_res=r["E"]["max_rel"], H_res=r["H"]["max_rel"], E_pred=pE, H_pred=pH,
                           E_res_per_rad=r["E"]["max_rel"] / wdt, H_res_per_rad=r["H"]["max_rel"] / wdt,
                           divE_rel=r["divE_rel"], divH_rel=r["divH_rel"],
                           ky_cont=st.ky_cont, ky_disc=st.ky)
                results["controls"].append(rec)
                print(f"{ctrl:>5} {pol:>3} {nl:>5} | {rec['E_res']:9.2e} {pE:9.2e} | {rec['H_res']:9.2e} "
                      f"{pH:9.2e} | {rec['E_res_per_rad']:12.3e} | {rec['divE_rel']:9.2e}")

    # Sensitivity: control residuals must dwarf the exact-case residuals, and scale as
    # (k0 Delta)^3 per step, i.e. (k0 Delta)^2 per radian of phase advance.
    slopes = {}
    for ctrl in ("ky", "full"):
        for pol in ("s", "p"):
            recs = [c for c in results["controls"] if c["control"] == ctrl and c["pol"] == pol]
            x = np.log([c["k0Delta"] for c in recs])
            y = np.log([c["E_res_per_rad"] for c in recs])
            slopes[f"{ctrl}_{pol}"] = float(np.polyfit(x, y, 1)[0])
    results["control_slopes_per_rad"] = slopes
    worst_exact = max(max(c["E"]["max_rel"], c["H"]["max_rel"]) for c in results["cases"])
    min_ctrl = min(c["E_res"] for c in results["controls"])
    sensitive = min_ctrl > 1e6 * worst_exact
    pred_ok = all(abs(c["E_res"] / c["E_pred"] - 1) < 0.02 for c in results["controls"])
    slope_ok = all(abs(s - 2.0) < 0.1 for s in slopes.values())
    results["sensitivity"] = dict(worst_exact=worst_exact, min_control=min_ctrl, ratio=min_ctrl / worst_exact,
                                  control_matches_prediction=pred_ok, control_slopes_ok=slope_ok,
                                  passed=bool(sensitive and pred_ok and slope_ok))
    print(f"\nlog-log slope of E residual per radian vs k0*Delta: {slopes}")
    print(f"sensitivity ratio min(control)/max(exact) = {min_ctrl / worst_exact:.2e}; "
          f"controls match closed form: {pred_ok}; slopes 2±0.1: {slope_ok}")
    ok &= results["sensitivity"]["passed"]
    results["passed"] = bool(ok)

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "level0.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    make_figures(results)
    print("\nLEVEL 0:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def make_figures(results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = os.path.join(ROOT, "figures", "level0")
    os.makedirs(out, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    labels, vals_e, vals_h = [], [], []
    for c in results["cases"]:
        if c["n_lambda"] == 20:
            labels.append(f"{c['pol']} ({c['m']},{c['n']})")
            vals_e.append(max(c["E"]["max_rel"], 1e-18))
            vals_h.append(max(c["H"]["max_rel"], 1e-18))
    for c in results["controls"]:
        if c["n_lambda"] == 20:
            labels.append(f"ctrl-{c['control']} {c['pol']}")
            vals_e.append(c["E_res"])
            vals_h.append(max(c["H_res"], 1e-18))
    xx = np.arange(len(labels))
    ax.bar(xx - 0.2, vals_e, 0.4, label="E residual")
    ax.bar(xx + 0.2, vals_h, 0.4, label="H residual")
    ax.axhline(TOL_RES, color="k", ls="--", lw=1, label="threshold 1e-12")
    ax.set_yscale("log")
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("max|update - analytic| / max|field|")
    ax.set_title("One Yee step, Δ = λ0/20, S = 0.5")
    ax.legend(fontsize=8)

    ax = axes[1]
    for ctrl, mk in (("ky", "o"), ("full", "s")):
        for pol, ls in (("s", "-"), ("p", "--")):
            recs = [c for c in results["controls"] if c["control"] == ctrl and c["pol"] == pol]
            ax.loglog([c["k0Delta"] for c in recs], [c["E_res_per_rad"] for c in recs], marker=mk, ls=ls,
                      label=f"control {ctrl}, {pol}")
    kd = np.array([c["k0Delta"] for c in results["controls"] if c["control"] == "ky" and c["pol"] == "s"])
    ref = [c["E_res_per_rad"] for c in results["controls"] if c["control"] == "ky" and c["pol"] == "s"]
    ax.loglog(kd, ref[1] * (kd / kd[1]) ** 2, "k:", label="slope 2")
    ax.set_xlabel("k0 Δ")
    ax.set_ylabel("E residual / (ω0 Δt)")
    ax.set_title("Continuous-ky controls: O((k0Δ)²) per radian")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "level0_residuals.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
