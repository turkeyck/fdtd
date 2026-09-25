"""Level 2: half-space of n = 1.5 inside the TF region -- Fresnel R, T, energy balance, Snell.

Geometry (cells scale with nl): near PML 20 | SF 1λ0 | TF vacuum 3λ0 | medium 5λ0 | far PML 20 (inside the medium).
Interface at y1 = j1*Δ on an integer (E-tangential) plane; Ex, Ez there use the arithmetic mean
(1+n²)/2, Ey (half-integer y) never lies on it (derivation.md §8). Fluxes use the discrete-conserved
Φ (derivation.md §9); Φ_inc is the exact discrete-plane-wave value (verified in L1-5):
    R = −Φ_SF / Φ_inc,   T = Φ_medium / Φ_inc.
Runs: (m,n)=(1,1) at Δ = λ0/10, λ0/20, λ0/40 for s and p (+ a staircase-interface variant at λ0/10, λ0/20);
optional Brewster sweep with n=0, Lz=2Δ, m=1 and Lx in [1.1, 2.0].
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import analyze as an  # noqa: E402
import fdtd_io  # noqa: E402
from fdtd_theory import discrete_fresnel, fresnel  # noqa: E402

RUNS = os.path.join(ROOT, "runs", "level2")
FIG = os.path.join(ROOT, "figures", "level2")
NMED = 1.5
T = []


def row(item, quantity, theory, measured, threshold, passed, note=""):
    T.append(dict(item=item, quantity=quantity, theory=theory, measured=measured, threshold=threshold,
                  passed=None if passed is None else bool(passed), note=note))
    flag = "INFO" if passed is None else ("PASS" if passed else "FAIL")
    print(f"[{item}] {quantity}: theory={theory} measured={measured} thr={threshold} -> {flag} {note}")


def layout(nl):
    npml, sf = 20, nl
    j0 = npml + sf
    j1 = j0 + 3 * nl
    return dict(npml=npml, sf=sf, tf=8 * nl, y1=3 * nl), j0, j1


def planes(nl):
    g, j0, j1 = layout(nl)
    jsf = g["npml"] + nl // 2
    vac = [j0 + nl, j0 + nl + 1, j0 + 2 * nl]
    med = [j1 + nl, j1 + nl + 1, j1 + 2 * nl, j1 + 3 * nl, j1 + 3 * nl + 1]
    return jsf, vac, med


def run_medium(nl, pol, ifmode="a", m=1, n=1, Lx=2.0, Lz=3.0, tag=""):
    per = 2 * nl
    g, j0, j1 = layout(nl)
    jsf, vac, med = planes(nl)
    out = os.path.join(RUNS, f"med_nl{nl}_{pol}_{ifmode}{tag}")
    Nz = int(round(Lz * nl))
    meta = fdtd_io.run(out, quiet=False, nl=nl, pol=pol, m=m, n=n, Lx=Lx, Lz=Lz, eps2=NMED ** 2, ifmode=ifmode,
                       inc="a", nsteps=110 * per, dft0=80 * per, dft1=110 * per, yplanes=[jsf] + vac + med,
                       zk=max(Nz // 3, 0), xi=nl // 2, energy_every=0, **g)
    return out, meta


def measure(out, meta, fit=True, plane_set=None):
    nl, D = meta["nl"], meta["Delta"]
    st = an.setup_from_meta(meta)
    stm = an.setup_from_meta(meta, n_ref=NMED)
    jsf, vac, med = plane_set or planes(nl)
    phi = lambda j: an.conserved_flux(fdtd_io.load_dft(out, f"y{j}", meta)[0], D)  # noqa: E731
    phi_inc = an.theory_flux(st, meta)
    phi_sf = phi(jsf)
    phi_vac = [phi(j) for j in vac]
    phi_med = [phi(j) for j in med]
    R = -phi_sf / phi_inc
    Tm = float(np.mean(phi_med)) / phi_inc
    cos1 = st.ky_cont / st.k0
    cos2 = stm.ky_cont / (NMED * stm.k0)
    F = fresnel(1.0, NMED, cos1, cos2)
    RF = F["Rs"] if meta["pol"] == "s" else F["Rp"]
    TF = 1.0 - RF
    eps_if = 0.5 * (1 + NMED ** 2) if meta["ifmode"] == "a" else NMED ** 2
    Dd = discrete_fresnel(st, stm, eps_if)
    res = dict(nl=nl, Delta=D, pol=meta["pol"], ifmode=meta["ifmode"], R=R, T=Tm, RT=R + Tm, R_fresnel=RF,
               T_fresnel=TF, R_relerr=(R - RF) / RF, T_relerr=(Tm - TF) / TF,
               R_disc=Dd["R"], T_disc=Dd["T"], R_vs_disc=(R - Dd["R"]) / Dd["R"], T_vs_disc=(Tm - Dd["T"]) / Dd["T"],
               R_disc_vs_fresnel=(Dd["R"] - RF) / RF,
               flux_var_medium=float((max(phi_med) - min(phi_med)) / abs(np.mean(phi_med))),
               flux_var_all=float((max(phi_vac + phi_med) - min(phi_vac + phi_med)) / abs(np.mean(phi_med))),
               theta1=math.degrees(math.acos(cos1)), theta2=math.degrees(math.acos(cos2)),
               ky_med_disc=stm.ky, ky_med_cont=stm.ky_cont)
    if fit:
        j1 = meta["j1"]
        sets, _ = an.phase_sets(out, meta, stm, ["xy", "yz"] + [f"y{j}" for j in med], j1 + 2,
                                meta["Ny"] - meta["npml"] - 1)
        k, rms, mx, _ = an.wavefront_fit(sets)
        res.update(k_med=k.tolist(), ky_med_relerr=abs(k[1] - stm.ky) / stm.ky, kx_med_relerr=abs(k[0] - st.kx) / st.kx,
                   kz_med_relerr=abs(k[2] - st.kz) / st.kz if st.kz else 0.0, fit_rms=rms)
    return res


def pml_study():
    """FDTD vs exact discrete Fresnel as the PML thickens (removes PML-echo contamination, SPEC-endorsed idea
    of lengthening the domain). Same geometry otherwise; the near and far PML both get N cells."""
    rows = []
    for nl in (10, 20):
        per = 2 * nl
        for pol in ("s", "p"):
            for N in (20, 40, 60):
                j0 = N + nl
                j1 = j0 + 3 * nl
                ps = (N + nl // 2, [j0 + nl, j0 + nl + 1, j0 + 2 * nl],
                      [j1 + nl, j1 + nl + 1, j1 + 2 * nl, j1 + 3 * nl, j1 + 3 * nl + 1])
                out = os.path.join(RUNS, f"pmlexp_nl{nl}_{pol}_N{N}")
                meta = fdtd_io.run(out, nl=nl, pol=pol, eps2=NMED ** 2, ifmode="a", inc="a", nsteps=110 * per,
                                   dft0=80 * per, dft1=110 * per, yplanes=[ps[0]] + ps[1] + ps[2], energy_every=0,
                                   npml=N, sf=nl, tf=8 * nl, y1=3 * nl)
                r = measure(out, meta, fit=False, plane_set=ps)
                rows.append(dict(nl=nl, pol=pol, N=N, R_vs_disc=r["R_vs_disc"], T_vs_disc=r["T_vs_disc"],
                                 RT=r["RT"]))
    return rows


def brewster_sweep():
    rows = []
    nl = 20
    for Lx in (1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.6, 2.0):
        for pol in ("p", "s"):
            out, meta = run_medium(nl, pol, m=1, n=0, Lx=Lx, Lz=2 * 1.0 / nl, tag=f"_brew_Lx{Lx}")
            r = measure(out, meta, fit=False)
            r["Lx"] = Lx
            rows.append(r)
    return rows


def main():
    fdtd_io.build()
    os.makedirs(FIG, exist_ok=True)
    res = {}
    for nl in (10, 20, 40):
        for pol in ("s", "p"):
            res[(nl, pol, "a")] = measure(*run_medium(nl, pol))
    for nl in (10, 20):
        for pol in ("s", "p"):
            res[(nl, pol, "s")] = measure(*run_medium(nl, pol, ifmode="s"), fit=False)

    for pol in ("s", "p"):
        r = res[(20, pol, "a")]
        tag = f"{pol}-pol, Δ=λ0/20, θ1={r['theta1']:.2f}°"
        row("L2", f"R vs continuous Fresnel ({tag})", f"{r['R_fresnel']:.6f}", f"{r['R']:.6f}", "rel err < 1%",
            abs(r["R_relerr"]) < 0.01, f"rel err {r['R_relerr']:+.3e}; absolute |ΔR| = {abs(r['R'] - r['R_fresnel']):.2e} "
            "(if '1%' meant 1 percentage point: INFO only, not used for the verdict)")
        row("L2", f"T vs continuous Fresnel ({tag})", f"{r['T_fresnel']:.6f}", f"{r['T']:.6f}", "rel err < 1%",
            abs(r["T_relerr"]) < 0.01, f"rel err {r['T_relerr']:+.3e}")
        row("L2", f"R, T vs exact discrete (Yee-lattice) Fresnel ({tag}, N_pml=20)", f"R {r['R_disc']:.10f}, T {r['T_disc']:.10f}",
            f"R {r['R']:.10f}, T {r['T']:.10f}", "see PML study", None,
            f"rel err R {r['R_vs_disc']:+.1e}, T {r['T_vs_disc']:+.1e} (PML-echo contaminated, see thickness study); "
            f"the discrete theory itself differs from continuous Fresnel by {r['R_disc_vs_fresnel']:+.2e} in R (derivation §8.1)")
        row("L2", f"R + T − 1 ({tag})", "0", f"{r['RT'] - 1:+.2e}", "|·| < 1e-4", abs(r["RT"] - 1) < 1e-4)
        row("L2", f"Snell: ky in medium vs discrete dispersion (n=1.5) ({tag})", f"{r['ky_med_disc']:.10f} /λ0",
            f"{r['k_med'][1]:.10f} /λ0", "rel err < 1e-6", r["ky_med_relerr"] < 1e-6,
            f"rel err {r['ky_med_relerr']:.1e}; continuous ky_m = {r['ky_med_cont']:.6f}; kx rel err "
            f"{r['kx_med_relerr']:.1e}; θ2 = {r['theta2']:.2f}°")
    for pol in ("s", "p"):
        for q in ("R", "T"):
            errs = np.array([abs(res[(nl, pol, "a")][f"{q}_relerr"]) for nl in (10, 20, 40)])
            sl = float(np.polyfit(np.log([1 / 10, 1 / 20, 1 / 40]), np.log(errs), 1)[0])
            row("L2", f"{q} convergence order ({pol}-pol, Δ=λ0/10,20,40)", "2 (O(Δ²))", f"{sl:.3f}",
                "2 ± 0.3", abs(sl - 2) < 0.3, "errors: " + ", ".join(f"{e:.2e}" for e in errs))
        for nl in (10, 20, 40):
            rr = res[(nl, pol, "a")]
            row("L2", f"R+T−1 ({pol}, Δ=λ0/{nl})", "0", f"{rr['RT'] - 1:+.2e}", "|·| < 1e-4", abs(rr["RT"] - 1) < 1e-4)
            if nl != 20:
                row("L2", f"R, T vs exact discrete Fresnel ({pol}, Δ=λ0/{nl}, N_pml=20)", f"R {rr['R_disc']:.8f}",
                    f"R {rr['R']:.8f}", "see PML study", None,
                    f"rel err R {rr['R_vs_disc']:+.1e}, T {rr['T_vs_disc']:+.1e}; vs continuous Fresnel R "
                    f"{rr['R_relerr']:+.2e}, T {rr['T_relerr']:+.2e}")
        for nl in (10, 20):
            rs = res[(nl, pol, "s")]
            row("L2", f"staircase interface (tangential ε = n² on the plane), {pol}, Δ=λ0/{nl}", "—",
                f"R err {rs['R_relerr']:+.2e}, T err {rs['T_relerr']:+.2e}", "—", None,
                "informational: shows why the arithmetic mean is used")
    study = pml_study()
    for nl in (10, 20):
        for pol in ("s", "p"):
            s = [q for q in study if q["nl"] == nl and q["pol"] == pol]
            eR = [abs(q["R_vs_disc"]) for q in s]
            eT = [abs(q["T_vs_disc"]) for q in s]
            mono = all(eR[i + 1] < eR[i] for i in range(2)) and all(eT[i + 1] < eT[i] for i in range(2))
            row("L2", f"FDTD → exact discrete Fresnel as N_pml = 20→40→60 ({pol}, Δ=λ0/{nl})", "decreasing to 0",
                "R: " + " → ".join(f"{e:.1e}" for e in eR) + "; T: " + " → ".join(f"{e:.1e}" for e in eT),
                "monotone decrease", mono, "R+T−1: " + " → ".join(f"{q['RT'] - 1:+.1e}" for q in s))
            if nl == 20:
                last = s[-1]
                worst = max(abs(last["R_vs_disc"]), abs(last["T_vs_disc"]))
                row("L2", f"FDTD vs exact discrete Fresnel, N_pml=60 ({pol}, Δ=λ0/20)", "0", f"{worst:.1e}",
                    "rel err < 1e-6", worst < 1e-6, "proves the code solves the lattice problem exactly up to PML echo")
    brew = brewster_sweep()
    thp = [(b["theta1"], b["R"]) for b in brew if b["pol"] == "p"]
    tmin = min(thp, key=lambda q: q[1])
    tb = math.degrees(math.atan(NMED))
    row("L2 (opt.)", "Brewster: angle of minimum R_p (p-pol sweep, Δ=λ0/20)", f"{tb:.2f}°",
        f"{tmin[0]:.2f}° (R_p = {tmin[1]:.2e})", "dip present", tmin[1] < 0.2 * min(b["R"] for b in brew if b["pol"] == "s"),
        "coarse sweep in Lx; minimum sampled R_p vs smallest R_s")
    figures(res, brew)
    gates = [t for t in T if t["passed"] is not None]
    ok = all(t["passed"] for t in gates)
    with open(os.path.join(ROOT, "results", "level2.json"), "w") as fh:
        json.dump(dict(table=T, runs={f"{k[0]}_{k[1]}_{k[2]}": v for k, v in res.items()}, brewster=brew,
                       pml_study=study, passed=ok),
                  fh, indent=2, default=float)
    print(f"\nLEVEL 2: {sum(t['passed'] for t in gates)}/{len(gates)} gates PASS ->", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def figures(res, brew):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    for pol, mk in (("s", "o"), ("p", "s")):
        for q, ls in (("R", "-"), ("T", "--")):
            e = [abs(res[(nl, pol, "a")][f"{q}_relerr"]) for nl in (10, 20, 40)]
            ax.loglog([0.1, 0.05, 0.025], e, marker=mk, ls=ls, label=f"{q}, {pol} (arith. mean)")
        e = [abs(res[(nl, pol, "s")]["R_relerr"]) for nl in (10, 20)]
        ax.loglog([0.1, 0.05], e, marker=mk, ls=":", color="gray", label=f"R, {pol} (staircase)")
    ax.loglog([0.1, 0.025], [0.02, 0.02 / 16], "k:", label="slope 2")
    ax.axhline(0.01, color="r", lw=0.8, ls="--", label="1%")
    ax.set_xlabel("Δ  [λ0]")
    ax.set_ylabel("|relative error| vs continuous Fresnel")
    ax.set_title("Fresnel convergence, n=1.5, θ1=36.9°")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    for pol in ("s", "p"):
        ax.semilogy([10, 20, 40], [max(abs(res[(nl, pol, "a")]["RT"] - 1), 1e-17) for nl in (10, 20, 40)], "o-",
                    label=f"{pol}-pol")
    ax.axhline(1e-4, color="r", ls="--", lw=0.8, label="1e-4")
    ax.set_xlabel("cells per λ0")
    ax.set_ylabel("|R + T − 1|")
    ax.set_title("Energy balance (discrete-conserved flux)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[2]
    th = np.linspace(0.5, 70, 300)
    c1 = np.cos(np.radians(th))
    c2 = np.sqrt(1 - (np.sin(np.radians(th)) / NMED) ** 2)
    F = fresnel(1.0, NMED, c1, c2)
    ax.plot(th, F["Rs"], "C0-", lw=1, label="Fresnel R_s")
    ax.plot(th, F["Rp"], "C1-", lw=1, label="Fresnel R_p")
    for pol, col in (("s", "C0"), ("p", "C1")):
        b = [q for q in brew if q["pol"] == pol]
        ax.plot([q["theta1"] for q in b], [q["R"] for q in b], "o", color=col, label=f"FDTD R_{pol} (Δ=λ0/20)")
    ax.axvline(math.degrees(math.atan(NMED)), color="k", ls=":", lw=0.8, label="θ_B = 56.31°")
    ax.set_xlabel("θ1  [deg]")
    ax.set_ylabel("R")
    ax.set_title("Angle sweep (m=1, n=0, varying Lx)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L2_fresnel.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
