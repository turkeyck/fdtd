"""Level 3: compare fdtd3d_oblique with Meep (numpy + matplotlib only; reads saved outputs).

Inputs (produced by tests/level3_meep.py):
  ours : runs/level1/main_nl20_{s,p}       (vacuum, TF 10 λ0)   runs/level2/med_nl20_{s,p}_a (medium)
  meep : runs/meep/vac10_{s,p}  vac8_{s,p} (normalization)  med_{s,p}_avg{1,0}
Items (SPEC §10 Level 3):
  (a) vacuum ky            |ky_meep − ky_ours|/ky < 1e-6  (echo-immune 3-term-recurrence estimator)
  (b) wavefront angle and phase residual (3D-style fit on the x–y slice + y-planes)
  (c) Fresnel R, T         |X_meep − X_ours|/X_ours < 0.5 %, each also vs continuous Fresnel
  (d) point-wise DFT field (vacuum, Yee-aligned, reference-point normalized) relative L2 < 1e-3
  (e) PML reflection       recorded only
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
import analyze as an  # noqa: E402
import fdtd_io  # noqa: E402
from fdtd_theory import ECOMP, HCOMP, Setup, fresnel  # noqa: E402

COMPS = ECOMP + HCOMP
NL = 20
D = 1.0 / NL
FIG = os.path.join(ROOT, "figures", "level3")
T = []


def row(item, quantity, ours, meep, theory, threshold, passed, note=""):
    T.append(dict(item=item, quantity=quantity, ours=ours, meep=meep, theory=theory, threshold=threshold,
                  passed=None if passed is None else bool(passed), note=note))
    flag = "INFO" if passed is None else ("PASS" if passed else "FAIL")
    print(f"[{item}] {quantity}: ours={ours} meep={meep} theory={theory} thr={threshold} -> {flag} {note}")


# ----------------------------------------------------------------------------- data access
def ours_xy(pol):
    out = os.path.join(ROOT, "runs", "level1", f"main_nl{NL}_{pol}")
    meta = fdtd_io.load_meta(out)
    F, sl = fdtd_io.load_dft(out, "xy", meta)
    return out, meta, F, sl


def meep_load(tag):
    d = os.path.join(ROOT, "runs", "meep", tag)
    with open(os.path.join(d, "meep_meta.json")) as fh:
        meta = json.load(fh)
    return meta, np.load(os.path.join(d, "meep.npz"))


def meep_xy_on_our_grid(mz, meta_o):
    """Meep x–y slab mapped onto our (i, j) index grid; returns dict comp -> (Nx, Ny+1) complex (nan where absent)."""
    Nx, Ny = meta_o["Nx"], meta_o["Ny"]
    out = {}
    for c in COMPS:
        ox, oy, _, _ = meta_o["offsets"][c]
        arr, xs, ys = mz[f"xy_{c}"], mz[f"xy_{c}_x"], mz[f"xy_{c}_y"]
        G = np.full((Nx, Ny + 1), np.nan + 0j)
        ii = np.rint(xs / D - ox).astype(int)
        jj = np.rint(ys / D - oy).astype(int)
        okx = (ii >= 0) & (ii < Nx) & (np.abs(xs / D - ox - ii) < 1e-6)
        oky = (jj >= 0) & (jj <= Ny) & (np.abs(ys / D - oy - jj) < 1e-6)
        G[np.ix_(ii[okx], jj[oky])] = arr[np.ix_(okx, oky)]
        out[c] = G
    return out


def meep_plane(mz, j, c, meta_o):
    """Meep y-plane j mapped to our (Nx, Nz) grid for comp c."""
    Nx, Nz = meta_o["Nx"], meta_o["Nz"]
    ox, _, oz, _ = meta_o["offsets"][c]
    arr, xs, zs = mz[f"y{j}_{c}"], mz[f"y{j}_{c}_x"], mz[f"y{j}_{c}_z"]
    ii = np.rint(xs / D - ox).astype(int)
    kk = np.rint(zs / D - oz).astype(int)
    okx = (ii >= 0) & (ii < Nx)
    okz = (kk >= 0) & (kk < Nz)
    G = np.full((Nx, Nz), np.nan + 0j)
    G[np.ix_(ii[okx], kk[okz])] = arr[np.ix_(okx, okz)]
    return G


# ----------------------------------------------------------------------------- estimators
def ky_recurrence(G, meta_o, st, jlo, jhi, comps):
    """cos(kyΔ) from u_{j+1}+u_{j-1} = 2cos(kyΔ) u_j, u_j = x-mode projection of each comp on row j.

    Exact for any superposition of the +ky and -ky waves (i.e. immune to PML echo)."""
    num, den = 0.0 + 0j, 0.0
    for c in comps:
        ox = meta_o["offsets"][c][0]
        x = (np.arange(meta_o["Nx"]) + ox) * D
        u = (G[c] * np.exp(-1j * st.kx * x)[:, None]).mean(axis=0)  # (Ny+1,)
        js = np.arange(jlo + 1, jhi)
        num += np.sum((u[js + 1] + u[js - 1]) * np.conj(u[js]))
        den += np.sum(np.abs(u[js]) ** 2)
    return float(math.acos((num / den).real / 2.0) / D)


def fwd_bwd(G, meta_o, st, ky, jlo, jhi, comps):
    """Least-squares split of each comp's x-projected profile into a e^{iky y} + b e^{-iky y}; returns |b|/|a| (E)."""
    a2 = b2 = 0.0
    for c in comps:
        ox, oy = meta_o["offsets"][c][:2]
        x = (np.arange(meta_o["Nx"]) + ox) * D
        u = (G[c] * np.exp(-1j * st.kx * x)[:, None]).mean(axis=0)
        js = np.arange(jlo, jhi + 1)
        y = (js + oy) * D
        A = np.stack([np.exp(1j * ky * y), np.exp(-1j * ky * y)], axis=1)
        sol, *_ = np.linalg.lstsq(A, u[js], rcond=None)
        a2 += abs(sol[0]) ** 2
        b2 += abs(sol[1]) ** 2
    return math.sqrt(b2 / a2)


def fit_xy(G, meta_o, st, jlo, jhi, comps):
    """Plane-wave phase fit on the x–y slice (kx, ky; kz fixed by the slice) with per-comp offsets."""
    sets = []
    for c in comps:
        ox, oy = meta_o["offsets"][c][:2]
        sub = G[c][:, jlo:jhi + 1]
        x = ((np.arange(meta_o["Nx"]) + ox) * D)[:, None] * np.ones_like(sub.real)
        y = ((np.arange(jlo, jhi + 1) + oy) * D)[None, :] * np.ones_like(sub.real)
        ph = an.unwrap2d(np.angle(sub))
        sets.append((x.ravel(), y.ravel(), np.zeros(x.size), ph.ravel()))
    k, rms, mx, _ = an.wavefront_fit(sets)
    return k, rms


def conserved_flux_planes(Fa):
    return an.conserved_flux(Fa, D)


# ----------------------------------------------------------------------------- items
def item_a_b_d_e(summary):
    for pol in ("s", "p"):
        out, meta_o, Fo, sl = ours_xy(pol)
        st = an.setup_from_meta(meta_o)
        mm, mz = meep_load(f"vac10_{pol}")
        Gm = meep_xy_on_our_grid(mz, meta_o)
        e, h = st.amplitudes()
        amp = dict(zip(ECOMP, np.abs(e))) | dict(zip(HCOMP, np.abs(h)))
        comps = [c for c in COMPS if amp[c] > 0.05 * max(amp.values())]
        jlo, jhi = meta_o["j0"] + 2, meta_o["Ny"] - meta_o["npml"] - 1
        # (a)
        ky_o = ky_recurrence(Fo, meta_o, st, jlo, jhi, comps)
        ky_m = ky_recurrence(Gm, meta_o, st, jlo, jhi, comps)
        rel = abs(ky_m - ky_o) / ky_o
        row("L3-a", f"vacuum ky, {pol}-pol, Δ=λ0/20, S=0.5 (recurrence estimator)", f"{ky_o:.12f}", f"{ky_m:.12f}",
            f"{st.ky:.12f} (discrete)", "|Δky|/ky < 1e-6", rel < 1e-6,
            f"ours−theory {abs(ky_o - st.ky) / st.ky:.1e}, meep−theory {abs(ky_m - st.ky) / st.ky:.1e}; "
            f"continuous ky = {st.ky_cont:.6f}")
        # (b)
        ko, rms_o = fit_xy(Fo, meta_o, st, jlo, jhi, comps)
        km, rms_m = fit_xy(Gm, meta_o, st, jlo, jhi, comps)
        kz = st.kz
        th = lambda k: math.degrees(math.atan2(math.hypot(k[0], kz), k[1]))  # noqa: E731
        ph = lambda k: math.degrees(math.atan2(kz, k[0]))  # noqa: E731
        row("L3-b", f"polar angle θ from phase fit, {pol}-pol", f"{th(ko):.8f}°", f"{th(km):.8f}°",
            f"{th([st.kx, st.ky]):.8f}° (discrete k)", "|Δθ| < 1e-4°", abs(th(ko) - th(km)) < 1e-4,
            f"phase residual RMS ours {rms_o:.1e} rad, meep {rms_m:.1e} rad (raw, echo included)")
        row("L3-b", f"Meep k_point check: fitted kx/(2π), {pol}-pol", f"{ko[0] / (2 * math.pi):.10f}",
            f"{km[0] / (2 * math.pi):.10f}", f"{mm['k_point'][0]:.10f} (k_point.x)", "|Δ| < 1e-8",
            abs(km[0] / (2 * math.pi) - mm["k_point"][0]) < 1e-8,
            "confirms k_point in units of 2π/a with Bloch phase exp(2πi k·r)")
        # kz from a Meep y-plane
        jpl = meta_o["j0"] + NL
        Gp = {c: meep_plane(mz, jpl, c, meta_o) for c in comps}
        kzs = []
        for c in comps:
            oz = meta_o["offsets"][c][2]
            z = (np.arange(meta_o["Nz"]) + oz) * D
            ph_ = np.unwrap(np.angle(Gp[c][0, :]))
            kzs.append(np.polyfit(z, ph_, 1)[0])
        kz_m = float(np.mean(kzs))
        row("L3-b", f"Meep k_point check: fitted kz/(2π), {pol}-pol", "—", f"{kz_m / (2 * math.pi):.10f}",
            f"{mm['k_point'][2]:.10f} (k_point.z)", "|Δ| < 1e-8", abs(kz_m / (2 * math.pi) - mm["k_point"][2]) < 1e-8)
        # (d) point-wise comparison, TF region, all significant comps, reference-point normalization
        jref = meta_o["j0"] + NL
        cref = "Ez" if pol == "s" else "Ey"
        scale = Fo[cref][0, jref] / Gm[cref][0, jref]
        num = den = numE = denE = numH = denH = 0.0
        diffmap = np.zeros((meta_o["Nx"], jhi - jlo + 1))
        for c in comps:
            a = Fo[c][:, jlo:jhi + 1]
            b = scale * Gm[c][:, jlo:jhi + 1]
            dd = np.abs(a - b) ** 2
            num += np.nansum(dd)
            den += np.nansum(np.abs(a) ** 2)
            if c in ECOMP:
                numE += np.nansum(dd)
                denE += np.nansum(np.abs(a) ** 2)
            else:
                numH += np.nansum(dd)
                denH += np.nansum(np.abs(a) ** 2)
            diffmap = np.maximum(diffmap, np.sqrt(np.nan_to_num(dd)))
        L2 = math.sqrt(num / den)
        # least-squares scalar (no single reference point) for attribution
        s_ls = sum(np.nansum(Fo[c][:, jlo:jhi + 1] * np.conj(Gm[c][:, jlo:jhi + 1])) for c in comps) / \
            sum(np.nansum(np.abs(Gm[c][:, jlo:jhi + 1]) ** 2) for c in comps)
        L2ls = math.sqrt(sum(np.nansum(np.abs(Fo[c][:, jlo:jhi + 1] - s_ls * Gm[c][:, jlo:jhi + 1]) ** 2)
                             for c in comps) / den)
        row("L3-d", f"point-wise DFT field, vacuum TF region, {pol}-pol (ref. point {cref} at y={jref * D:.2f}λ0)",
            "—", f"rel. L2 = {L2:.2e}", "0", "< 1e-3", L2 < 1e-3,
            f"E-only {math.sqrt(numE / denE):.1e}, H-only {math.sqrt(numH / denH):.1e}; LS-scaled {L2ls:.1e}")
        # (e) echo
        ro = fwd_bwd(Fo, meta_o, st, ky_o, jlo, jhi, comps)
        rm = fwd_bwd(Gm, meta_o, st, ky_m, jlo, jhi, comps)
        row("L3-e", f"PML reflection (N=20 cells, {pol}-pol, θ=36.9°)", f"{20 * math.log10(ro):.1f} dB",
            f"{20 * math.log10(rm):.1f} dB", "—", "record only", None,
            "ours: CPML m=3, κ_max=5, α_max=0.1π; Meep: default PML (quadratic profile, R_asymptotic 1e-15)")
        summary[pol] = dict(ky_ours=ky_o, ky_meep=ky_m, ky_theory=st.ky, L2=L2, L2_ls=L2ls, echo_ours=ro,
                            echo_meep=rm, diffmap=diffmap, jlo=jlo, jhi=jhi)


def item_c(summary):
    res = {}
    for pol in ("s", "p"):
        ours = json.load(open(os.path.join(ROOT, "results", "level2.json")))["runs"][f"{NL}_{pol}_a"]
        st = Setup(n_lambda=NL, pol=pol)
        stm = Setup(n_lambda=NL, pol=pol, n_ref=1.5)
        F = fresnel(1.0, 1.5, st.ky_cont / st.k0, stm.ky_cont / (1.5 * stm.k0))
        RF = F["Rs"] if pol == "s" else F["Rp"]
        _, zn = meep_load(f"vac8_{pol}")
        out = {}
        for avg in (1, 0):
            mm, zm = meep_load(f"med_{pol}_avg{avg}")
            inc = float(zn["flux_R"][0])
            R_native = -float(zm["flux_R"][0]) / inc
            T_native = float(zm["flux_T"][0]) / inc
            # same conserved-flux definition as ours, on Meep's y-planes (E at j, H at j+1/2)
            meta_o = fdtd_io.load_meta(os.path.join(ROOT, "runs", "level2", f"med_nl{NL}_{pol}_a"))
            jr, jt = mm["j0"] + NL, mm["j1"] + 2 * NL
            Fn = {c: meep_plane(zn, jr, c, meta_o) for c in COMPS}
            Fr = {c: meep_plane(zm, jr, c, meta_o) - Fn[c] for c in COMPS}
            Ft = {c: meep_plane(zm, jt, c, meta_o) for c in COMPS}
            Pinc = an.conserved_flux(Fn, D)
            R_phi = -an.conserved_flux(Fr, D) / Pinc
            T_phi = an.conserved_flux(Ft, D) / Pinc
            out[avg] = dict(R_native=R_native, T_native=T_native, R_phi=R_phi, T_phi=T_phi)
        res[pol] = dict(ours=ours, meep=out, R_fresnel=RF)
        m1 = out[1]
        for q, qo in (("R", ours["R"]), ("T", ours["T"])):
            for kind in ("native", "phi"):
                v = m1[f"{q}_{kind}"]
                d = abs(v - qo) / qo
                row("L3-c", f"{q}, {pol}-pol, Meep eps_averaging=True, Meep flux = {kind}", f"{qo:.6f}", f"{v:.6f}",
                    f"{RF if q == 'R' else 1 - RF:.6f} (Fresnel)", "rel diff < 0.5%",
                    d < 0.005 if kind == "native" else None,
                    f"rel diff {d:.2e}; Meep vs Fresnel {(v - (RF if q == 'R' else 1 - RF)) / (RF if q == 'R' else 1 - RF):+.2e}; "
                    f"ours vs Fresnel {(qo - (RF if q == 'R' else 1 - RF)) / (RF if q == 'R' else 1 - RF):+.2e}")
        m0 = out[0]
        row("L3-c", f"R, T with Meep eps_averaging=False, {pol}-pol", f"{ours['R']:.6f}, {ours['T']:.6f}",
            f"{m0['R_native']:.6f}, {m0['T_native']:.6f}", "—", "attribution", None,
            f"R rel diff {abs(m0['R_native'] - ours['R']) / ours['R']:.2e}: interface ε treatment differs "
            "(no averaging on the tangential interface plane)")
        row("L3-c", f"Meep R+T−1, {pol}-pol (native add_flux / conserved Φ)", f"{ours['RT'] - 1:+.1e}",
            f"{m1['R_native'] + m1['T_native'] - 1:+.1e} / {m1['R_phi'] + m1['T_phi'] - 1:+.1e}", "0", "record", None)
    summary["c"] = res


def figures(summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6))
    for ax, pol in zip(axes, ("s", "p")):
        s = summary[pol]
        im = ax.imshow(np.log10(np.maximum(s["diffmap"], 1e-12)).T, origin="lower", aspect="auto", cmap="viridis",
                       extent=[0, 2, s["jlo"] * D, s["jhi"] * D])
        ax.set_xlabel("x  [λ0]")
        ax.set_ylabel("y  [λ0]")
        ax.set_title(f"log10 |F_ours − c·F_meep| (max over comps), {pol}-pol\nrel. L2 = {s['L2']:.1e}", fontsize=9)
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L3_d_fielddiff.png"), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    labels, vo, vm1, vm0, vf = [], [], [], [], []
    for pol in ("s", "p"):
        c = summary["c"][pol]
        for q in ("R", "T"):
            labels.append(f"{q}_{pol}")
            vo.append(c["ours"][q])
            vm1.append(c["meep"][1][f"{q}_native"])
            vm0.append(c["meep"][0][f"{q}_native"])
            vf.append(c["R_fresnel"] if q == "R" else 1 - c["R_fresnel"])
    x = np.arange(len(labels))
    for off, v, lab in ((-0.3, vf, "Fresnel"), (-0.1, vo, "ours"), (0.1, vm1, "Meep (avg)"), (0.3, vm0, "Meep (no avg)")):
        ax.bar(x + off, np.array(v) / np.array(vf) - 1, 0.2, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.axhspan(-0.005, 0.005, color="g", alpha=0.1)
    ax.set_ylabel("relative deviation from Fresnel")
    ax.set_title("Level 3 (c): R, T at Δ=λ0/20, n=1.5, θ1=36.9°")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L3_c_RT.png"), dpi=130)
    plt.close(fig)


def main():
    summary = {}
    item_a_b_d_e(summary)
    item_c(summary)
    figures(summary)
    gates = [t for t in T if t["passed"] is not None]
    ok = all(t["passed"] for t in gates)
    clean = {k: ({kk: vv for kk, vv in v.items() if kk != "diffmap"} if isinstance(v, dict) else v)
             for k, v in summary.items()}
    with open(os.path.join(ROOT, "results", "level3.json"), "w") as fh:
        json.dump(dict(table=T, summary=clean, passed=ok), fh, indent=2, default=float)
    print(f"\nLEVEL 3: {sum(t['passed'] for t in gates)}/{len(gates)} gates PASS ->", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
