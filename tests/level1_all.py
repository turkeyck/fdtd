"""Level 1 items 2-8: vacuum self-validation (item 1, leakage, is tests/level1_leakage.py).

Runs (cached by fdtd_io.run), all S = 0.5, SF = 1 lambda0, TF = 10 lambda0, N_pml = 20 unless stated:
  main_nl{10,20,40}_s, main_nl20_p      steady-state DFT over periods [80, 110]
  main_nl20_{s,p}_mneg                  (m, n) = (-1, 1) for the symmetry test
  pml_N{5,10,20,30}_{s,p}, pml_kappa*   PML reflection (TF = 4 lambda0)
  stab_nl20_s                           20000 steps, source switched off
Results: results/level1.json, figures/level1/*.png, figures/level1/L1_wave.gif
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
from fdtd_theory import ECOMP, HCOMP  # noqa: E402

RUNS = os.path.join(ROOT, "runs", "level1")
FIG = os.path.join(ROOT, "figures", "level1")
T = []  # validation table rows


def row(item, quantity, theory, measured, threshold, passed, note=""):
    T.append(dict(item=item, quantity=quantity, theory=theory, measured=measured, threshold=threshold,
                  passed=None if passed is None else bool(passed), note=note))
    flag = "INFO" if passed is None else ("PASS" if passed else "FAIL")
    print(f"[{item}] {quantity}: theory={theory} measured={measured} thr={threshold} -> {flag} {note}")


# ----------------------------------------------------------------------------- runs
def geom(nl, tf_l=10):
    return dict(nl=nl, npml=20, sf=nl, tf=tf_l * nl)


def plane_layout(nl):
    j0 = 20 + nl
    jsf = 20 + nl // 2
    pairs = [j0 + int(round(f * 10 * nl)) for f in (0.05, 0.25, 0.5, 0.75, 0.97)]
    return j0, jsf, pairs


def main_run(nl, pol="s", m=1, tag=""):
    per = 2 * nl
    j0, jsf, pairs = plane_layout(nl)
    planes = [jsf] + [p for q in pairs for p in (q, q + 1)]
    out = os.path.join(RUNS, f"main_nl{nl}_{pol}{tag}")
    meta = fdtd_io.run(out, quiet=False, pol=pol, m=m, n=1, inc="a", nsteps=110 * per, dft0=80 * per,
                       dft1=110 * per, yplanes=planes, zk=nl, xi=nl // 2, div_every=20, energy_every=0,
                       **geom(nl))
    return out, meta


# ----------------------------------------------------------------------------- per-run analysis
def analyze_main(out, meta):
    st = an.setup_from_meta(meta)
    nl, D = meta["nl"], meta["Delta"]
    j0, jsf, pairs = plane_layout(nl)
    jlo, jhi = j0 + 2, meta["Ny"] - meta["npml"] - 1
    B, echo_res, echo_rel = an.echo_from_sf(out, meta, st, jsf)
    names = ["xy", "yz"] + [f"y{p}" for q in pairs for p in (q, q + 1)]
    fits = {}
    for label, Bx in (("raw", None), ("echo_subtracted", B)):
        sets, labels = an.phase_sets(out, meta, st, names, jlo, jhi, B=Bx)
        k, rms, mx, _ = an.wavefront_fit(sets)
        fits[label] = dict(k=k.tolist(), rms=rms, max=mx, npts=int(sum(len(s[3]) for s in sets)),
                           rel_err=(np.abs(k - an.k_forward(st)) / np.abs(an.k_forward(st))).tolist())
    # forward amplitudes, divergence, flux, Poynting on each TF plane pair
    kf = an.k_forward(st)
    Kt = st.Kt
    wt = st.wt
    amps, divs, fluxes, ratio_raw, ratio_int, Sang = [], [], [], [], [], []
    Ei_th, Hi_th, S_th, ratio_int_th = an.interp_theory(st)
    for q in pairs:
        Fa, sla = fdtd_io.load_dft(out, f"y{q}", meta)
        Fb, slb = fdtd_io.load_dft(out, f"y{q + 1}", meta)
        fluxes += [an.conserved_flux(Fa, D), an.conserved_flux(Fb, D)]
        dE, dH = an.div_from_planes(Fa, Fb, D)
        Emax = max(np.abs(Fb[c]).max() for c in ECOMP)
        Hmax = max(np.abs(Fa[c]).max() for c in HCOMP)
        divs.append((float(np.abs(dE).max() / (np.linalg.norm(Kt) * Emax)),
                     float(np.abs(dH).max() / (np.linalg.norm(Kt) * Hmax))))
        Fa_f = an.subtract_echo(Fa, meta, sla, st, B)
        Fb_f = an.subtract_echo(Fb, meta, slb, st, B)
        A = {c: an.mode_amplitude(Fa_f[c], meta, sla, c, kf) for c in ECOMP + HCOMP}
        amps.append(A)
        AE = np.array([A[c] for c in ECOMP])
        AH = np.array([A[c] for c in HCOMP])
        ratio_raw.append(np.linalg.norm(AH) / np.linalg.norm(AE))
        Fc = an.to_cell_centres(Fa_f, Fb_f)
        Ec = np.sqrt(sum(np.abs(Fc[c]) ** 2 for c in ECOMP))
        Hc = np.sqrt(sum(np.abs(Fc[c]) ** 2 for c in HCOMP))
        ratio_int.append(float(np.mean(Hc / Ec)))
        S = an.poynting(Fc).reshape(3, -1).mean(axis=1)
        Sang.append(S)
    AE = np.array([[a[c] for c in ECOMP] for a in amps])
    AH = np.array([[a[c] for c in HCOMP] for a in amps])
    nE = np.linalg.norm(AE, axis=1)
    nH = np.linalg.norm(AH, axis=1)
    EK = np.abs(AE @ Kt) / (nE * np.linalg.norm(Kt))
    HK = np.abs(AH @ Kt) / (nH * np.linalg.norm(Kt))
    EH = np.abs(np.sum(AE * np.conj(AH), axis=1)) / (nE * nH)
    Smean = np.mean(np.array(Sang), axis=0)
    log = fdtd_io.load_log(out)
    dv = ~np.isnan(log["divE_TF"])
    tdiv_E = log["divE_TF"][dv] / (np.linalg.norm(Kt) * np.linalg.norm(st.amplitudes()[0]))
    tdiv_H = log["divH_TF"][dv] / (np.linalg.norm(Kt) * np.linalg.norm(st.amplitudes()[1]))
    return dict(
        nl=nl, Delta=D, dt=meta["dt"], pol=meta["pol"], m=meta["m"], ky_disc=st.ky, ky_cont=st.ky_cont,
        echo_rel=echo_rel, echo_dB=20 * math.log10(echo_rel) if echo_rel > 0 else -math.inf,
        echo_model_residual=echo_res, fits=fits,
        div_dft=[max(d[0] for d in divs), max(d[1] for d in divs)],
        div_time=dict(E_max=float(tdiv_E.max()), H_max=float(tdiv_H.max()), t=log["t"][dv].tolist(),
                      E=tdiv_E.tolist(), H=tdiv_H.tolist()),
        EK=float(EK.max()), HK=float(HK.max()), EH=float(EH.max()),
        ratio_raw=[float(r) for r in ratio_raw], ratio_theory=float(np.linalg.norm(Kt) / wt),
        ratio_int=ratio_int, ratio_int_theory=float(ratio_int_th),
        flux=fluxes, flux_theory=an.theory_flux(st, meta),
        S_meas=Smean.tolist(), S_theory_interp=S_th.tolist(),
        ang_S_k=an.angle_deg(Smean, kf), ang_Sth_k=an.angle_deg(S_th, kf), ang_S_Sth=an.angle_deg(Smean, S_th),
        ang_Kt_k=an.angle_deg(Kt, kf), ang_vg_k=an.angle_deg(an.vg_direction(st), kf),
        jlo=jlo, jhi=jhi, jsf=jsf, pairs=pairs,
    )


# ----------------------------------------------------------------------------- items
def item2_3_4_5(res):
    D20 = {p: res[(20, p)] for p in ("s", "p")}
    for p, r in D20.items():
        tag = f"{p}-pol, Δ=λ0/20"
        fr = r["fits"]["raw"]
        fs = r["fits"]["echo_subtracted"]
        row("L1-2", f"far-PML echo |B|/|A| ({tag})", "quantify", f"{r['echo_rel']:.3e} ({r['echo_dB']:.1f} dB)",
            "—", None, f"single-plane-wave model residual {r['echo_model_residual']:.1e}")
        row("L1-2", f"k fit, raw TF phase ({tag})", "k_disc", f"max rel err {max(fr['rel_err']):.2e}, RMS {fr['rms']:.2e} rad",
            "—", None, "echo NOT removed; shows the echo effect")
        kth = [math.pi, r["ky_disc"], 2 * math.pi / 3]
        for q, ax in enumerate("xyz"):
            row("L1-2", f"k{ax} fit rel. error ({tag}, echo subtracted)", f"{kth[q]:.12f} /λ0",
                f"{fs['rel_err'][q]:.2e}", "< 1e-6", fs["rel_err"][q] < 1e-6,
                f"fit {fs['k'][q]:.12f} /λ0")
        row("L1-2", f"phase residual RMS ({tag})", "0", f"{fs['rms']:.2e} rad", "< 1e-3 rad", fs["rms"] < 1e-3,
            f"{fs['npts']} samples, max {fs['max']:.1e} rad")
    # item 3: convergence
    nls = [10, 20, 40]
    kyf = np.array([res[(nl, "s")]["fits"]["echo_subtracted"]["k"][1] for nl in nls])
    kyd = np.array([res[(nl, "s")]["ky_disc"] for nl in nls])
    kyc = res[(20, "s")]["ky_cont"]
    Ds = np.array([1.0 / nl for nl in nls])
    err_c = np.abs(kyf - kyc)
    err_d = np.abs(kyf - kyd) / kyd
    slope = float(np.polyfit(np.log(Ds), np.log(err_c), 1)[0])
    for nl, ec, ed in zip(nls, err_c, err_d):
        row("L1-3", f"|ky_meas − ky_disc|/ky_disc (Δ=λ0/{nl})", "0", f"{ed:.2e}", "< 1e-6", ed < 1e-6)
        row("L1-3", f"|ky_meas − ky_cont| (Δ=λ0/{nl})", f"{3.013 / nl ** 2:.3e} (leading term)", f"{ec:.4e} /λ0",
            "—", None)
    row("L1-3", "log-log slope of |ky_meas − ky_cont| vs Δ", "2", f"{slope:.4f}", "2 ± 0.1", abs(slope - 2) < 0.1)
    # item 4
    for p, r in D20.items():
        tag = f"{p}-pol, Δ=λ0/20"
        row("L1-4", f"max|∇·E|/(|K̃||E|), DFT, TF interior ({tag})", "0", f"{r['div_dft'][0]:.2e}", "< 1e-10",
            r["div_dft"][0] < 1e-10)
        row("L1-4", f"max|∇·H|/(|K̃||H|), DFT, TF interior ({tag})", "0", f"{r['div_dft'][1]:.2e}", "< 1e-10",
            r["div_dft"][1] < 1e-10)
        tE = np.array(r["div_time"]["E"])
        half = len(tE) // 2
        grow = float(tE[half:].max() / max(tE[:half].max(), 1e-300))
        row("L1-4", f"time-domain max|∇·E| over run ({tag})", "0, not growing",
            f"{r['div_time']['E_max']:.2e} (late/early max ratio {grow:.2f})", "< 1e-10",
            r["div_time"]["E_max"] < 1e-10 and r["div_time"]["H_max"] < 1e-10,
            f"∇·H max {r['div_time']['H_max']:.2e}")
        row("L1-4", f"|E·K̃|, |H·K̃|, |E·H*| normalized ({tag})", "0",
            f"{r['EK']:.1e}, {r['HK']:.1e}, {r['EH']:.1e}", "< 1e-10",
            max(r["EK"], r["HK"], r["EH"]) < 1e-10)
    # item 5
    for p, r in D20.items():
        tag = f"{p}-pol, Δ=λ0/20"
        rr = np.array(r["ratio_raw"])
        err = float(np.abs(rr / r["ratio_theory"] - 1).max())
        row("L1-5", f"|H|/|E| (own Yee points) vs |K̃|/(μ0ω̃) ({tag})", f"{r['ratio_theory']:.15f}",
            f"{rr.mean():.15f}", "rel < 1e-8", err < 1e-8, f"max rel err {err:.1e}; |K̃|/(μ0ω̃) ≡ 1/η0 (derivation §2.5)")
        ri = np.array(r["ratio_int"])
        err_i = float(np.abs(ri / r["ratio_int_theory"] - 1).max())
        row("L1-5", f"|H'|/|E'| cell-centre interpolated ({tag})", f"{r['ratio_int_theory']:.10f}", f"{ri.mean():.10f}",
            "—", None, f"rel err vs interp theory {err_i:.1e}; deviation from 1/η0 = {ri.mean() - 1:+.3e} "
            f"(O((k0Δ)²), (k0Δ)²={(2 * math.pi / 20) ** 2:.3f})")
        fl = np.array(r["flux"])
        var = float((fl.max() - fl.min()) / abs(fl.mean()))
        row("L1-5", f"S_y (conserved flux Φ) variation over 10 TF planes ({tag})", "0", f"{var:.2e}", "< 1e-4",
            var < 1e-4, f"Φ/Φ_inc,theory − 1 = {fl.mean() / r['flux_theory'] - 1:+.2e}")
        row("L1-5", f"angle(S, k) ({tag})", f"{r['ang_Sth_k']:.4f}° (interp. theory)", f"{r['ang_S_k']:.4f}°", "report",
            None, f"angle(S, S_theory) = {r['ang_S_Sth']:.1e}°; angle(K̃,k) = {r['ang_Kt_k']:.4f}°, "
            f"angle(v_g,k) = {r['ang_vg_k']:.4f}°")


def item6():
    rows = []
    for pol in ("s", "p"):
        for N in (5, 10, 20, 30):
            rows.append(pml_run(pol, N))
    for kap, alp in ((1.0, 0.0), (1.0, None), (10.0, None)):
        rows.append(pml_run("s", 20, kappa=kap, alpha=alp))
    for r in rows:
        default = r["kappa_max"] == 5.0 and abs(r["alpha_max"] - 0.1 * math.pi) < 1e-12
        if default:
            row("L1-6", f"PML reflection ({r['pol']}-pol, N_pml={r['N']}, Δ=λ0/20)", "—", f"{r['R_dB']:.1f} dB",
                "< −40 dB", r["R_dB"] < -40.0 if r["N"] == 20 else None,
                "target applies to N_pml=20" if r["N"] != 20 else "")
        else:
            row("L1-6", f"PML reflection (s, N=20, κ_max={r['kappa_max']}, α_max={r['alpha_max']:.3f})", "—",
                f"{r['R_dB']:.1f} dB", "—", None, "sensitivity, not a gate")
    return rows


def pml_run(pol, N, kappa=5.0, alpha=None):
    nl, per = 20, 40
    alpha = 0.1 * math.pi if alpha is None else alpha
    tag = f"pml_N{N}_{pol}" + ("" if (kappa == 5.0 and abs(alpha - 0.1 * math.pi) < 1e-12) else f"_k{kappa}_a{alpha:.3f}")
    out = os.path.join(RUNS, tag)
    jsf = N + 10
    meta = fdtd_io.run(out, pol=pol, inc="a", npml=N, sf=20, tf=80, nl=nl, nsteps=90 * per, dft0=60 * per,
                       dft1=90 * per, yplanes=[jsf], energy_every=0, kappa_max=kappa, alpha_max=alpha)
    st = an.setup_from_meta(meta)
    _, res, rel = an.echo_from_sf(out, meta, st, jsf)
    return dict(pol=pol, N=N, kappa_max=kappa, alpha_max=alpha, R=rel, R_dB=20 * math.log10(rel), model_res=res)


def item7():
    nl, per = 20, 40
    out = os.path.join(RUNS, "stab_nl20_s")
    meta = fdtd_io.run(out, quiet=False, pol="s", inc="a", nsteps=20000, off_t=40, energy_every=10,
                       div_every=200, **geom(nl))
    log = fdtd_io.load_log(out)
    ok = ~np.isnan(log["W_total"])
    t, W, Wp = log["t"][ok], log["W_total"][ok], log["W_phys"][ok]
    t_off = 40 + 20 + 5 * 4  # off ramp: erf centred at off_t + t0, width tau -> 1e-12 after +5 tau
    peak = W.max()
    after = t >= t_off
    # SPEC: "monotonically decays to < 1e-8 x peak" -> monotone over the whole descent, from the energy peak
    # (reached while the off-ramp starts) down to 1e-8 x peak
    ipk = int(np.argmax(W))
    desc = (np.arange(len(W)) >= ipk) & (W >= 1e-8 * peak)
    t_cross = float(t[(np.arange(len(W)) >= ipk) & (W < 1e-8 * peak)][0])
    dW = np.diff(W[desc]) / W[desc][:-1]
    mono = bool(np.all(dW <= 0.0))
    final = float(W[-1] / peak)
    # SPEC: "no growth late" -> over the second half the energy never exceeds its value at the start of that half
    late = t >= t[-1] / 2
    late_max_ratio = float(W[late].max() / W[late][0])
    # strict per-sample monotonicity over the entire tail, including the double-precision floor (INFO)
    r = np.diff(W[after]) / W[after][:-1]
    ninc = int(np.sum(r > 0))
    floor_at_inc = float(W[after][1:][r > 0].max() / peak) if ninc else 0.0
    row("L1-7", "steps run (Δ=λ0/20, S=0.5)", "≥ 20000", f"{meta['nsteps']}", "≥ 20000", meta["nsteps"] >= 20000)
    row("L1-7", f"W_total monotone from its peak (t={t[ipk]:.2f} T0, off-ramp) down to 1e-8·peak",
        "monotone", f"{int(desc.sum())} samples, max ΔW/W = {dW.max():.2e}", "ΔW ≤ 0", mono,
        f"source fully off at t = {t_off} T0; 1e-8·peak crossed at t = {t_cross:.1f} T0")
    row("L1-7", "final W_total / peak (t = 500 T0)", "→ 0", f"{final:.2e}", "< 1e-8", final < 1e-8)
    row("L1-7", "late-time growth: max W over 2nd half / W at its start", "≤ 1", f"{late_max_ratio:.6f}", "≤ 1",
        late_max_ratio <= 1.0)
    row("L1-7", "strict per-sample monotonicity over the whole tail (incl. round-off floor)", "—",
        f"{ninc} increases, max +{r.max():.2e} relative, all at W/peak ≤ {floor_at_inc:.1e}", "—", None,
        "increases occur only at W/peak ~ 4e-28 (fields ~ 2e-14 ≈ 100 ulp of O(1)); no upward trend")
    return dict(t=t.tolist(), W=W.tolist(), Wp=Wp.tolist(), t_off=t_off, final=final, mono=mono)


def mirror(F, meta, pol):
    """Field of the (-m, n) problem predicted from the (m, n) field by x -> -x (derivation.md §5, L1-8)."""
    sig = -1.0 if pol == "s" else 1.0
    sign = {"Ex": -sig, "Ey": sig, "Ez": sig, "Hx": sig, "Hy": -sig, "Hz": -sig}
    Nx = meta["Nx"]
    out = {}
    for c, f in F.items():
        half = meta["offsets"][c][0] == 0.5
        idx = np.arange(Nx)
        src = (Nx - 1 - idx) if half else (Nx - idx) % Nx
        out[c] = sign[c] * f[src]
    return out


def item8(res_runs):
    for pol in ("s", "p"):
        o1, m1 = res_runs[(20, pol)]
        o2, m2 = res_runs[(20, pol, "mneg")]
        worst, scale = 0.0, 0.0
        for name in ["xy"] + [s["name"] for s in m1["slices"] if s["type"] == 0]:
            F1, sl = fdtd_io.load_dft(o1, name, m1)
            F2, _ = fdtd_io.load_dft(o2, name, m2)
            if sl["type"] == 2:
                continue
            P = mirror(F1, m1, pol)
            for c in F1:
                v = fdtd_io.valid_j_mask(m1, sl, c)
                worst = max(worst, float(np.abs(F2[c] - P[c])[v].max()))
                scale = max(scale, float(np.abs(F1[c])[v].max()))
        rel = worst / scale
        row("L1-8", f"(m,n)=(−1,1) vs x-mirror of (1,1), {pol}-pol, Δ=λ0/20", "0", f"{rel:.2e}", "< 1e-6", rel < 1e-6,
            "all six components, x–y slice + all y-planes (incl. SF)")


# ----------------------------------------------------------------------------- figures
def figures(res, runs, stab, pml):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    out, meta = runs[(20, "s")]
    r = res[(20, "s")]
    st = an.setup_from_meta(meta)
    # (a) x-y slice instantaneous field with theoretical equiphase lines
    F, sl = fdtd_io.load_dft(out, "xy", meta)
    comp = "Ez"
    x, y, z = fdtd_io.slice_coords(meta, sl, comp)
    fld = np.real(F[comp])
    fig, ax = plt.subplots(figsize=(11, 3.6))
    im = ax.pcolormesh(y.T, x.T, fld.T, shading="nearest", cmap="RdBu_r", vmin=-1, vmax=1)
    ph = st.kx * x + st.ky * y + st.kz * z
    ax.contour(y.T, x.T, np.cos(ph).T, levels=[0.0], colors="k", linewidths=0.6)
    for yy, lab in ((meta["j0"] * meta["Delta"], "TF/SF"), (meta["npml"] * meta["Delta"], "PML"),
                    ((meta["Ny"] - meta["npml"]) * meta["Delta"], "PML")):
        ax.axvline(yy, color="g", lw=1, ls="--")
        ax.text(yy + 0.05, 0.05, lab, color="g", fontsize=8, ha="left", va="bottom",
                bbox=dict(facecolor="white", alpha=0.7, lw=0))
    ax.set_xlabel("y  [λ0]")
    ax.set_ylabel("x  [λ0]")
    ax.set_title(f"Re Ez (DFT phasor, t=0) on x–y slice z={z[0, 0]:.3f}λ0, s-pol, Δ=λ0/20; black: theoretical "
                 "zero-phase lines cos(k·r)=0 (parallel slanted lines)", fontsize=9)
    fig.colorbar(im, ax=ax, pad=0.01)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_2_xy_field_equiphase.png"), dpi=130)
    plt.close(fig)
    # (b) x-z phase residual map on a TF y-plane (after the fit)
    j0, jsf, pairs = plane_layout(20)
    name = f"y{pairs[2]}"
    B, _, _ = an.echo_from_sf(out, meta, st, jsf)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for axx, (lab, Bx) in zip(axes, (("raw", None), ("echo subtracted", B))):
        Fp, slp = fdtd_io.load_dft(out, name, meta)
        if Bx is not None:
            Fp = an.subtract_echo(Fp, meta, slp, st, Bx)
        key = "raw" if Bx is None else "echo_subtracted"
        k = np.array(r["fits"][key]["k"])
        xx, yy, zz = fdtd_io.slice_coords(meta, slp, comp)
        resid = np.angle(Fp[comp] * np.exp(-1j * (k[0] * xx + k[1] * yy + k[2] * zz)))
        resid -= np.angle(np.mean(np.exp(1j * resid)))
        imx = axx.pcolormesh(xx, zz, resid, shading="nearest", cmap="coolwarm")
        axx.set_xlabel("x  [λ0]")
        axx.set_ylabel("z  [λ0]")
        axx.set_title(f"Ez phase residual, plane y={yy[0, 0]:.3f}λ0 ({lab})\nRMS {np.sqrt(np.mean(resid ** 2)):.1e} rad",
                      fontsize=9)
        fig.colorbar(imx, ax=axx, label="rad")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_2_xz_phase_residual.png"), dpi=130)
    plt.close(fig)
    # (c) dispersion convergence
    nls = np.array([10, 20, 40])
    kyf = np.array([res[(nl, "s")]["fits"]["echo_subtracted"]["k"][1] for nl in nls])
    kyd = np.array([res[(nl, "s")]["ky_disc"] for nl in nls])
    kyc = res[(20, "s")]["ky_cont"]
    fig, ax = plt.subplots(figsize=(6, 4.3))
    ax.loglog(1 / nls, np.abs(kyf - kyc), "o-", label="|ky_meas − ky_cont|")
    ax.loglog(1 / nls, np.maximum(np.abs(kyf - kyd), 1e-16), "s-", label="|ky_meas − ky_disc|")
    ax.loglog(1 / nls, 3.013 / nls ** 2, "k:", label="3.013 Δ² (leading-order theory)")
    ax.axhline(1e-6 * kyd.mean(), color="r", ls="--", lw=0.8, label="1e-6 · ky")
    ax.set_xlabel("Δ  [λ0]")
    ax.set_ylabel("ky error  [1/λ0]")
    ax.set_title("Dispersion convergence (s-pol, S=0.5)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_3_dispersion.png"), dpi=130)
    plt.close(fig)
    # (d) divergence in time
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for p in ("s", "p"):
        dt_ = res[(20, p)]["div_time"]
        ax.semilogy(dt_["t"], np.maximum(dt_["E"], 1e-20), lw=1, label=f"∇·E, {p}")
        ax.semilogy(dt_["t"], np.maximum(dt_["H"], 1e-20), lw=1, ls="--", label=f"∇·H, {p}")
    ax.axhline(1e-10, color="r", ls=":", label="1e-10")
    ax.set_xlabel("t  [λ0/c]")
    ax.set_ylabel("max |∇·F| / (|K̃||F0|), TF interior")
    ax.set_title("Discrete Gauss law in time (Δ=λ0/20)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_4_divergence_time.png"), dpi=130)
    plt.close(fig)
    # (e) flux planes
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for p in ("s", "p"):
        rr = res[(20, p)]
        ys = [q * 0.05 for pr in rr["pairs"] for q in (pr, pr + 1)]
        ax.plot(ys, np.array(rr["flux"]) / rr["flux_theory"] - 1, "o-", label=f"{p}-pol")
    ax.axhspan(-1e-4, 1e-4, color="g", alpha=0.1, label="±1e-4")
    ax.set_xlabel("plane y  [λ0]")
    ax.set_ylabel("Φ(y)/Φ_inc,theory − 1")
    ax.set_title("Conserved Poynting flux across TF planes")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_5_flux.png"), dpi=130)
    plt.close(fig)
    # (f) PML reflection
    fig, ax = plt.subplots(figsize=(6, 4))
    for p in ("s", "p"):
        rs = [q for q in pml if q["pol"] == p and q["kappa_max"] == 5.0 and abs(q["alpha_max"] - 0.1 * math.pi) < 1e-12]
        ax.plot([q["N"] for q in rs], [q["R_dB"] for q in rs], "o-", label=f"{p}-pol")
    ax.axhline(-40, color="r", ls="--", label="target −40 dB")
    ax.set_xlabel("N_pml")
    ax.set_ylabel("reflection  [dB]")
    ax.set_title("CPML reflection vs thickness (Δ=λ0/20, θ=36.9°)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_6_pml.png"), dpi=130)
    plt.close(fig)
    # (g) stability
    fig, ax = plt.subplots(figsize=(6.5, 4))
    W = np.array(stab["W"])
    ax.semilogy(stab["t"], W / W.max(), lw=1, label="W_total")
    ax.semilogy(stab["t"], np.array(stab["Wp"]) / W.max(), lw=1, ls="--", label="W_phys")
    ax.axvline(stab["t_off"], color="k", ls=":", label="source fully off")
    ax.axhline(1e-8, color="r", ls="--", lw=0.8, label="1e-8")
    ax.set_xlabel("t  [λ0/c]")
    ax.set_ylabel("W / peak")
    ax.set_title("20000 steps, Δ=λ0/20, S=0.5")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "L1_7_stability.png"), dpi=130)
    plt.close(fig)


def animation():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    out = os.path.join(RUNS, "anim_nl20_s")
    meta = fdtd_io.run(out, pol="s", inc="a", nsteps=2000, zk=20, snap=10, snapcomp="Ez", energy_every=0,
                       **geom(20))
    S = fdtd_io.load_snapshots(out, meta)
    D = meta["Delta"]
    fig, ax = plt.subplots(figsize=(9, 2.6))
    im = ax.imshow(S[0], origin="lower", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
                   extent=[0, (meta["Ny"] + 1) * D, 0, meta["Lx"]])
    for yy in (meta["j0"] * D, meta["npml"] * D, (meta["Ny"] - meta["npml"]) * D):
        ax.axvline(yy, color="g", lw=0.8, ls="--")
    ax.set_xlabel("y  [λ0]")
    ax.set_ylabel("x  [λ0]")
    ttl = ax.set_title("")

    def upd(q):
        im.set_data(S[q])
        ttl.set_text(f"Ez on x–y slice, t = {(q + 1) * 10 * meta['dt']:.2f} T0 (s-pol, Δ=λ0/20)")
        return im, ttl

    frames = range(0, len(S), 2)
    FuncAnimation(fig, upd, frames=frames, blit=False).save(os.path.join(FIG, "L1_wave.gif"),
                                                            writer=PillowWriter(fps=15), dpi=70)
    plt.close(fig)


def main():
    fdtd_io.build()
    os.makedirs(FIG, exist_ok=True)
    runs = {}
    for nl in (10, 20, 40):
        runs[(nl, "s")] = main_run(nl, "s")
    runs[(20, "p")] = main_run(20, "p")
    for pol in ("s", "p"):
        runs[(20, pol, "mneg")] = main_run(20, pol, m=-1, tag="_mneg")
    res = {key: analyze_main(*runs[key]) for key in runs if len(key) == 2}
    item2_3_4_5(res)
    pml = item6()
    stab = item7()
    item8(runs)
    figures(res, runs, stab, pml)
    animation()
    gates = [t for t in T if t["passed"] is not None]
    ok = all(t["passed"] for t in gates)
    summary = {f"{k[0]}_{k[1]}": {kk: vv for kk, vv in v.items() if kk != "div_time"} for k, v in res.items()}
    with open(os.path.join(ROOT, "results", "level1.json"), "w") as fh:
        json.dump(dict(table=T, runs=summary, pml=pml, passed=ok), fh, indent=2, default=float)
    print(f"\nLEVEL 1 (items 2-8): {sum(t['passed'] for t in gates)}/{len(gates)} gates PASS ->",
          "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
