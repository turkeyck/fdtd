"""Gate D: dispersive materials (Lorentz ADE) and subpixel averaging (CMT-like) of planar layers.

Thresholds: tests/thresholds_disp.json. References: disp_tmm.lattice_rt (exact Yee lattice with the ADE's
discrete-time permittivity and the subpixel/staircase weights) and disp_tmm.stack_cont (continuous TMM with the
fitted model or the raw data, i.e. what Lumerical STACK computes).

    python3 tests/disp_gate.py          (FDTD_FRESH=1 to re-run every simulation)
"""
import concurrent.futures as cf
import json
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ["OMP_NUM_THREADS"] = "1"
import nu_common as nc  # noqa: E402
import pulse_analysis as pa  # noqa: E402
import disp_tmm as D  # noqa: E402

TH = json.load(open(os.path.join(nc.ROOT, "tests", "thresholds_disp.json"), encoding="utf-8"))
CM = TH["common"]
L0 = CM["lambda0_um"]
RUNS = os.path.join(nc.ROOT, "runs", "disp")
FIGS = os.path.join(nc.ROOT, "figures", "disp")
PU = CM["pulse"]
KX = CM["k"]["kxn"]
Y0 = CM["layer_start_lambda0"]
NM = lambda v: v / (L0 * 1000.0)        # noqa: E731  nm -> lambda0 units

STRUCT = {
    "A": ([("Si3N4", NM(75))], "Si"),
    "B": ([("GaAs", NM(50))], "SiO2"),
    "C": ([("SiO2", NM(103)), ("Si3N4", NM(73.5))] * 3, "Si"),
}


def stack_file(name, layers=None, sub=None):
    layers, sub = (layers, sub) if layers is not None else STRUCT[name]
    y, out = Y0, []
    for m, d in layers:
        out.append((m, y, y + d))
        y += d
    out.append((sub, y, 1000.0))
    return D.write_stack(os.path.join(RUNS, "stacks", f"{name}.json"), L0, out)


def args(nl, pol, stk, smode="s", nsteps=None, **extra):
    a = dict(bloch=1, Lx=repr(4 / nl), Lz=repr(4 / nl), kxn=KX, kzn=0, nl=nl, npml=40, sf=40, tf=int(3 * nl),
             pol=pol, src="g", fc=PU["fc"], fw=PU["fw"], nfreq=PU["nfreq"], fmin=PU["fmin"], fmax=PU["fmax"],
             energy_every=0, stack=stk, smode=smode,
             nsteps=nsteps if nsteps is not None else int(round(CM["t_end_T0"] * nl / CM["S"])))
    a.update(extra)
    return a


def run(name, **kw):
    out = os.path.join(RUNS, name)
    meta, _, err = nc.run(out, **kw)
    return out, meta, err


def par(jobs):
    with cf.ThreadPoolExecutor(max_workers=os.cpu_count() or 2) as ex:
        fut = {k: ex.submit(run, k, **kw) for k, kw in jobs.items()}
        return {k: f.result() for k, f in fut.items()}


def geom_only(name, nl, stk, smode):
    """Solver with nsteps=0: writes meta/grid/stack_used only (for the lattice reference at many meshes)."""
    out = os.path.join(RUNS, "geom", f"{name}_{smode}{nl}")
    os.makedirs(out, exist_ok=True)
    kw = args(nl, "s", stk, smode, nsteps=0)
    argv = [f"{k}={v}" for k, v in kw.items()]
    subprocess.run([nc.BIN] + argv + [f"out={out}"], check=True, capture_output=True)
    return out


def band(f):
    return np.exp(-0.5 * ((f - PU["fc"]) / PU["fw"]) ** 2) >= 1e-2


def main():
    nc.build()
    os.makedirs(FIGS, exist_ok=True)
    T = nc.Table("disp")
    gi = lambda nl: f"Δ={L0 * 1000 / nl:.2f} nm (λ0/{nl}), Δt={CM['S'] / nl:.5g} T0, Bloch 4x4, kx/k0={KX}"  # noqa: E731

    # ------------------------------------------------------------ D0-1 non-dispersive parity with the legacy medium
    half = os.path.join(RUNS, "stacks", "half.json")
    os.makedirs(os.path.dirname(half), exist_ok=True)
    json.dump({"lambda0_um": L0, "nmat": 1, "mat1": {"eps_inf": 2.25, "de": [], "nu0": [], "g": []},
               "layer_mat": [1], "layer_y0": [3.0], "layer_y1": [1000.0]}, open(half, "w"))
    base = dict(pol="p", nsteps=600, dump=1, Lx=2.0, Lz=3.0, tf=120, energy_every=0)
    r = par({"d01_legacy": dict(base, eps2=2.25, y1=60), "d01_stack": dict(base, stack=half)})
    A = np.fromfile(os.path.join(r["d01_legacy"][0], "fields_final.bin"))
    B = np.fromfile(os.path.join(r["d01_stack"][0], "fields_final.bin"))
    d = float(np.abs(A - B).max() / np.abs(A).max())
    thr = TH["D0-1"]["max_rel"]
    T.row("D0-1", "non-dispersive stack ε=2.25 half-space (interface on a primal node) vs legacy eps2/y1, p, 600 steps",
          "0", f"{d:.1e}", f"< {thr:g}", d < thr, grid="Δ=λ0/20, 40x60 cell", note="roundoff of the layer coordinate")

    # ------------------------------------------------------------ D1 FDTD vs exact lattice
    nl = CM["nl_fdtd"]
    stk = {k: stack_file(k) for k in STRUCT}
    jobs = {}
    for k in STRUCT:
        for pol in ("s", "p"):
            jobs[f"d1_{k}_{pol}_s"] = args(nl, pol, stk[k])
    jobs["d1_A_s_c"] = args(nl, "s", stk["A"], smode="c")
    jobs["d1_A_p_c"] = args(nl, "p", stk["A"], smode="c")
    jobs["d1_A_s_s_nl80"] = args(80, "s", stk["A"])
    R = par(jobs)
    detail = {}
    for name, (out, meta, _) in R.items():
        me, th = pa.measure(out), pa.theory(out)
        ok = band(me["f"]) & np.isfinite(th["R"])
        dR = float(np.abs(me["R"] - th["R"])[ok].max())
        parts = name.split("_")
        lab = f"structure {parts[1]} ({TH['common']['structures'][parts[1]]}), {parts[2]}, " \
              f"{'subpixel' if parts[3] == 's' else 'staircase'}"
        nlr = 80 if name.endswith("nl80") else nl
        thr = TH["D1-1"]["abs"]
        T.row("D1-1", f"{lab}: max |R−R_lat| over 400–800 nm ({ok.sum()} freqs)", "0", f"{dR:.2e}", f"< {thr:g}",
              dR < thr, grid=gi(nlr))
        if "T" in me and np.isfinite(th["T"]).any():
            okT = ok & np.isfinite(th["T"])
            dT = float(np.abs(me["T"] - th["T"])[okT].max())
            thr = TH["D1-2"]["abs"]
            T.row("D1-2", f"{lab}: max |T−T_lat|", "0", f"{dT:.2e}", f"< {thr:g}", dT < thr, grid=gi(nlr))
        detail[name] = dict(f=me["f"], R=me["R"], T=me.get("T"), R_lat=th["R"], T_lat=th["T"])

    # ------------------------------------------------------------ D2 convergence vs continuous TMM (model)
    f = np.linspace(PU["fmin"], PU["fmax"], PU["nfreq"])
    layA, subA = STRUCT["A"]
    Rc = D.stack_cont(layA, subA, f, L0, 2 * math.pi * KX, "s")[0]
    conv = {}
    for mode in ("s", "c"):
        nls = [10, 14, 20, 28, 40, 56, 80, 113, 160, 226, 320, 452, 640]
        errs = []
        for n in nls:
            Rl, _ = D.lattice_rt(geom_only("A", n, stk["A"], mode), f)
            errs.append(float(np.sqrt(np.mean((Rl - Rc) ** 2))))
        conv[mode] = (nls, errs)
    sel = TH["D2-1"]["nl"]
    for mode, item in (("s", "D2-1"), ("c", "D2-2")):
        nls, errs = conv[mode]
        idx = [nls.index(n) for n in sel]
        order = nc.order_fit([L0 * 1000 / nls[q] for q in idx], [errs[q] for q in idx])
        if mode == "s":
            w = TH["D2-1"]["order_window"]
            T.row("D2-1", f"subpixel: order of RMS_band |R−R_cont| (A, s) over nl={sel}", "2", f"{order:.2f}",
                  f"in {w}", w[0] <= order <= w[1], grid=f"Δ = {L0 * 1000 / sel[0]:.2f} … {L0 * 1000 / sel[-1]:.2f} nm",
                  note="errors " + ", ".join(f"{errs[q]:.1e}" for q in idx))
        else:
            T.row("D2-2", f"staircase: same metric", "1 (expected)", f"{order:.2f}", "INFO", None,
                  grid=f"Δ = {L0 * 1000 / sel[0]:.2f} … {L0 * 1000 / sel[-1]:.2f} nm", note="errors " + ", ".join(f"{errs[q]:.1e}" for q in idx))

    # ------------------------------------------------------------ D3 thickness sweep at fixed mesh (CMT demonstration)
    ds_nm = np.arange(48.0, 102.01, 2.5)
    sweep = {}
    for n in (40, 80):
        for mode in ("s", "c"):
            Rl, Rcs = [], []
            for dn in ds_nm:
                sf = stack_file(f"sweep_{dn:.1f}", [("Si3N4", NM(dn))], "Si")
                Rl.append(D.lattice_rt(geom_only(f"sw{dn:.1f}", n, sf, mode), [1.0])[0][0])
                Rcs.append(D.stack_cont([("Si3N4", NM(dn))], "Si", [1.0], L0, 2 * math.pi * KX, "s")[0][0])
            sweep[(n, mode)] = (np.array(Rl), np.array(Rcs))
        es = np.abs(sweep[(n, "s")][0] - sweep[(n, "s")][1]).max()
        ec = np.abs(sweep[(n, "c")][0] - sweep[(n, "c")][1]).max()
        T.row("D3", f"Si3N4 48–102 nm on Si, R at 600 nm, Δ={L0 * 1000 / n:.1f} nm: max|R−R_cont| subpixel vs staircase",
              "INFO", f"{es:.2e} vs {ec:.2e} (×{ec / es:.0f})", "INFO", None, grid=gi(n))

    # ------------------------------------------------------------ D4 mesh equivalence (Lumerical CMT claim analog)
    # error metric over arbitrary film thickness and the whole band: RMS over Si3N4 48-102 nm (5 nm) x 11 freqs
    f11 = np.linspace(PU["fmin"], PU["fmax"], 11)
    ds4 = np.arange(48.0, 102.01, 5.0)
    Rc4 = {dn: D.stack_cont([("Si3N4", NM(dn))], "Si", f11, L0, 2 * math.pi * KX, "s")[0] for dn in ds4}
    nl4 = [20, 28, 40, 56, 80, 113, 160, 226, 320, 452, 640]
    e4 = {}
    for mode in ("s", "c"):
        e4[mode] = []
        for n in nl4:
            sq = []
            for dn in ds4:
                sf = stack_file(f"sweep_{dn:.1f}", [("Si3N4", NM(dn))], "Si")
                Rl, _ = D.lattice_rt(geom_only(f"sw{dn:.1f}", n, sf, mode), f11)
                sq.append(np.mean((Rl - Rc4[dn]) ** 2))
            e4[mode].append(float(np.sqrt(np.mean(sq))))
    n_si = D.eps_cont_model("Si", [1.5], L0)[0] ** 0.5
    n_film = D.eps_cont_model("Si3N4", [1.5], L0)[0] ** 0.5
    lc, le = np.log(nl4), np.log(e4["c"])
    for ref in (20, 40, 80):
        e_ref = e4["s"][nl4.index(ref)]
        # staircase error is not monotone in nl: first nl beyond which it stays below e_ref, log-log interpolated
        k = next((q for q in range(len(nl4)) if all(e <= e_ref for e in e4["c"][q:])), None)
        if k is None:
            meas = f"not reached by λ0/{nl4[-1]}"
        elif k == 0:
            meas = f"≤ λ0/{nl4[0]}"
        else:
            t = (math.log(e_ref) - le[k - 1]) / (le[k] - le[k - 1])
            need = math.exp(lc[k - 1] + t * (lc[k] - lc[k - 1]))
            ppw = 400.0 / (L0 * 1000 / need) / n_film.real
            meas = f"Δ ≈ {L0 * 1000 / need:.2f} nm (×{need / ref:.1f} finer, {ppw:.1f} pts/λ in Si3N4 at 400 nm)"
        ppw_ref = 400.0 / (L0 * 1000 / ref) / n_film.real
        T.row("D4", f"staircase mesh matching the subpixel error at Δ={L0 * 1000 / ref:.1f} nm "
              f"({ppw_ref:.1f} pts/λ in Si3N4, {400.0 / (L0 * 1000 / ref) / n_si.real:.1f} in Si at 400 nm); "
              f"RMS over thickness 48–102 nm × band = {e_ref:.2e}", "INFO", meas, "INFO", None)
    conv["thk_s"], conv["thk_c"] = (nl4, e4["s"]), (nl4, e4["c"])

    # ------------------------------------------------------------ D5 material fits vs data
    md = D.models()
    for nm_, e in md.items():
        T.row("D5", f"{nm_}: Lorentz fit ({len(e['poles'])} poles) vs {e['source']}, 400–800 nm: RMS rel ε error",
              "INFO", f"{e['fit_rms_rel_eps']:.2e} (max |Δn| {e['fit_max_abs_dn']:.1e}, max |Δk| {e['fit_max_abs_dk']:.1e})",
              "INFO", None)
    for k, (lay, sub) in STRUCT.items():
        for pol in ("s", "p"):
            Rm = D.stack_cont(lay, sub, f, L0, 2 * math.pi * KX, pol)[0]
            Rd = D.stack_cont(lay, sub, f, L0, 2 * math.pi * KX, pol, source="data")[0]
            T.row("D5", f"structure {k}, {pol}: max |R_TMM(model) − R_TMM(data)| over 400–800 nm", "INFO",
                  f"{np.abs(Rm - Rd).max():.2e}", "INFO", None)
            detail[f"tmm_{k}_{pol}"] = dict(f=f, R_model=Rm, R_data=Rd)

    figures(detail, conv, sweep, ds_nm, f)
    T.save(extra=dict(thresholds_status=TH["_status"], detail=detail,
                      convergence={m: dict(nl=v[0], rms=v[1]) for m, v in conv.items()}), name="disp_gate")
    return 0 if T.ok() else 1


def figures(detail, conv, sweep, ds_nm, f):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lam = L0 * 1000 / f
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, k in zip(axs, "ABC"):
        for pol, c in (("s", "C0"), ("p", "C3")):
            d = detail[f"d1_{k}_{pol}_s"]
            t = detail[f"tmm_{k}_{pol}"]
            ax.plot(lam, d["R"], "o", ms=3, color=c, label=f"FDTD {pol} (λ0/40, subpixel)")
            ax.plot(lam, t["R_model"], "-", lw=1, color=c, label=f"TMM {pol}, fitted model")
            ax.plot(lam, t["R_data"], "--", lw=1, color=c, alpha=0.7, label=f"TMM {pol}, raw data")
        ax.set_title(f"{k}: {TH['common']['structures'][k]}", fontsize=8)
        ax.set_xlabel("λ (nm)"); ax.set_ylabel("R")
    axs[0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "D1_D5_spectra.png"), dpi=130); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.2))
    for m, lab in (("s", "subpixel, A (75 nm)"), ("c", "staircase, A (75 nm)"),
                   ("thk_s", "subpixel, RMS over 48–102 nm"), ("thk_c", "staircase, RMS over 48–102 nm")):
        nls, e = conv[m]
        ax.loglog([L0 * 1000 / n for n in nls], e, "o-" if "thk" not in m else "s--", ms=3, label=lab)
    xx = np.array([L0 * 1000 / 640, L0 * 1000 / 160])
    ax.loglog(xx, conv["s"][1][-1] * (xx / xx[0]) ** 2, "k:", lw=0.8, label="slope 2")
    ax.set_xlabel("Δ (nm)"); ax.set_ylabel("RMS over 400–800 nm of |R − R_TMM|")
    ax.set_title("Si3N4 film on Si, s, kx/k0 = 0.2, 400–800 nm", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "D2_convergence.png"), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, n in zip(axs, (40, 80)):
        ax.plot(ds_nm, sweep[(n, "s")][1], "k-", lw=1.2, label="TMM (continuous)")
        ax.plot(ds_nm, sweep[(n, "s")][0], "o-", ms=3, label="subpixel (CMT-like)")
        ax.plot(ds_nm, sweep[(n, "c")][0], "s-", ms=3, label="staircase")
        ax.set_title(f"R at 600 nm vs Si3N4 thickness, Δ = {L0 * 1000 / n:.1f} nm", fontsize=9)
        ax.set_xlabel("Si3N4 thickness (nm)"); ax.set_ylabel("R"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "D3_thickness_sweep.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
