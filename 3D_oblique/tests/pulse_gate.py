"""Gate P: Bloch periodic boundaries (bloch=1) and broadband Gaussian pulse (src=g, nfreq>0).

Thresholds: tests/thresholds_pulse.json (separate from the frozen thresholds_nu.json). Theory: tmm.tmm_discrete,
the exact Yee-lattice R(f), T(f) for the run's y grid, dt and fixed tangential K~ (pulse_analysis.theory).

P0-1  bloch=1 with an integer Bloch phase reproduces the legacy real solver (copy 0), fields at several steps
P0-2  Bloch energy conservation (random init, PEC y walls, non-commensurate kx, kz)
P0-3  pulse spectra of a 4x4 Bloch cell == legacy real supercell with the same (kx, kz)
P1    half-space, non-commensurate (kx, kz), s and p: R(f), T(f) vs exact discrete theory; R + T = 1
P2    normal incidence (kx = kz = 0): R(f), T(f) vs theory
P3    nonuniform-y film F1 on a 4x4 Bloch cell: R(f), T(f) vs theory
P4    cutoff-leakage control (INFO): pulse with energy at f_cut must trigger the solver warning

    python3 tests/pulse_gate.py          (FDTD_FRESH=1 to re-run every simulation)
"""
import concurrent.futures as cf
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ["OMP_NUM_THREADS"] = "1"          # small x-z cells; the gate runs several solvers in parallel instead
import nu_common as nc  # noqa: E402
import grid_gen  # noqa: E402
import pulse_analysis as pa  # noqa: E402

TH = json.load(open(os.path.join(nc.ROOT, "tests", "thresholds_pulse.json"), encoding="utf-8"))
CM = TH["common"]
RUNS = os.path.join(nc.ROOT, "runs", "pulse")
FIGS = os.path.join(nc.ROOT, "figures", "pulse")
KN = CM["k_noncommensurate"]
PU = CM["pulse"]
DT_U = CM["S"] / CM["nl"]                    # uniform mesh: dt = S * Delta
NST_U = int(round(CM["t_end_T0"] / DT_U))
CELL = dict(bloch=1, Lx=0.2, Lz=0.2)
PULSE = dict(src="g", fc=PU["fc"], fw=PU["fw"], nfreq=PU["nfreq"], fmin=PU["fmin"], fmax=PU["fmax"], energy_every=0)
HALF = dict(eps2=2.25, y1=60, tf=120, npml=40)


def run(name, **kw):
    out = os.path.join(RUNS, name)
    meta, _, err = nc.run(out, **kw)
    return out, meta, err


def par(jobs):
    """Run {name: kwargs} concurrently (one thread per solver)."""
    with cf.ThreadPoolExecutor(max_workers=os.cpu_count() or 2) as ex:
        fut = {k: ex.submit(run, k, **kw) for k, kw in jobs.items()}
        return {k: f.result() for k, f in fut.items()}


def band(f):
    a = np.exp(-0.5 * ((f - PU["fc"]) / PU["fw"]) ** 2)
    return a >= 1e-2


def compare(out):
    me, th = pa.measure(out), pa.theory(out)
    ok = band(me["f"]) & np.isfinite(th["R"])
    dR = np.abs(me["R"] - th["R"])[ok]
    dT = np.abs(me["T"] - th["T"])[ok] if "T" in me else np.array([np.nan])
    cons = np.abs(me["R"] + me["T"] - 1)[ok] if "T" in me else np.array([np.nan])
    return me, th, ok, float(dR.max()), float(dT.max()), float(cons.max())


def film_grid():
    """L2_F1_ppw40 (vacuum | n=2.0, 0.08 lambda0 | n=1.46 substrate) with a 4x4 x-z cell for Bloch runs."""
    g = grid_gen.level2("F1", ppw=40)
    g["x"] = grid_gen.uniform_axis(0.2, 0.05)
    g["z"] = grid_gen.uniform_axis(0.2, 0.05)
    g["name"] = "P_F1_ppw40_cell4"
    g["generator"]["note"] = "gate P: L2_F1_ppw40 y grid on a 4x4 Bloch cell"
    g["hash"] = grid_gen.canonical_hash(g)
    path = nc.grid_path(g["name"])
    if not os.path.exists(path) or json.load(open(path)).get("hash") != g["hash"]:
        grid_gen.write(g)
    return g["name"]


def main():
    nc.build()
    os.makedirs(FIGS, exist_ok=True)
    T = nc.Table("pulse")
    gi_u = f"Δ=λ0/{CM['nl']} uniform, Δt={DT_U:g} T0, Bloch cell 4x4"

    # ---------------------------------------------------------------- P0-1 parity at an integer Bloch phase
    steps = TH["P0-1"]["steps"]
    base = dict(pol="p", eps2=2.25, y1=60, tf=120, nsteps=max(steps), dump_at=steps, energy_every=0, Lx=2.0, Lz=3.0)
    r = par({"p01_legacy": base, "p01_bloch": dict(base, bloch=1)})
    (oa, ma, _), (ob, mb, _) = r["p01_legacy"], r["p01_bloch"]
    worst = 0.0
    for N in steps:
        A = nc.load_fields(oa, f"fields_n{N}.bin", ma)
        B = nc.load_fields(ob, f"fields_n{N}.bin", mb)
        for c in nc.COMPS:
            s = max(np.abs(A[c]).max(), 1e-300)
            worst = max(worst, float(np.abs(A[c] - B[c]).max() / s))
    thr = TH["P0-1"]["max_rel"]
    T.row("P0-1", f"bloch=1, kx=2πm/Lx: copy 0 vs legacy, max|ΔF|/max|F| at n={steps}", "0", f"{worst:.1e}",
          f"< {thr:g}", worst < thr, grid="Δ=λ0/20 uniform, 40x60 cell", note="bit-identical" if worst == 0 else "")

    # ---------------------------------------------------------------- P0-2 energy
    n2 = TH["P0-2"]["nsteps"]
    o, m, _ = run("p02_energy", init="r", inc="0", npml=0, sf=20, tf=60, nsteps=n2, energy_every=50,
                  kxn=KN["kxn"], kzn=KN["kzn"], **CELL)
    L = nc.load_log(o)
    W = L["W_mod"][~np.isnan(L["W_mod"])]
    nW = L["n"][~np.isnan(L["W_mod"])]
    drift = float(np.max(np.abs(W / W[0] - 1)))
    thr = TH["P0-2"]["rel_drift"]
    T.row("P0-2", f"Bloch |F|² energy W_mod drift over {n2} steps (PEC y, kx/k0={KN['kxn']}, kz/k0={KN['kzn']}, "
          f"phase e^(ikxLx)={m['bloch_phase_x'][0]:.4f}{m['bloch_phase_x'][1]:+.4f}i)", "0", f"{drift:.2e}",
          f"< {thr:g}", drift < thr, grid=gi_u)

    # ---------------------------------------------------------------- simulations of P0-3 .. P4
    fg = film_grid()
    geo = nc.Geo(nc.load_grid(fg))
    dtf = geo.dt_run()
    nst_f = int(math.ceil(CM["t_end_T0"] / dtf))
    jobs = {}
    for pol in ("s", "p"):
        jobs[f"p03_super_{pol}"] = dict(pol=pol, Lx=2.0, Lz=3.0, m=1, n=1, nsteps=NST_U, **HALF, **PULSE)
        jobs[f"p03_cell_{pol}"] = dict(pol=pol, kxn=0.5, kzn=1.0 / 3.0, nsteps=NST_U, **CELL, **HALF, **PULSE)
        jobs[f"p1_half_{pol}"] = dict(pol=pol, kxn=KN["kxn"], kzn=KN["kzn"], nsteps=NST_U, **CELL, **HALF, **PULSE)
        jobs[f"p2_normal_{pol}"] = dict(pol=pol, kxn=0.0, kzn=0.0, nsteps=NST_U, **CELL, **HALF, **PULSE)
        jobs[f"p3_film_{pol}"] = dict(grid=fg, pol=pol, kxn=KN["kxn"], kzn=KN["kzn"], dt=dtf, nsteps=nst_f,
                                      **CELL, **PULSE)
    jobs["p4_cutoff"] = dict(pol="s", kxn=0.5, kzn=1.0 / 3.0, nsteps=NST_U, **CELL, **HALF,
                             **dict(PULSE, fw=0.15, fmin=0.7, fmax=1.3, nfreq=61))
    R = par(jobs)

    # ---------------------------------------------------------------- P0-3 Bloch cell vs supercell
    thr = TH["P0-3"]["max_rel"]
    for pol in ("s", "p"):
        a, b = pa.load(R[f"p03_super_{pol}"][0]), pa.load(R[f"p03_cell_{pol}"][0])
        d = max(float(np.abs(a[k] - b[k]).max() / np.abs(a[k]).max()) for k in ("proj", "aux"))
        T.row("P0-3", f"pulse spectra, Bloch 4x4 cell vs real 40x60 supercell (m=n=1), {pol}", "0", f"{d:.1e}",
              f"< {thr:g}", d < thr, grid=gi_u)

    # ---------------------------------------------------------------- P1 .. P3 vs exact discrete theory
    detail = {}

    def judge(item, label, out, gi, keys=("R", "T", "RT")):
        me, th, ok, dR, dT, cons = compare(out)
        meta = pa.load(out)["meta"]
        detail[os.path.basename(out)] = dict(f=me["f"], R=me["R"], T=me.get("T"), R_th=th["R"], T_th=th["T"],
                                             judged=ok, theta_deg=pa.angle_deg(meta, me["f"]))
        fr = f"f ∈ [{me['f'][ok].min():.3f}, {me['f'][ok].max():.3f}] c/λ0 ({ok.sum()} freqs)"
        lim = TH[item if item in TH else item + "-1"]["abs"] if item in ("P2", "P3") else None
        if "R" in keys:
            t = lim or TH["P1-1"]["abs"]
            T.row(f"{item}-1" if item != "P1" else "P1-1", f"{label}: max |R−R_disc| over {fr}", "0", f"{dR:.2e}",
                  f"< {t:g}", dR < t, grid=gi)
        if "T" in keys:
            t = lim or TH["P1-2"]["abs"]
            T.row(f"{item}-1" if item != "P1" else "P1-2", f"{label}: max |T−T_disc|", "0", f"{dT:.2e}", f"< {t:g}",
                  dT < t, grid=gi)
        if "RT" in keys and item == "P1":
            t = TH["P1-3"]["abs"]
            T.row("P1-3", f"{label}: max |R+T−1|", "0", f"{cons:.2e}", f"< {t:g}", cons < t, grid=gi)
        return me, th, ok

    for pol in ("s", "p"):
        judge("P1", f"half-space ε=2.25, kx/k0={KN['kxn']}, kz/k0={KN['kzn']}, {pol}", R[f"p1_half_{pol}"][0], gi_u)
    for pol in ("s", "p"):
        judge("P2", f"normal incidence, half-space ε=2.25, {pol}", R[f"p2_normal_{pol}"][0], gi_u, keys=("R", "T"))
    for pol in ("s", "p"):
        judge("P3", f"film F1 (nonuniform y), kx/k0={KN['kxn']}, kz/k0={KN['kzn']}, {pol}", R[f"p3_film_{pol}"][0],
              nc.grid_info(geo, dtf) + ", Bloch cell 4x4", keys=("R", "T"))

    for pol in ("s", "p"):   # INFO: the superseded SF-flux definition (nu_failures.json, gate P3-1)
        o = R[f"p3_film_{pol}"][0]
        me, th = pa.measure(o), pa.theory(o)
        ok = band(me["f"])
        e = float(np.max(np.abs(me["R_sf"] - th["R"])[ok]))
        T.row("P3-info", f"film F1, {pol}: SF-flux R_sf (= |r − r_grade|², aux grading echo) vs R_disc", "INFO",
              f"{e:.2e}", "INFO", None, note="superseded definition, kept for the record")

    # ---------------------------------------------------------------- P4 cutoff leakage control
    o4, m4, err4 = R["p4_cutoff"]
    if not err4:   # cached run: stderr not kept, re-derive the warning condition the solver uses
        fcut = math.hypot(m4["kx"], m4["kz"]) / (2 * math.pi)
        warned = math.exp(-0.5 * ((m4["pulse"]["fc"] - fcut) / m4["pulse"]["fw"]) ** 2) > 1e-3
        src = "(cached run: condition re-evaluated)"
    else:
        warned = "WARNING: pulse amplitude at the cutoff" in err4
        src = ""
    me4, th4 = pa.measure(o4), pa.theory(o4)
    e4 = float(np.nanmax(np.abs(me4["R"] - th4["R"])))
    T.row("P4", "cutoff control kx/k0=0.5, kz/k0=1/3 (f_cut=0.601), fw=0.15: solver prints the cutoff warning", "yes",
          "yes" if warned else "no", "warning", warned, grid=gi_u, note=src)
    T.row("P4", "same: max |R−R_disc| over 0.7..1.3 (expected large: truncated cutoff ringing)", "INFO", f"{e4:.2e}",
          "INFO", None, grid=gi_u)

    figures(detail, W, nW)
    T.save(extra=dict(thresholds_status=TH["_status"], detail=detail), name="pulse_gate")
    return 0 if T.ok() else 1


def figures(detail, W, nW):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases = [k for k in detail if k.startswith(("p1_", "p2_", "p3_"))]
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for k in cases:
        d = detail[k]
        ax[0].plot(d["f"], d["R"], "o", ms=3, label=f"{k} FDTD")
        ax[0].plot(d["f"], d["R_th"], "-", lw=0.8, color="k")
        ax[1].semilogy(d["f"], np.abs(d["R"] - d["R_th"]) + 1e-16, label=f"{k} |ΔR|")
    ax[1].axhline(TH["P1-1"]["abs"], color="#b91c1c", ls="--", lw=0.8, label="threshold")
    ax[0].set_ylabel("R"); ax[0].legend(fontsize=7); ax[0].set_title("Broadband pulse, fixed (kx, kz): FDTD vs exact discrete TMM (lines)")
    ax[1].set_ylabel("|R − R_disc|"); ax[1].set_xlabel("f (c/λ0)"); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "P1_P3_R_spectra.png"), dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.semilogy(nW, np.abs(W / W[0] - 1) + 1e-18)
    ax.axhline(TH["P0-2"]["rel_drift"], color="#b91c1c", ls="--", lw=0.8, label="threshold")
    ax.set_xlabel("step"); ax.set_ylabel("|W(n)/W(0) − 1|"); ax.set_title("P0-2 Bloch energy (PEC y, random init)")
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIGS, "P0_2_energy.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
