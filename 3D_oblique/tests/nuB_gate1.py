"""Stage B gates B1-1, B1-2, B1-3, B2-1, B2-2 (SPEC_nonuniform §20.5, amendments D16, D20, D22).

Families B_{vac,film}_k{1,2,4} (x graded, s-pol, classical incidence m=1, n=0). All quantities from the weighted
Floquet projection (order 0) of every y row (dft_proj.bin) or from y-plane DFT slices.
B1-1 leakage of the analytic injection (i): backward amplitude in the SF region minus the far-PML echo (the
      backward amplitude in the TF region, which crosses the TF/SF plane unchanged), relative to the forward amplitude.
B1-2 R, T by (i) TF/SF-analytic and (ii) current sheet + vacuum normalization run (amendment D20): at every level
      |R_i − R_ii| ≤ 2|r||L_i| + |L_i|² (R_i carries the coherent sum of the film reflection and (i)'s own SF leakage
      L_i); order of |T_i − T_ii| ≥ 1.8; |R_i − R_ii| after subtracting (i)'s leakage: INFO.
B1-3 (amendment D22) leakage of the modal injection (iii) inc=m: x–y DFT slice, every row projected on each
      transverse mode profile; per mode a forward/backward fit (the mode's discrete ky) in the TF region and a
      backward fit in the SF region; leak = max_q |L_SF,q − b_TF,q| / |a_q|. INFO: inc=m vs inc=a on a uniform grid.
B2-1 x–z phase residual (max over TF planes of the RMS phase residual vs the Floquet mode).
B2-2 power in non-specular orders on a TF plane of the vacuum run (INFO + order).
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_runs as RN  # noqa: E402
import analyze_nu as an  # noqa: E402
import tmm  # noqa: E402

TH = nc.thresholds()["stageB"]
T = nc.Table("B_gate1")
KS = (1, 2, 4)
POL = "s"


def fit_back(u, y, ky):
    A = np.exp(-1j * ky * y)[:, None]
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    return complex(sol[0])


def two_wave(u, y, ky):
    A = np.stack([np.exp(1j * ky * y), np.exp(-1j * ky * y)], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    return complex(sol[0]), complex(sol[1])


def flux_rows(a, used, rows):
    yv, yd = used["y"], used["y_dual"]
    S = []
    for j in rows:
        wa = (yd[j] - yv[j]) / (yd[j] - yd[j - 1])
        wb = (yv[j] - yd[j - 1]) / (yd[j] - yd[j - 1])
        Hx = wa * a["Hx"][j - 1] + wb * a["Hx"][j]
        Hz = wa * a["Hz"][j - 1] + wb * a["Hz"][j]
        S.append(0.5 * np.real(a["Ez"][j] * np.conj(Hx) - a["Ex"][j] * np.conj(Hz)))
    return np.array(S)


def comp(a):
    return "Ez" if np.nanmax(np.abs(a["Ez"])) >= np.nanmax(np.abs(a["Ex"])) else "Ex"


def modal_leak(out, geo):
    """B1-3 (D22): per-mode leakage of the modal injection from the x–y DFT slice and modes.json."""
    used, meta = nc.load_used(out), nc.load_meta(out)
    F, _ = nc.load_dft(out, "xy", meta)
    with open(os.path.join(out, "modes.json")) as fh:
        M = json.load(fh)
    c = "Ez" if np.nanmax(np.abs(F["Ez"])) >= np.nanmax(np.abs(F["Ex"])) else "Ex"
    prim = c == "Ez"                                   # Ez: x primal (weights dx); Ex: x dual (weights hx)
    w = used["dx"] if prim else used["hx"]
    y = used["y"]
    sf, tf, _ = regions(geo)
    rows = []
    for P in M["product"]:
        pr = np.array(M["x"]["modes"][P["a"]]["primal" if prim else "dual"])
        prof = pr[:, 0] + 1j * pr[:, 1]
        u = (np.conj(prof) * w) @ F[c] / np.sum(w * np.abs(prof) ** 2)
        a, b = two_wave(u[tf], y[tf], P["ky"])
        L = fit_back(u[sf], y[sf], P["ky"])
        rows.append(dict(kappa=P["K"][0], leak=abs(L - b) / abs(a), echo=abs(b) / abs(a)))
    return rows


def uniform_m_vs_a():
    """D22 INFO: on an x/z-uniform grid the modal injection reduces to the aux line (one exact mode)."""
    g = "P_uniform_nl10"
    geo = nc.Geo(nc.load_grid(g))
    dt = geo.dt_run()
    Np = int(round(1.0 / dt))
    F = {}
    for inc in ("a", "m"):
        out = os.path.join(nc.RUNS, "stageB", f"{g}_p_{inc}_cmp")
        meta, _, _ = nc.run(out, grid=g, pol="p", inc=inc, m=1, n=1, dt=repr(dt), energy_every=0, nsteps=40 * Np,
                            dft0=30 * Np, dft1=40 * Np, zk=0)
        F[inc], _ = nc.load_dft(out, "xy", meta)
    scale = max(np.nanmax(np.abs(F["a"][c])) for c in nc.COMPS)
    return max(np.nanmax(np.abs(F["a"][c] - F["m"][c])) for c in nc.COMPS) / scale


def regions(geo):
    sf = np.arange(geo.npml_lo + 2, geo.j0 - 1)
    tfv = geo.marks["tf_vac"]
    tf = np.arange(geo.j0 + 2, tfv[0] + tfv[1] - 1)
    sub = geo.marks["sub"]
    subr = np.arange(sub[0] + 2, sub[0] + sub[1] - 1)
    return sf, tf, subr


def main():
    nc.build()
    res = {}
    for k in KS:
        gv, gf = f"B_vac_k{k}", f"B_film_k{k}"
        geo = nc.Geo(nc.load_grid(gv))
        mid = geo.j0 + geo.marks["tf_vac"][1] // 2
        o_vi, m_vi, _ = RN.stageB(gv, POL, "p", planar=1, yplanes=[mid])
        o_fi, m_fi, _ = RN.stageB(gf, POL, "p")
        o_vj, m_vj, _ = RN.stageB(gv, POL, "j")
        o_fj, m_fj, _ = RN.stageB(gf, POL, "j")
        used = nc.load_used(o_vi)
        a_vi, a_fi = RN.load_proj(o_vi, m_vi), RN.load_proj(o_fi, m_fi)
        a_vj, a_fj = RN.load_proj(o_vj, m_vj), RN.load_proj(o_fj, m_fj)
        sf, tf, sub = regions(geo)
        c = comp(a_vi)
        y = used["y"]
        h = float(geo.hy[geo.j0])
        ky, _ = tmm.ky_three_point(a_vi[c][tf], h)
        # B1-1 leakage
        L = fit_back(a_vi[c][sf], y[sf], ky)
        fa, fb = two_wave(a_vi[c][tf], y[tf], ky)
        leak = abs(L - fb) / abs(fa)
        # B1-2 R, T
        gfi = nc.Geo(nc.load_grid(gf))
        sff, tff, subf = regions(gfi)
        Lf = fit_back(a_fi[c][sff], y[sff], ky)
        R_i = abs(Lf / fa) ** 2
        R_ic = abs((Lf - L) / fa) ** 2                 # (i) with its own SF leakage removed
        S_inc_i = flux_rows(a_vi, used, tf).mean()
        T_i = flux_rows(a_fi, used, subf).mean() / S_inc_i
        cj = comp(a_vj)
        fa_j, _ = two_wave(a_vj[cj][tf], y[tf], ky)
        diff = {q: a_fj[q] - a_vj[q] for q in a_fj}
        film0 = gfi.marks["film"][0]
        rr = np.arange(geo.j0 + 2, film0 - 1)
        b_j = fit_back(diff[cj][rr], y[rr], ky)
        R_ii = abs(b_j / fa_j) ** 2
        S_inc_j = flux_rows(a_vj, used, tf).mean()
        T_ii = flux_rows(a_fj, used, subf).mean() / S_inc_j
        # B2-1 planarity (vacuum run, TF planes)
        pl = RN.load_planar(o_vi, m_vi)
        top = geo.Ny - geo.npml_hi - 1
        pres = max(float(np.nanmax(pl[q][geo.j0 + 1:top])) for q in ("Ez", "Hx", "Hy"))
        # B2-2 Floquet purity on the mid TF plane
        D, _ = nc.load_dft(o_vi, f"y{mid}", m_vi)
        kxs = np.array([m_vi["kx"] + 2 * math.pi * p / used["x"][-1] for p in range(-4, 5)])
        S, _ = an.order_fluxes(D, used, kxs, m_vi["kz"])
        prop = np.abs(kxs) < 2 * math.pi
        S0 = S[4]
        frac = float(np.sum(np.abs(S[prop])) - abs(S0)) / abs(S0)
        res[k] = dict(leak=leak, L_i=abs(L / fa), R_ic=R_ic, R_i=R_i, T_i=T_i, R_ii=R_ii, T_ii=T_ii, planar=pres, floquet=frac,
                      grid=nc.grid_info(geo, m_vi["dt"]))
        print(k, res[k], flush=True)
    hs = [1.0 / k for k in KS]
    w = TH["B1-1"]["order_window"]

    def orow(item, what, vals, window):
        o = nc.order_fit(hs, vals)
        T.row(item, what, "2", f"{o:.3f} (values {', '.join(f'{v:.2e}' for v in vals)})", f"∈ {window}",
              window[0] <= o <= window[1], grid="B_{vac,film}_k{1,2,4} (x graded; k = refinement)")
        return o
    orow("B1-1", "order of the analytic-injection leakage |L_SF − echo|/|a| (i)", [res[k]["leak"] for k in KS], w)
    dR = [abs(res[k]["R_i"] - res[k]["R_ii"]) for k in KS]
    dT = [abs(res[k]["T_i"] - res[k]["T_ii"]) for k in KS]
    th12 = TH["B1-2"]
    for k, d in zip(KS, dR):
        r = res[k]
        bound = 2 * math.sqrt(r["R_ii"]) * r["L_i"] + r["L_i"] ** 2
        T.row("B1-2", f"|R_(i) − R_(ii)| at k={k} vs the leakage bound 2|r||L_(i)| + |L_(i)|² (D20)", f"≤ {bound:.2e}",
              f"{d:.2e}", "≤ bound", d <= bound, grid=r["grid"], note=f"|L_(i)| = {r['L_i']:.2e}")
    oT = nc.order_fit(hs, dT)
    T.row("B1-2", "order of |T_(i) − T_(ii)| (D20: lower bound only)", "≥ 2",
          f"{oT:.3f} (values {', '.join(f'{v:.2e}' for v in dT)})", f"≥ {th12['T_order_min']}", oT >= th12["T_order_min"],
          grid="B_{vac,film}_k{1,2,4} (x graded; k = refinement)")
    T.row("B1-2", "|R_(i) − R_(ii)| with (i)'s own SF leakage subtracted", "→ 0",
          ", ".join(f"k={k}: {abs(res[k]['R_ic'] - res[k]['R_ii']):.1e}" for k in KS), "INFO", None,
          note="raw: " + ", ".join(f"{v:.2e}" for v in dR))
    th13 = TH["B1-3"]["max_rel"]
    for k in KS:
        gv = f"B_vac_k{k}"
        geo = nc.Geo(nc.load_grid(gv))
        o_m, m_m, _ = RN.stageB(gv, POL, "m", zk=0)
        rows = modal_leak(o_m, geo)
        lk = max(r["leak"] for r in rows)
        T.row("B1-3", f"modal injection (iii) inc=m: max over modes |L_SF − b_TF|/|a| at k={k} (D22)", "0",
              f"{lk:.2e}", f"< {th13:g}", lk < th13, grid=nc.grid_info(geo, m_m["dt"]),
              note="modes κx = " + ", ".join(f"{r['kappa']:.6f}" for r in rows)
                   + f"; far-PML echo {max(r['echo'] for r in rows):.1e} (crosses the TF/SF plane unchanged)")
    dma = uniform_m_vs_a()
    T.row("B1-3", "inc=m vs inc=a on an x/z-uniform grid (P_uniform_nl10, p-pol), max |ΔF|/max|F| (D22)", "0",
          f"{dma:.1e}", "INFO", None, grid="P_uniform_nl10")
    orow("B2-1", "order of the x–z phase residual (max over TF planes, vacuum, inc=p)", [res[k]["planar"] for k in KS],
         TH["B2-1"]["order_window"])
    fr = [res[k]["floquet"] for k in KS]
    T.row("B2-2", "power fraction in non-specular propagating orders (vacuum, TF plane)", "→ 0",
          ", ".join(f"k={k}: {res[k]['floquet']:.2e}" for k in KS), "INFO", None,
          note=f"order {nc.order_fit(hs, fr):.2f}" if min(fr) > 0 else "")
    Rt, Tt = nc.pred(f"B/film/{POL}/R_TMM"), nc.pred(f"B/film/{POL}/T_TMM")
    for k in KS:
        r = res[k]
        T.row("B1-2", f"R, T at k={k} (s): (i) / (ii) vs TMM", f"R {Rt:.6f}, T {Tt:.6f}",
              f"(i) R {r['R_i']:.6f} T {r['T_i']:.6f}; (ii) R {r['R_ii']:.6f} T {r['T_ii']:.6f}", "INFO", None,
              grid=r["grid"], note=f"R+T−1: (i) {r['R_i'] + r['T_i'] - 1:+.1e}, (ii) {r['R_ii'] + r['T_ii'] - 1:+.1e}")
    T.save(extra=dict(levels={str(k): v for k, v in res.items()}))
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
