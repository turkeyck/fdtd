"""Gate 1 (stage A, vacuum): SPEC_nonuniform §20.2, items 1-0 .. 1-7 (+ amendments D8, D9, D12).

Every measurement is made on the solver output; theory values come from results/nu_predictions.json (v2/*).
The per-row reduced phasors a_c(j) (dft_proj.bin) are the x-z Floquet projections of the main-grid DFT; they carry
no negative-frequency image because the DFT window is an integer number of periods (D12) and the projection is
orthogonal to exp(-i(kx x + kz z)).
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_runs as R  # noqa: E402
import tmm  # noqa: E402

TH = nc.thresholds()
G1 = TH["gate1"]
AM = TH["amendments"]
T = nc.Table("gate1")
FIG = os.path.join(nc.FIGS, "gate1")
POLS = ("s", "p")


# ----------------------------------------------------------------------------- helpers
def ecomp(a):
    """The tangential-E component with the larger amplitude (both carry the same reduced profile)."""
    return "Ex" if np.nanmax(np.abs(a["Ex"])) >= np.nanmax(np.abs(a["Ez"])) else "Ez"


def block_nodes(geo, tag, lo=None, hi=None):
    """Primal nodes whose 3-point stencil (cells j-1, j) lies inside block `tag`."""
    c0, n = geo.marks[tag]
    js = np.arange(c0 + 1, c0 + n)
    if lo is not None:
        js = js[js >= lo]
    if hi is not None:
        js = js[js <= hi]
    return js


def graded_nodes(geo):
    js = set()
    for tag in ("down", "up"):
        if tag in geo.marks and geo.marks[tag][1]:
            c0, n = geo.marks[tag]
            js |= set(range(c0, c0 + n + 1))
    return sorted(js)


def region_tf(geo):
    return geo.j0 + 2, geo.Ny - geo.npml_hi - 2


def measure_R_uniform_region(a, geo, used, tag, lo, hi):
    """Two-wave fit (forward + backward) in a uniform block; ky from the three-point estimator. Returns
    (R = |b/a|^2, ky, fit residual)."""
    c = ecomp(a)
    js = block_nodes(geo, tag, lo, hi)
    js = np.arange(js[0] - 1, js[-1] + 2)
    u = a[c][js]
    h = float(geo.hy[js[0]])
    ky, _ = tmm.ky_three_point(u, h)
    fa, fb, res = two_wave(u, used["y"][js], ky)
    return abs(fb / fa) ** 2, ky, res


def two_wave(u, y, ky):
    A = np.stack([np.exp(1j * ky * y), np.exp(-1j * ky * y)], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    return complex(sol[0]), complex(sol[1]), float(np.linalg.norm(A @ sol - u) / np.linalg.norm(u))


def theta_series(a, geo, used, meta):
    c = ecomp(a)
    lo, hi = region_tf(geo)
    js = np.arange(lo, hi + 1)
    dt = meta["dt"]
    Kt = math.hypot(meta["Ktilde"][0], meta["Ktilde"][2])
    yy, d, ky = tmm.local_ky_forward(a[c][js], used["y"][js], dt, Kt)
    idx = js[1:-1]
    th_m = np.arctan2(Kt, ky)
    th_p = np.array([math.atan2(Kt, tmm.ky_local(dd, 1.0, dt, Kt)) for dd in d])
    kt = math.hypot(meta["kx"], meta["kz"])
    th_c = math.atan2(kt, tmm.ky_continuous(kt))
    return idx, yy, th_m, th_p, th_c


def db(x):
    return 20 * math.log10(max(x, 1e-300))


def legacy_pml_db(pol):
    """Legacy uniform result L1-6, N_pml = 20, Δ = λ0/20 (results/level1.json)."""
    with open(os.path.join(nc.RESULTS, "level1.json"), encoding="utf-8") as fh:
        L = json.load(fh)
    for r in L["table"]:
        if r["item"] == "L1-6" and f"({pol}-pol, N_pml=20, Δ=λ0/20)" in r["quantity"]:
            return float(str(r["measured"]).split()[0])
    raise KeyError("legacy PML row not found")


# ----------------------------------------------------------------------------- items
def item_10_11_12(grid, pol, out, meta, geo, used, full):
    L = nc.load_log(out)
    H0 = float(np.linalg.norm(meta["H0"]))
    gi = nc.grid_info(geo, meta["dt"])
    dE = np.nanmax(L["devE_ref"])
    dH = np.nanmax(L["devH_ref"]) / H0
    T.row("1-0", f"main field vs Re[aux_ref e^(i(kx x+kz z))], all y, every 10 T0 ({grid}, {pol})", "0",
          f"E {dE:.1e}, H·η0 {dH:.1e}", f"< {G1['1-0']['max_rel']:g}", max(dE, dH) < G1["1-0"]["max_rel"], grid=gi)
    nwin = R.front_window_steps(meta, used)
    w = L["n"] <= nwin
    lE = np.nanmax(L["leakE_ref"][w])
    lH = np.nanmax(L["leakH_ref"][w]) / H0
    T.row("1-1", f"SF leakage = |F − aux_ref scattered field|/E0, n ≤ {nwin} ({grid}, {pol})", "0",
          f"E {lE:.1e}, H·η0 {lH:.1e}", f"< {G1['1-1']['max_rel']:g}", max(lE, lH) < G1["1-1"]["max_rel"], grid=gi,
          note="D8")
    raw = max(np.max(L["maxE_SF"][w]), np.max(L["maxH_SF"][w]) / H0)
    T.row("1-1", f"literal SF max / E0 (same window) ({grid}, {pol})", "—", f"{raw:.2e}", "INFO (D8)", None, grid=gi,
          note="includes the physical numerical reflection of the grid grading")
    if full:
        pl = R.load_planar(out, meta)
        a = R.load_proj(out, meta)
        lo, hi = geo.j0 + 1, geo.Ny - geo.npml_hi - 1
        amax = max(np.nanmax(np.abs(a[c][lo:hi + 1])) for c in nc.COMPS)
        worst, wc = 0.0, ""
        for c in nc.COMPS:
            if np.nanmax(np.abs(a[c][lo:hi + 1])) < 1e-3 * amax:
                continue                                    # identically zero component (e.g. Ey for s)
            top = hi if c in ("Ex", "Ez", "Hy") else min(hi, geo.Ny - 1)
            v = float(np.nanmax(pl[c][lo:top + 1]))
            if v > worst:
                worst, wc = v, c
        T.row("1-2", f"x–z planarity: max over TF planes of RMS phase residual vs kx·x+kz·z ({grid}, {pol})", "0",
              f"{worst:.2e} rad ({wc})", f"< {G1['1-2']['rms_rad']:g} rad", worst < G1["1-2"]["rms_rad"], grid=gi,
              note=f"{hi - lo + 1} planes, all non-zero components")


def item_13a(grid, pol, a, geo, used, meta):
    c = ecomp(a)
    gi = nc.grid_info(geo, meta["dt"])
    for tag, lo, hi in (("tf_a", geo.j0 + 2, None), ("hold", None, None), ("tf_b", None, geo.Ny - geo.npml_hi - 1)):
        if tag not in geo.marks or geo.marks[tag][1] == 0:
            continue
        js = block_nodes(geo, tag, lo, hi)
        u = a[c][np.arange(js[0] - 1, js[-1] + 2)]
        h = float(geo.hy[js[0]])
        ky, imr = tmm.ky_three_point(u, h)
        kp = nc.pred(f"L1/{grid}/ky/{tag}")
        rel = abs(ky - kp) / kp
        T.row("1-3a", f"ky in sub-region {tag} (h=λ0/{1 / h:.4g}) ({grid}, {pol})", f"{kp:.12f} /λ0", f"{ky:.12f} /λ0",
              f"rel < {G1['1-3a']['rel']:g}", rel < G1["1-3a"]["rel"], grid=gi,
              note=f"rel {rel:.1e}; three-point over {len(js)} nodes; imag residual {imr:.1e}")


def item_13b(grid, pol, a, geo, used, meta, store):
    idx, yy, tm, tp, tc = theta_series(a, geo, used, meta)
    gset = set(graded_nodes(geo))
    sel = np.array([j in gset for j in idx])
    ratio = float(np.max(np.abs(tm[sel] - tp[sel]) / np.abs(tp[sel] - tc)))
    store[(grid, pol)] = dict(y=yy, idx=idx, tm=tm, tp=tp, tc=tc, sel=sel, err=float(np.mean(np.abs(tm[sel] - tc))))
    T.row("1-3b", f"taper angle: max |θ_meas−θ_pred| / |θ_pred−θ_cont| over graded nodes ({grid}, {pol})", "0",
          f"{ratio:.2e}", f"≤ {G1['1-3b']['ratio']}", ratio <= G1["1-3b"]["ratio"], grid=nc.grid_info(geo, meta["dt"]),
          note=f"{int(sel.sum())} nodes; mean |θ_meas−θ_cont| = {store[(grid, pol)]['err']:.3e} rad")


# ----------------------------------------------------------------------------- figures
def fig_theta(store):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, names in ((axes[0], ["L1_taper_r1.1_base20", "L1_taper_r1.1_base40", "L1_taper_r1.1_base80"]),
                      (axes[1], [f"L1_taper_r1.1_base{b}_xz" for b in (10, 20, 40)])):
        for q, g in enumerate(names):
            if (g, "s") not in store:
                continue
            s = store[(g, "s")]
            col = ["#1d4ed8", "#15803d", "#7c3aed"][q]
            ax.plot(s["y"], np.degrees(s["tp"] - s["tc"]), color="#475569", ls="--", lw=1)
            ax.plot(s["y"][s["sel"]], np.degrees(s["tm"][s["sel"]] - s["tc"]), "o", ms=3, color=col,
                    label=g.replace("L1_taper_", ""))
            ax.plot(s["y"][~s["sel"]], np.degrees(s["tm"][~s["sel"]] - s["tc"]), ".", ms=1.5, color=col)
        ax.set_xlabel("y [λ0]")
        ax.set_ylabel("θ − θ_cont [deg]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_title("θ(y), x/z fixed λ0/20 (dashed: θ_pred(h_loc))")
    axes[1].set_title("θ(y), x/z refined with the base (D9)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "theta_y.png"), dpi=130)
    plt.close(fig)


def fig_snapshot(out, meta, geo, used):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    raw = np.fromfile(os.path.join(out, "snapshots.bin"), dtype=np.float64)
    S = raw.reshape(-1, meta["Nx"], meta["Ny"] + 1)[-1]           # Ez at (x_i, y_j), z = z_{zk+1/2}
    x = used["x"][:-1]
    y = used["y"]
    dt = meta["dt"]
    Kt = math.hypot(meta["Ktilde"][0], meta["Ktilde"][2])
    Phi = np.concatenate([[0.0], np.cumsum([tmm.ky_local(h, 1.0, dt, Kt) * h for h in used["hy"]])])
    X, Y = np.meshgrid(x, y, indexing="ij")
    ph = meta["kx"] * X + Phi[None, :]
    fig, ax = plt.subplots(figsize=(12, 3.2))
    j_lo, j_hi = meta["j0"] - 10, meta["Ny"] - meta["npml_hi"]
    sub = (slice(None), slice(j_lo, j_hi + 1))
    ax.pcolormesh(Y[sub], X[sub], S[sub], cmap="RdBu_r", shading="nearest", vmin=-1, vmax=1)
    ax.contour(Y[sub], X[sub], np.cos(ph[sub] - ph[0, meta["j0"]]), levels=[0.999], colors="#111827",
               linewidths=0.6)
    for tag in ("down", "hold", "up"):
        if tag in geo.marks and geo.marks[tag][1]:
            c0, n = geo.marks[tag]
            ax.axvspan(y[c0], y[c0 + n], color="#f59e0b", alpha=0.08)
    ax.set_aspect("equal")
    ax.set_xlabel("y [λ0]")
    ax.set_ylabel("x [λ0]")
    ax.set_title("Ez snapshot (x–y, real coordinates, equal aspect); black: predicted phase fronts; shaded: taper/hold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "snapshot_xy.png"), dpi=130)
    plt.close(fig)


def fig_lines(xs, ys, labels, xlabel, ylabel, title, path, logy=True, hline=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for x, y, lab in zip(xs, ys, labels):
        ax.plot(x, y, "o-", label=lab)
    if logy:
        ax.set_yscale("log")
    if hline is not None:
        ax.axhline(hline, color="#b91c1c", ls="--", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    nc.build()
    os.makedirs(FIG, exist_ok=True)
    store = {}
    Rm = {}
    # --- main runs: base-20 taper and abrupt (full diagnostics), s and p
    for grid in ("L1_taper_r1.1_base20", "L1_abrupt_r4_base20"):
        for pol in POLS:
            out, meta, geo = R.level1(grid, pol, full=True)
            used = nc.load_used(out)
            a = R.load_proj(out, meta)
            item_10_11_12(grid, pol, out, meta, geo, used, True)
            item_13a(grid, pol, a, geo, used, meta)
            if grid.startswith("L1_taper"):
                item_13b(grid, pol, a, geo, used, meta, store)
            Rm[(grid, pol)] = measure_R_uniform_region(a, geo, used, "tf_a", geo.j0 + 2, None)
            if grid == "L1_taper_r1.1_base20" and pol == "s":
                fig_snapshot(out, meta, geo, used)
            # 1-5 energy flux
            lo, hi = geo.j0 + 1, geo.Ny - geo.npml_hi - 1
            yv, yd = used["y"], used["y_dual"]
            S, Sh = [], []
            for j in range(lo, hi + 1):
                wa = (yd[j] - yv[j]) / (yd[j] - yd[j - 1])
                wb = (yv[j] - yd[j - 1]) / (yd[j] - yd[j - 1])
                Hx = wa * a["Hx"][j - 1] + wb * a["Hx"][j]
                Hz = wa * a["Hz"][j - 1] + wb * a["Hz"][j]
                S.append(0.5 * np.real(a["Ez"][j] * np.conj(Hx) - a["Ex"][j] * np.conj(Hz)))
                Hx2 = 0.5 * (a["Hx"][j - 1] + a["Hx"][j])
                Hz2 = 0.5 * (a["Hz"][j - 1] + a["Hz"][j])
                Sh.append(0.5 * np.real(a["Ez"][j] * np.conj(Hx2) - a["Ex"][j] * np.conj(Hz2)))
            S, Sh = np.array(S), np.array(Sh)
            var = float((S.max() - S.min()) / abs(S.mean()))
            T.row("1-5", f"time-averaged S_y over {len(S)} TF planes, distance-weighted H ({grid}, {pol})", "constant",
                  f"(max−min)/mean = {var:.2e}", f"< {G1['1-5']['rel_var']:g}", var < G1["1-5"]["rel_var"],
                  grid=nc.grid_info(geo, meta["dt"]))
            T.row("1-5", f"control: ½-average interpolation ({grid}, {pol})", "—",
                  f"{(Sh.max() - Sh.min()) / abs(Sh.mean()):.2e}", "INFO", None, grid=grid)
            if grid.startswith("L1_taper") and pol == "s":
                # 0-5 divergence (gate 0) on the same run: C log (nonuniform operator) + dumped fields, both operators
                import analyze_nu as an
                L = nc.load_log(out)
                F = nc.load_fields(out, f"fields_n{meta['nsteps']}.bin", meta)
                lo5, hi5 = geo.j0 + 2, geo.Ny - geo.npml_hi - 2
                K = meta["omega_tilde"]
                Emax = max(np.abs(F[c]).max() for c in ("Ex", "Ey", "Ez"))
                Hmax = max(np.abs(F[c]).max() for c in ("Hx", "Hy", "Hz"))
                dE = np.abs(an.div_E(F, used)[:, lo5:hi5 + 1, :]).max() / (K * Emax)
                dH = np.abs(an.div_H(F, used)[:, lo5:hi5 + 1, :]).max() / (K * Hmax)
                dEu = np.abs(an.div_E(F, used, uniform_delta=used["hx"][0])[:, lo5:hi5 + 1, :]).max() / (K * Emax)
                dHu = np.abs(an.div_H(F, used, uniform_delta=used["hx"][0])[:, lo5:hi5 + 1, :]).max() / (K * Hmax)
                logE = np.nanmax(L["divE_TF"]) / (K * Emax)
                G0 = nc.Table("gate0_divergence")
                th5 = TH["gate0"]["0-5"]
                G0.row("0-5", f"max|div E|, |div H| / (|K~| max|F|), TF j∈[{lo5},{hi5}] ({grid}, s, final step)", "0",
                       f"E {dE:.1e}, H {dH:.1e} (C log, all periods: E {logE:.1e})", f"< {th5['rel']:g}",
                       max(dE, dH, logE) < th5["rel"], grid=nc.grid_info(geo, meta["dt"]))
                G0.row("0-5", "control: uniform operator (every difference / Δx)", "clearly ≠ 0",
                       f"E {dEu:.2e}, H {dHu:.2e}", f"≥ {th5['control_min']:g}", min(dEu, dHu) >= th5["control_min"],
                       grid=grid, note="proves the test is sensitive to the metric")
                G0.save()
    # --- 1-3a / 1-3b / 1-4 on the refinement families
    for grid in ("L1_taper_r1.1_base40", "L1_taper_r1.1_base80"):
        for pol in POLS:
            out, meta, geo = R.level1(grid, pol)
            used = nc.load_used(out)
            a = R.load_proj(out, meta)
            item_10_11_12(grid, pol, out, meta, geo, used, False)
            item_13a(grid, pol, a, geo, used, meta)
            item_13b(grid, pol, a, geo, used, meta, store)
            Rm[(grid, pol)] = measure_R_uniform_region(a, geo, used, "tf_a", geo.j0 + 2, None)
    fam = [f"L1_taper_r1.1_base{b}_xz" for b in (10, 20, 40)]
    for grid in fam:
        for pol in POLS:
            out, meta, geo = R.level1(grid, pol)
            used = nc.load_used(out)
            a = R.load_proj(out, meta)
            item_10_11_12(grid, pol, out, meta, geo, used, False)
            item_13b(grid, pol, a, geo, used, meta, store)
    for pol in POLS:
        e = [store[(g, pol)]["err"] for g in fam]
        o = nc.order_fit([1 / 10, 1 / 20, 1 / 40], e)
        op = nc.pred(f"L1/theta_order_D9rev/{pol}")
        w = G1["1-3b"]["order_window"]
        T.row("1-3b", f"order of mean|θ_meas−θ_cont| under Δ-doubling, D9 family 10/20/40 ({pol})", f"2 (pred. {op:.3f})",
              f"{o:.3f} (errors {', '.join(f'{x:.2e}' for x in e)} rad)", f"∈ {w}", w[0] <= o <= w[1],
              grid="x/z and y refined together", note="amendment D9")
    fig_theta(store)
    # --- 1-4 reflection
    for grid in ("L1_abrupt_r4_base20",):
        for pol in POLS:
            Rmeas, ky, res = Rm[(grid, pol)]
            Rd = nc.pred(f"L1/{grid}/{pol}/R_disc")
            rel = abs(Rmeas - Rd) / Rd
            T.row("1-4a", f"abrupt r=4 reflectance vs exact discrete R_disc ({pol})", f"{Rd:.8e}", f"{Rmeas:.8e}",
                  f"rel < {G1['1-4a']['rel_vs_Rdisc']:g}", rel < G1["1-4a"]["rel_vs_Rdisc"], grid=grid,
                  note=f"rel {rel:.1e} ({db(math.sqrt(Rmeas)):.1f} dB amplitude); fit residual {res:.1e}; D1")
            est = nc.pred(f"L1/{grid}/{pol}/R_fresnel_est")
            T.row("1-4a", f"numerical-index Fresnel estimate of the fine slab ({pol})", f"{est:.4e}",
                  f"measured/estimate = {Rmeas / est:.2f} (power), {math.sqrt(Rmeas / est):.2f} (amplitude)", "INFO (D1)",
                  None, grid=grid)
    for pol in POLS:
        seq = []
        for b in (20, 40, 80):
            grid = f"L1_taper_r1.1_base{b}"
            Rmeas, ky, res = Rm[(grid, pol)]
            Rd = nc.pred(f"L1/{grid}/{pol}/R_disc")
            rel = abs(Rmeas - Rd) / Rd
            seq.append(Rmeas)
            T.row("1-4b", f"r=1.1 taper reflectance vs R_disc, base λ0/{b} ({pol})", f"{Rd:.8e}", f"{Rmeas:.8e}",
                  f"rel < {G1['1-4b']['rel_vs_Rdisc']:g}", rel < G1["1-4b"]["rel_vs_Rdisc"], grid=grid,
                  note=f"rel {rel:.1e}; {db(math.sqrt(Rmeas)):.1f} dB (amplitude; −60 dB INFO: "
                       f"{'below' if db(math.sqrt(Rmeas)) < -60 else 'above'})")
        dec = seq[0] > seq[1] > seq[2]
        T.row("1-4b", f"|R| strictly decreasing with refinement 20→40→80 ({pol})", "decreasing",
              " → ".join(f"{db(math.sqrt(x)):.1f} dB" for x in seq), "strict decrease", dec, grid="D2")
    # 1-4c scan (INFO): FDTD at base 20/40 (s), exact discrete theory for every (r, base)
    scan = {}
    for r in (1.05, 1.1, 1.2, 1.5, 4.0):
        for b in (20, 40, 80):
            if r == 1.1:
                g = f"L1_taper_r1.1_base{b}"
            elif r == 4.0:
                g = "L1_abrupt_r4_base20" if b == 20 else f"L1_scan_r4_base{b}"
            else:
                g = f"L1_scan_r{r:g}_base{b}"
            th = nc.pred(f"L1/{g}/s/R_disc")
            meas = Rm.get((g, "s"), (None,))[0]
            if meas is None and b in (20, 40):
                out, meta, geo = R.level1(g, "s")
                used = nc.load_used(out)
                meas = measure_R_uniform_region(R.load_proj(out, meta), geo, used, "tf_a", geo.j0 + 2, None)[0]
            scan[(r, b)] = (th, meas)
            T.row("1-4c", f"reflectance vs r_max and base: r={r:g}, base λ0/{b} (s)", f"{db(math.sqrt(th)):.2f} dB",
                  "—" if meas is None else f"{db(math.sqrt(meas)):.2f} dB", "INFO", None, grid=g,
                  note="FDTD not run at base 80 for the scan (theory only)" if meas is None else "")
    fig_lines([[20, 40, 80]] * 5, [[math.sqrt(scan[(r, b)][0]) for b in (20, 40, 80)] for r in (1.05, 1.1, 1.2, 1.5, 4.0)],
              [f"r={r:g} (R_disc)" for r in (1.05, 1.1, 1.2, 1.5, 4.0)], "base cells per λ0", "|r| (amplitude)",
              "1-4c: grading reflection vs r_max and resolution (exact discrete; FDTD agrees, table)",
              os.path.join(FIG, "reflection_scan.png"), hline=1e-3)
    # --- 1-1 literal on the uniform control grid
    for pol in POLS:
        out, meta, geo = R.level1("L1_uniform_ctrl", pol)
        used = nc.load_used(out)
        L = nc.load_log(out)
        nwin = R.front_window_steps(meta, used)
        w = L["n"] <= nwin
        H0 = float(np.linalg.norm(meta["H0"]))
        raw = max(np.max(L["maxE_SF"][w]), np.max(L["maxH_SF"][w]) / H0)
        T.row("1-1", f"uniform control grid: literal SF max / E0, n ≤ {nwin} ({pol})", "0", f"{raw:.2e}",
              f"< {G1['1-1']['uniform_ctrl_literal_max']:g}", raw < G1["1-1"]["uniform_ctrl_literal_max"],
              grid=nc.grid_info(geo, meta["dt"]))
    # --- 1-6 PML
    pml = {}
    for N in G1["1-6"]["npml_list"]:
        for pol in POLS if N == 20 else ("s",):
            grid = "L1_taper_r1.1_base20" if N == 60 else f"L1_taper_r1.1_base20_N{N}"
            out, meta, geo = R.level1(grid, pol, full=(N == 60))
            used = nc.load_used(out)
            a = R.load_proj(out, meta)
            Rp, ky, res = measure_R_uniform_region(a, geo, used, "tf_b", None, geo.Ny - geo.npml_hi - 1)
            pml[(N, pol)] = db(math.sqrt(Rp))
    for pol in POLS:
        leg = legacy_pml_db(pol)
        v = pml[(20, pol)]
        ok = abs(v - leg) <= G1["1-6"]["max_db_diff_vs_uniform"] and v < G1["1-6"]["max_db"]
        T.row("1-6", f"far-PML reflection, N_pml=20 ({pol})", f"legacy uniform {leg:.1f} dB", f"{v:.1f} dB",
              "|Δ| ≤ 3 dB and < −40 dB", ok, grid="L1_taper_r1.1_base20_N20 (PML at λ0/20)")
    T.row("1-6", "far-PML reflection vs N_pml (s)", "—",
          ", ".join(f"N={N}: {pml[(N, 's')]:.1f} dB" for N in G1["1-6"]["npml_list"]), "INFO", None, grid="base 20")
    fig_lines([G1["1-6"]["npml_list"]], [[10 ** (pml[(N, "s")] / 20) for N in G1["1-6"]["npml_list"]]], ["s"],
              "N_pml", "|r_PML| (amplitude)", "1-6: far-PML reflection vs thickness (nonuniform run, PML at λ0/20)",
              os.path.join(FIG, "pml.png"))
    # --- 1-7 stability
    out, meta, geo = R.level1("L1_taper_r1.1_base20", "s", tag="_stab", stability=True)
    L = nc.load_log(out)
    ok = ~np.isnan(L["W_total"])
    n, W = L["n"][ok], L["W_total"][ok]
    t_off = (40.0 + 20.0 + 3 * 4.0) / meta["dt"]           # turn-off ramp complete (erf centre + 3 tau)
    after = n >= t_off
    Wa, na = W[after], n[after]
    k0 = int(np.argmax(Wa))
    tail = Wa[k0:]
    end = np.where(tail < 1e-8 * tail[0])[0]
    kend = end[0] if len(end) else len(tail) - 1
    mono = bool(np.all(np.diff(tail[: kend + 1]) < 0)) and len(end) > 0
    half = Wa[len(Wa) // 2:]
    nogrow = bool(half.max() <= half[0])
    T.row("1-7", f"{meta['nsteps']} steps; after switch-off: monotone decay from peak to < 1e-8 peak, no late growth",
          "decay", f"{kend} samples monotone to {tail[kend] / tail[0]:.1e}; 2nd-half max/start {half.max() / half[0]:.6f}",
          "monotone & no growth", mono and nogrow, grid=nc.grid_info(geo, meta["dt"]))
    fig_lines([na], [Wa / Wa.max()], ["W_total (legacy-form energy)"], "step", "W / W_peak",
              "1-7: energy after switch-off (base-20 taper, s)", os.path.join(FIG, "energy_decay.png"))
    T.save()
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
