"""Gate 2 (stage A): thin film F1 and 5-layer stack F5 vs TMM (SPEC_nonuniform §20.3 + amendments D10, D11, D12).

R: two-wave (forward/backward) fit of the reduced tangential E in the uniform vacuum region after the TF/SF plane.
T: conserved, distance-weighted S_y in the uniform substrate region / exact incident flux of the injected discrete
   plane wave  S_inc = 1/2 cos(ky Δa/2) (E0z H0x − E0x H0z)  (Δa, ky, E0, H0 from meta.json).
R and T are measured independently, so R + T − 1 is a real energy test.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_runs as RN  # noqa: E402
import grid_gen as gg  # noqa: E402
import tmm  # noqa: E402

TH = nc.thresholds()
G2 = TH["gate2"]
T = nc.Table("gate2")
FIG = os.path.join(nc.FIGS, "gate2")
POLS = ("s", "p")
STRUCTS = ("F1", "F5")


def uniform_run_from(geo, j_start, direction=+1, stop=None):
    """Cells of constant spacing and eps starting at cell j_start going up (+1) or down (-1)."""
    h0, e0 = geo.hy[j_start], geo.eps_cell[j_start]
    j = j_start
    while 0 <= j + direction < geo.Ny and geo.hy[j + direction] == h0 and geo.eps_cell[j + direction] == e0:
        if stop is not None and j + direction == stop:
            break
        j += direction
    return (j_start, j) if direction > 0 else (j, j_start)


def measure(out, meta, geo):
    used = nc.load_used(out)
    a = RN.load_proj(out, meta)
    c = "Ex" if np.nanmax(np.abs(a["Ex"])) >= np.nanmax(np.abs(a["Ez"])) else "Ez"
    # vacuum region: cells [j0, c1] uniform -> nodes j0+2 .. c1 (stencils inside)
    c0, c1 = uniform_run_from(geo, geo.j0, +1)
    js = np.arange(geo.j0 + 1, c1 + 2)
    h = float(geo.hy[geo.j0])
    u = a[c][js]
    ky, _ = tmm.ky_three_point(u[1:], h)
    A = np.stack([np.exp(1j * ky * used["y"][js]), np.exp(-1j * ky * used["y"][js])], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    fit_res = float(np.linalg.norm(A @ sol - u) / np.linalg.norm(u))
    R = abs(sol[1] / sol[0]) ** 2
    # substrate region: uniform cells ending before the far PML
    top = geo.Ny - geo.npml_hi - 1
    s0, s1 = uniform_run_from(geo, top, -1)
    yv, yd = used["y"], used["y_dual"]

    def S_at(j):
        wa = (yd[j] - yv[j]) / (yd[j] - yd[j - 1])
        wb = (yv[j] - yd[j - 1]) / (yd[j] - yd[j - 1])
        Hx = wa * a["Hx"][j - 1] + wb * a["Hx"][j]
        Hz = wa * a["Hz"][j - 1] + wb * a["Hz"][j]
        return 0.5 * float(np.real(a["Ez"][j] * np.conj(Hx) - a["Ex"][j] * np.conj(Hz)))
    Ssub = np.array([S_at(j) for j in range(s0 + 1, s1 + 1)])
    Svac = np.array([S_at(j) for j in range(geo.j0 + 2, c1 + 1)])
    E0, H0 = np.array(meta["E0"]), np.array(meta["H0"])
    Sinc = 0.5 * math.cos(meta["ky"] * meta["Delta_a"] / 2) * (E0[2] * H0[0] - E0[0] * H0[2])
    Tm = float(Ssub.mean() / Sinc)
    L = nc.load_log(out)
    dev = max(np.nanmax(L["devE_ref"]), np.nanmax(L["devH_ref"]) / np.linalg.norm(H0))
    return dict(R=R, T=Tm, fit_res=fit_res, S_var=float(np.ptp(np.concatenate([Svac, Ssub])) / abs(Sinc)),
                R_flux=float(1 - Svac.mean() / Sinc), dev_ref=float(dev), runtime=meta["runtime_s"],
                cells=meta["Nx"] * (meta["Ny"] + 1) * meta["Nz"], steps=meta["nsteps"])


def run_measure(grid, pol, **kw):
    out, meta, geo = RN.level2(grid, pol, **kw)
    return measure(out, meta, geo), geo, meta


def uniform_aligned_grid(M, struct="F1"):
    """2-6: uniform grid with spacing D_FILM/M everywhere (all interfaces on nodes)."""
    h = gg.D_FILM / M
    layers = [(gg.N_FILM, gg.D_FILM)] if struct == "F1" else gg.STACK5
    blocks = [dict(kind="uniform", n=20, h=h, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=int(round(1.0 / h)), h=h, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(3.0 / h)), h=h, eps=1.0, tag="tf_vac")]
    for q, (n, d) in enumerate(layers):
        blocks.append(dict(kind="uniform", n=int(round(d / h)), h=h, eps=n * n, tag=f"layer{q}"))
    blocks += [dict(kind="uniform", n=int(round(3.0 / h)), h=h, eps=gg.N_SUB ** 2, tag="sub"),
               dict(kind="uniform", n=20, h=h, eps=gg.N_SUB ** 2, tag="pml_hi")]
    g = gg.make_grid(f"L2_{struct}_uniform_M{M}", blocks, gg.AX_X(), gg.AX_Z(), 1.1, 20, 20,
                     20 + int(round(1.0 / h)), ppw=None, extra=dict(M=M))
    gg.write(g)
    return g["name"]


def figure(conv, ctrl, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    cols = {"F1s": "#1d4ed8", "F1p": "#60a5fa", "F5s": "#15803d", "F5p": "#4ade80"}
    for key, (hs, err) in conv.items():
        ax.loglog(hs, np.abs(err), "o-", color=cols[key], label=f"aligned nonuniform {key} (D10)")
    for pol, (hs, err) in ctrl.items():
        ax.loglog(hs, np.abs(err), "s--", color="#b91c1c" if pol == "s" else "#f87171",
                  label=f"unaligned uniform control F1{pol} (D11)")
    x = np.array([0.003, 0.03])
    ax.loglog(x, 2e-3 * (x / 0.03) ** 2, "k:", lw=0.8, label="slope 2")
    ax.loglog(x, 2e-2 * (x / 0.03), "k-.", lw=0.8, label="slope 1")
    ax.set_xlabel("Δy_min [λ0]")
    ax.set_ylabel("|R − R_∞(Δx)|")
    ax.set_title("Gate 2: convergence of R with Δy (x/z fixed at λ0/20)", fontsize=9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    nc.build()
    os.makedirs(FIG, exist_ok=True)
    res = {}
    # ---- 2-1, 2-2, 2-7 default grid (D6)
    for st in STRUCTS:
        for pol in POLS:
            grid = f"L2_{st}_ppw40"
            m, geo, meta = run_measure(grid, pol)
            res[(grid, pol)] = m
            Rt, Tt = nc.pred(f"L2/{st}/{pol}/R_TMM"), nc.pred(f"L2/{st}/{pol}/T_TMM")
            gi = nc.grid_info(geo, meta["dt"])
            T.row("2-1", f"|R − R_TMM| ({st}, {pol}, default grid D6)", f"{Rt:.8f}", f"{m['R']:.8f}",
                  f"< {G2['2-1']['abs']:g}", abs(m["R"] - Rt) < G2["2-1"]["abs"], grid=gi,
                  note=f"Δ = {m['R'] - Rt:+.2e}")
            T.row("2-1", f"|T − T_TMM| ({st}, {pol}, default grid D6)", f"{Tt:.8f}", f"{m['T']:.8f}",
                  f"< {G2['2-1']['abs']:g}", abs(m["T"] - Tt) < G2["2-1"]["abs"], grid=gi,
                  note=f"Δ = {m['T'] - Tt:+.2e}")
            T.row("2-2", f"|R + T − 1| ({st}, {pol})", "0", f"{m['R'] + m['T'] - 1:+.2e}", f"< {G2['2-2']['abs']:g}",
                  abs(m["R"] + m["T"] - 1) < G2["2-2"]["abs"], grid=gi,
                  note=f"R by 2-wave fit, T by substrate flux; S_y spread {m['S_var']:.1e}")
            Rd = nc.pred(f"L2/{grid}/{pol}/R_disc")
            T.row("2-7", f"FDTD vs exact discrete R_disc, same grid ({st}, {pol})", f"{Rd:.10f}", f"{m['R']:.10f}",
                  "INFO", None, grid=grid, note=f"rel {abs(m['R'] - Rd) / Rd:.1e}; main vs aux_ref {m['dev_ref']:.1e}")
    # ---- 2-3 D10 bisection family (k = 1, 2, 4), 3-level order
    conv = {}
    rinf = {}
    for st in STRUCTS:
        for pol in POLS:
            Rs = []
            for k in (1, 2, 4):
                m, geo, meta = run_measure(f"L2_{st}_bis{k}", pol)
                res[(f"L2_{st}_bis{k}", pol)] = m
                Rs.append(m["R"])
            p3 = math.log2((Rs[0] - Rs[1]) / (Rs[1] - Rs[2]))
            Ri = Rs[2] + (Rs[2] - Rs[1]) / (2 ** p3 - 1)
            rinf[(st, pol)] = Ri
            hs = [0.02, 0.01, 0.005]
            conv[st + pol] = (hs, np.array(Rs) - Ri)
            pp = nc.pred(f"L2/order_y3/{st}/{pol}")
            w = G2["2-3"]["order_window"]
            T.row("2-3", f"order of R vs Δy, aligned nonuniform (D10: bis1,2,4; {st}, {pol})", f"2 (pred. {pp:.3f})",
                  f"{p3:.3f}", f"∈ {w}", w[0] <= p3 <= w[1], grid="L2_*_bis{1,2,4}, dt ∝ h, x/z λ0/20",
                  note="R = " + ", ".join(f"{r:.8f}" for r in Rs) + f"; R_∞(Δx) = {Ri:.8f}")
    # ---- 2-4 D11 control
    ctrl = {}
    for pol in POLS:
        Cs, hs = [], []
        for mh in (4.5, 9.5, 18.5):
            m, geo, meta = run_measure(f"L2_F1_ctrl_mh{mh:g}", pol)
            res[(f"L2_F1_ctrl_mh{mh:g}", pol)] = m
            Cs.append(m["R"])
            hs.append(gg.D_FILM / mh)
        o = nc.order_fit(hs, np.array(Cs) - rinf[("F1", pol)])
        ctrl[pol] = (hs, np.array(Cs) - rinf[("F1", pol)])
        w = G2["2-4"]["order_window"]
        T.row("2-4", f"control: unaligned uniform grid, order of R (D11; F1, {pol})",
              f"1 (pred. {nc.pred(f'L2/order_control3/{pol}'):.3f})", f"{o:.3f}", f"∈ {w}", w[0] <= o <= w[1],
              grid="L2_F1_ctrl_mh{4.5,9.5,18.5}", note="errors " + ", ".join(f"{e:+.2e}" for e in ctrl[pol][1]))
    figure(conv, ctrl, os.path.join(FIG, "convergence.png"))
    # ---- 2-5 error budget (INFO): y part, transverse part, x/z-only refinement
    for st in STRUCTS:
        for pol in POLS:
            Rt = nc.pred(f"L2/{st}/{pol}/R_TMM")
            Rdef = res[(f"L2_{st}_ppw40", pol)]["R"]
            ey = Rdef - rinf[(st, pol)]
            ex = rinf[(st, pol)] - Rt
            note = ""
            if st == "F1":
                m40, _, _ = run_measure("L2_F1_ppw40_dxz40", pol)
                d40 = m40["R"] - Rdef
                pr = nc.pred(f"L2/L2_F1_ppw40_dxz40/{pol}/R_disc") - nc.pred(f"L2/L2_F1_ppw40/{pol}/R_disc")
                note = f"x/z λ0/20→λ0/40 at fixed y grid: ΔR = {d40:+.2e} (theory {pr:+.2e})"
            T.row("2-5", f"error budget ({st}, {pol}): total = y-discretization + transverse (K~t vs k_t)",
                  f"total {Rdef - Rt:+.2e}", f"y {ey:+.2e}, transverse {ex:+.2e}", "INFO", None,
                  grid="default grid / D10 limit", note=note)
    # ---- 2-6 efficiency (INFO): uniform aligned grid reaching the same |R − R_TMM| as the default grid
    for pol in POLS:
        Rt = nc.pred(f"L2/F1/{pol}/R_TMM")
        target = abs(res[("L2_F1_ppw40", pol)]["R"] - Rt)
        chosen = None
        for M in range(7, 80):
            name = uniform_aligned_grid(M)
            geo = nc.Geo(nc.load_grid(name))
            Rd = tmm.tmm_discrete(geo.hy, geo.eps_cell, geo.dt_run(), nc.Kt_of(geo), pol, geo.j0 + 2,
                                  geo.Ny - geo.npml_hi - 2)["R"]
            if abs(Rd - Rt) <= target:
                chosen = name
                break
        if chosen is None:
            T.row("2-6", f"efficiency ({pol})", "—", "no uniform grid up to M=79 reaches the target", "INFO", None)
            continue
        m, geo, meta = run_measure(chosen, pol)
        d = res[("L2_F1_ppw40", pol)]
        T.row("2-6", f"efficiency: same |R−R_TMM| ({pol}): nonuniform default vs uniform aligned {chosen}",
              f"target {target:.2e}",
              f"nonuniform {d['cells'] / 1e6:.2f}M cells × {d['steps']} steps, {d['runtime']:.0f} s; uniform "
              f"{m['cells'] / 1e6:.2f}M cells × {m['steps']} steps, {m['runtime']:.0f} s (err {abs(m['R'] - Rt):.2e})",
              "INFO", None, grid="F1")
    T.save(extra=dict(measurements={f"{k[0]}|{k[1]}": v for k, v in res.items()}))
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
