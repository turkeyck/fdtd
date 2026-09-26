"""Stage B gate B3: lamellar grating (x nonuniform, nodes on the ridge edges), conical incidence (A2), vs RCWA.
SPEC_nonuniform §20.5 + amendment D15.

Per-order efficiencies from y-plane DFT slices (E at y_j, H at y_{j+1/2}: discrete-conserved pairing), Floquet
coefficients with the quadrature weights of the nonuniform x nodes:
  (i)  inc=p (TF/SF with the analytic continuous wave): R_p on an SF plane (scattered field), T_p on a substrate plane;
  (ii) inc=j (current sheet): R_p from (run − vacuum normalization run) between source and grating, T_p in the substrate.
Incident flux: the order-0 flux of the vacuum normalization run (same grid, same source) on a TF plane.
Judged: (i). (ii) and |(i) − (ii)| reported.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_runs as RN  # noqa: E402
import analyze_nu as an  # noqa: E402

TH = nc.thresholds()["stageB"]
T = nc.Table("B_gate3")
G, GV = "B_grating_ppw40", "B_grating_ppw40_vac"
FIG = os.path.join(nc.FIGS, "stageB")


def planes(geo):
    jR_i = geo.npml_lo + (geo.j0 - geo.npml_lo) // 2          # SF plane (inc=p)
    jR_j = geo.j0 + 20                                           # between the sheet (j0) and the grating (inc=j)
    s0, n = geo.marks["sub"]
    jT = s0 + n // 2
    return jR_i, jR_j, jT


def main():
    nc.build()
    os.makedirs(FIG, exist_ok=True)
    geo = nc.Geo(nc.load_grid(G))
    jR_i, jR_j, jT = planes(geo)
    yp = [jR_i, jR_j, jT]
    kxs_all = None
    out = {}
    for pol in ("s", "p"):
        rc_R = nc.pred(f"B/rcwa/{pol}/R_p")
        rc_T = nc.pred(f"B/rcwa/{pol}/T_p")
        oR = nc.pred(f"B/rcwa/{pol}/orders_R")
        oT = nc.pred(f"B/rcwa/{pol}/orders_T")
        eff = {}
        for inc in ("p", "j"):
            if inc == "p":     # D19: the grating rings (guided-mode resonance near f = 1.039): long run, long window
                o_g, m_g, _ = RN.stageB(G, pol, inc, tag="_long450", m=1, n=1, periods=450, dft_periods=150, yplanes=yp)
            else:
                o_g, m_g, _ = RN.stageB(G, pol, inc, m=1, n=1, periods=90, dft_periods=30, yplanes=yp)
            o_v, m_v, _ = RN.stageB(GV, pol, inc, m=1, n=1, periods=90, dft_periods=30, yplanes=yp)
            used = nc.load_used(o_g)
            Lx = used["x"][-1]
            kxs = np.array([m_g["kx"] + 2 * math.pi * p / Lx for p in range(-6, 7)])
            kxs_all = kxs
            idx = {p: q for q, p in enumerate(range(-6, 7))}
            # incident flux (order 0, vacuum normalization run, TF plane jR_j)
            Dv, _ = nc.load_dft(o_v, f"y{jR_j}", m_v)
            Sv, _ = an.order_fluxes(Dv, used, kxs, m_g["kz"])
            S_inc = Sv[idx[0]]
            if inc == "p":
                DR, _ = nc.load_dft(o_g, f"y{jR_i}", m_g)
            else:
                Dg, _ = nc.load_dft(o_g, f"y{jR_j}", m_g)
                DR = {c: Dg[c] - Dv[c] for c in Dg}
            SR, _ = an.order_fluxes(DR, used, kxs, m_g["kz"])
            DT, _ = nc.load_dft(o_g, f"y{jT}", m_g)
            ST, _ = an.order_fluxes(DT, used, kxs, m_g["kz"])
            R = {p: -SR[idx[p]] / S_inc for p in oR}
            Tt = {p: ST[idx[p]] / S_inc for p in oT}
            allR = -SR / S_inc
            allT = ST / S_inc
            eff[inc] = dict(R=R, T=Tt, sum=float(sum(R.values()) + sum(Tt.values())),
                            evan=float(np.sum(np.abs(allR)) + np.sum(np.abs(allT)) - sum(abs(v) for v in R.values())
                                       - sum(abs(v) for v in Tt.values())), runtime=m_g["runtime_s"])
        out[pol] = eff
        e = eff["p"]
        dmax = max([abs(e["R"][p] - rc_R[q]) for q, p in enumerate(oR)] + [abs(e["T"][p] - rc_T[q]) for q, p in enumerate(oT)])
        T.row("B3-1", f"per-order efficiency |η_FDTD − η_RCWA|, max over propagating orders ({pol}, injection (i))",
              "0", f"{dmax:.2e}", f"< {TH['B3-1']['abs']:g}", dmax < TH["B3-1"]["abs"], grid=nc.grid_info(geo, geo.dt_run()),
              note="R_p: " + ", ".join(f"{p}: {e['R'][p]:.5f}/{rc_R[q]:.5f}" for q, p in enumerate(oR))
                   + "; T_p: " + ", ".join(f"{p}: {e['T'][p]:.5f}/{rc_T[q]:.5f}" for q, p in enumerate(oT)))
        T.row("B3-2", f"Σ R_p + Σ T_p − 1 ({pol}, (i))", "0", f"{e['sum'] - 1:+.2e}", f"|·| < {TH['B3-2']['abs']:g}",
              abs(e["sum"] - 1) < TH["B3-2"]["abs"], grid=G)
        e2 = eff["j"]
        d2 = max([abs(e2["R"][p] - rc_R[q]) for q, p in enumerate(oR)] + [abs(e2["T"][p] - rc_T[q]) for q, p in enumerate(oT)])
        dij = max([abs(e2["R"][p] - e["R"][p]) for p in oR] + [abs(e2["T"][p] - e["T"][p]) for p in oT])
        T.row("B3-1", f"injection (ii) current sheet + normalization ({pol})", "RCWA", f"max |Δη| {d2:.2e}; Σ−1 "
              f"{e2['sum'] - 1:+.1e}; max |η_(i) − η_(ii)| {dij:.2e}", "INFO", None, grid=G)
    for pol in ("s", "p"):
        T.row("B3-3", f"RCWA order convergence 161 → 321 orders ({pol}) (D15)", "< 1e-5",
              f"{nc.pred(f'B/rcwa/{pol}/maxdiff_321'):.2e}", "< 1e-5", nc.pred(f"B/rcwa/{pol}/maxdiff_321") < 1e-5,
              grid="rcwa.py")
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, pol in zip(axes, ("s", "p")):
        oR, oT = nc.pred(f"B/rcwa/{pol}/orders_R"), nc.pred(f"B/rcwa/{pol}/orders_T")
        rc_R, rc_T = nc.pred(f"B/rcwa/{pol}/R_p"), nc.pred(f"B/rcwa/{pol}/T_p")
        xr = np.arange(len(oR))
        xt = np.arange(len(oT)) + len(oR) + 1
        ax.bar(xr - 0.2, rc_R, 0.4, color="#475569", label="RCWA (321 orders)")
        ax.bar(xr + 0.2, [out[pol]["p"]["R"][p] for p in oR], 0.4, color="#1d4ed8", label="FDTD (i)")
        ax.bar(xt - 0.2, rc_T, 0.4, color="#475569")
        ax.bar(xt + 0.2, [out[pol]["p"]["T"][p] for p in oT], 0.4, color="#1d4ed8")
        ax.set_xticks(list(xr) + list(xt))
        ax.set_xticklabels([f"R{p}" for p in oR] + [f"T{p}" for p in oT], fontsize=7)
        ax.set_title(f"lamellar grating, conical A2, {pol}-pol", fontsize=9)
        ax.set_ylabel("efficiency")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "grating_orders.png"), dpi=130)
    plt.close(fig)
    T.save(extra=dict(efficiencies={pol: {inc: {"R": {str(k): v for k, v in e["R"].items()},
                                                "T": {str(k): v for k, v in e["T"].items()}, "sum": e["sum"]}
                                          for inc, e in out[pol].items()} for pol in out}))
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
