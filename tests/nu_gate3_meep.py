"""Gate 3: comparison with Meep (SPEC_nonuniform §20.4 + amendments D13, D14). Exit 3 (BLOCKED) without Meep.

(A) converged values: thin Bloch cell, resolutions 40/80/160/320, Richardson; vs TMM, vs our D14 limit, vs our
    default grid (gate-2 measurement). 3A-1: thin vs full cell at resolution 40 (D13).
(B) transformation optics on the mapped grid: tangential DFT profiles (Bloch-projected, normalized by the forward
    amplitude of the incident wave) vs our run on L3B_mapped_{vac,film}; continuous metric (judged) and the
    discrete metric of our grid (3B-3 diagnostic).
"""
import json
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_runs as RN  # noqa: E402
import grid_gen as gg  # noqa: E402
import tmm  # noqa: E402

TH = nc.thresholds()
G3 = TH["gate3"]
T = nc.Table("gate3")
MEEP_PY = os.environ.get("MEEP_PYTHON", os.path.expanduser("~/micromamba/envs/mp/bin/python"))
MROOT = os.path.join(nc.RUNS, "meep")
FIG = os.path.join(nc.FIGS, "gate3")


def meep(script, out, args, result):
    """Run a Meep script (cached on its argument list)."""
    os.makedirs(out, exist_ok=True)
    argfile = os.path.join(out, "args.json")
    res = os.path.join(out, result)
    if not nc.FRESH and os.path.exists(res) and os.path.exists(argfile):
        with open(argfile) as fh:
            if json.load(fh) == args:
                return
    cmd = [MEEP_PY, os.path.join(nc.ROOT, script)] + args + ["--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{script} failed ({r.returncode}):\n{r.stderr[-3000:]}")
    with open(argfile, "w") as fh:
        json.dump(args, fh)


def meep_A(struct, pol, res, cell="thin"):
    norm = os.path.join(MROOT, f"A_norm_{pol}_r{res}_{cell}")
    meep("meep_ref_uniform.py", norm, ["--res", str(res), "--pol", pol, "--struct", "vac", "--cell", cell],
         "meep_A.json")
    out = os.path.join(MROOT, f"A_{struct}_{pol}_r{res}_{cell}")
    meep("meep_ref_uniform.py", out, ["--res", str(res), "--pol", pol, "--struct", struct, "--cell", cell,
                                      "--norm", norm], "meep_A.json")
    with open(os.path.join(out, "meep_A.json")) as fh:
        return json.load(fh)


# ----------------------------------------------------------------------------- (A)
def part_A():
    g2 = {}
    p2 = os.path.join(nc.RESULTS, "nu_gate2.json")
    if os.path.exists(p2):
        with open(p2, encoding="utf-8") as fh:
            g2 = json.load(fh).get("measurements", {})
    conv = {}
    for st in ("F1", "F5"):
        for pol in ("s", "p"):
            ress = G3["3A-2"].get("resolutions_D18", G3["3A-2"]["resolutions"])      # D18
            Rs = [meep_A(st, pol, r)["R"] for r in ress]
            Ts = [meep_A(st, pol, r)["T"] for r in ress]
            p = math.log2(abs((Rs[-3] - Rs[-2]) / (Rs[-2] - Rs[-1])))
            Rinf = Rs[-1] + (Rs[-1] - Rs[-2]) / (2 ** p - 1)               # 3-level Richardson (observed order)
            Tinf = Ts[-1] + (Ts[-1] - Ts[-2]) / (2 ** p - 1)
            Rt, Tt = nc.pred(f"L2/{st}/{pol}/R_TMM"), nc.pred(f"L2/{st}/{pol}/T_TMM")
            conv[(st, pol)] = (ress, Rs, Rinf, Rt)
            T.row("3A-2", f"Meep Richardson limit vs TMM ({st}, {pol}), res {ress}", f"R {Rt:.8f}, T {Tt:.8f}",
                  f"R {Rinf:.8f}, T {Tinf:.8f}", f"< {G3['3A-2']['abs']:g}",
                  max(abs(Rinf - Rt), abs(Tinf - Tt)) < G3["3A-2"]["abs"], grid="Meep thin Bloch cell",
                  note=f"ΔR {Rinf - Rt:+.1e}, ΔT {Tinf - Tt:+.1e}; observed order {p:.2f}; R(res) = "
                       + ", ".join(f"{r:.7f}" for r in Rs))
            # 3A-3: our D14 two-level limit
            own = []
            for k in (1, 2):
                out, meta, geo = RN.level2(f"L2_{st}_bisxz{k}", pol)
                import nu_gate2_film as g2m
                own.append(g2m.measure(out, meta, geo)["R"])
            Rown = own[1] + (own[1] - own[0]) / 3.0
            T.row("3A-3", f"own limit (D14: y + x/z halved, 2-level) vs Meep limit ({st}, {pol})", f"{Rinf:.8f}",
                  f"{Rown:.8f}", f"< {G3['3A-3']['abs']:g}", abs(Rown - Rinf) < G3["3A-3"]["abs"],
                  grid="L2_*_bisxz{1,2}", note=f"own R(k=1,2) = {own[0]:.8f}, {own[1]:.8f}; own limit − TMM "
                                              f"{Rown - Rt:+.1e} (pred. {nc.pred(f'L2/Rinf_xz2/{st}/{pol}') - Rt:+.1e})")
            key = f"L2_{st}_ppw40|{pol}"
            if key in g2:
                Rd = g2[key]["R"]
                T.row("3A-4", f"own default grid (D6) vs Meep limit ({st}, {pol})", f"{Rinf:.8f}", f"{Rd:.8f}",
                      f"< {G3['3A-4']['abs']:g}", abs(Rd - Rinf) < G3["3A-4"]["abs"], grid="L2_*_ppw40")
            else:
                T.row("3A-4", f"own default grid vs Meep limit ({st}, {pol})", "—", "gate-2 result missing",
                      "< 1e-3", False)
    # 3A-1 thin vs full (D13)
    r = G3["3A-1"]["resolution"] if "D13" not in TH["amendments"] else 40
    thin = meep_A("F1", "s", r, "thin")
    full = meep_A("F1", "s", r, "full")
    T.row("3A-1", f"Meep thin Bloch cell vs full Lx×Lz cell, res {r} (F1, s)", f"R {full['R']:.10f}",
          f"R {thin['R']:.10f}", f"< {G3['3A-1']['abs']:g}", abs(thin["R"] - full["R"]) < G3["3A-1"]["abs"],
          grid="Meep", note="D13 (resolution 40)")
    return conv


# ----------------------------------------------------------------------------- (B)
def reduced_profiles_own(out, meta):
    a = RN.load_proj(out, meta)
    return {c: a[c] for c in ("Ex", "Ez", "Hx", "Hz")}


def normalise(prof, geo, js_fit, ky, conj_ok=True):
    c = "Ex" if np.nanmax(np.abs(prof["Ex"])) >= np.nanmax(np.abs(prof["Ez"])) else "Ez"
    y = geo.y[js_fit]
    u = prof[c][js_fit]
    A = np.stack([np.exp(1j * ky * y), np.exp(-1j * ky * y)], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    conj = False
    if conj_ok and abs(sol[1]) > abs(sol[0]):          # opposite phasor convention: conjugate
        prof = {k: np.conj(v) for k, v in prof.items()}
        return normalise(prof, geo, js_fit, ky, conj_ok=False)
    return {k: v / sol[0] for k, v in prof.items()}, conj


def part_B():
    rows = {}
    for case in ("vac", "film"):
        grid = f"L3B_mapped_{case}"
        geo = nc.Geo(nc.load_grid(grid))
        pars = geo.g["generator"]
        u1 = pars["u14"][0]
        j_fit = np.arange(geo.j0 + 2, int(round(u1 / pars["du"])) - 1)
        for pol in ("s", "p"):
            out, meta, _ = RN.level2(grid, pol)
            own = reduced_profiles_own(out, meta)
            cmax = "Ex" if np.nanmax(np.abs(own["Ex"])) >= np.nanmax(np.abs(own["Ez"])) else "Ez"
            ky, _ = tmm.ky_three_point(own[cmax][j_fit], geo.hy[geo.j0])
            own_n, _ = normalise(own, geo, j_fit, ky)
            jlo, jhi = geo.j0 + 2, geo.Ny - geo.npml_hi - 2
            for variant in ("cont", "disc"):
                mo = os.path.join(MROOT, f"B_{case}_{pol}_{variant}")
                meep("meep_ref_transform.py", mo, ["--grid", nc.grid_path(grid), "--pol", pol, "--dt", repr(meta["dt"]),
                                                   "--variant", variant], "meep_B.json")
                with open(os.path.join(mo, "meep_B.json")) as fh:
                    mmeta = json.load(fh)
                z = np.load(os.path.join(mo, "meep_B.npz"))
                mp_prof = {}
                for c in ("Ex", "Ez", "Hx", "Hz"):
                    arr = np.full(geo.Ny + 1, np.nan + 0j)
                    arr[z[c + "_j"]] = z[c]
                    mp_prof[c] = arr
                mp_n, conj = normalise(mp_prof, geo, j_fit, ky)
                num = den = 0.0
                for c in ("Ex", "Ez", "Hx", "Hz"):
                    top = jhi if c in ("Ex", "Ez") else jhi - 1
                    a_, b_ = own_n[c][jlo:top + 1], mp_n[c][jlo:top + 1]
                    if np.nanmax(np.abs(a_)) < 1e-3:
                        continue
                    num += np.nansum(np.abs(a_ - b_) ** 2)
                    den += np.nansum(np.abs(a_) ** 2)
                rel = math.sqrt(num / den)
                rows[(case, pol, variant)] = rel
                probe_ok = all(v["position_error_cells"] < 1e-9 for v in mmeta["probe"].values())
                if variant == "cont":
                    T.row(f"3B-{1 if case == 'vac' else 2}", f"tangential DFT profiles, own vs Meep (transformation "
                          f"optics, continuous metric; {case}, {pol})", "0", f"rel L2 {rel:.2e}",
                          f"< {G3['3B-1']['rel_L2']:g}", rel < G3["3B-1"]["rel_L2"], grid=nc.grid_info(geo, meta["dt"]),
                          note=f"Meep Courant {mmeta['courant']:.4f} (= Δt_own/Δu); position probe "
                               f"{'OK' if probe_ok else 'MISMATCH'}; conj {conj}")
                    T.row("3B-4", f"both runs finite ({case}, {pol})", "finite",
                          f"own exit 0, Meep finite={mmeta['finite']}", "finite", bool(mmeta["finite"]), grid=grid)
                else:
                    T.row("3B-3", f"discrete-metric variant (s = d/Δu at nodes, h/Δu at half nodes; {case}, {pol})",
                          "—", f"rel L2 {rel:.2e} (continuous metric {rows[(case, pol, 'cont')]:.2e})", "INFO", None,
                          grid=grid)
        # metric sampling difference (3B-3)
        du = pars["du"]
        un = np.arange(geo.Ny + 1) * du
        s_node = gg.s_of_u(un, dict(du=du, amp=pars["amp"], u14=tuple(pars["u14"])))
        s_half = gg.s_of_u(un[:-1] + du / 2, dict(du=du, amp=pars["amp"], u14=tuple(pars["u14"])))
        dn = geo.dy / du
        hh = geo.hy / du
        T.row("3B-3", f"metric sampling difference ({case})", "O(Δu²)",
              f"max|s(u_j) − d_j/Δu| = {np.max(np.abs(s_node - dn)):.2e}, max|s(u_j+½) − h_j/Δu| = "
              f"{np.max(np.abs(s_half - hh)):.2e}", "INFO", None, grid=grid,
              note="our dual nodes are primal midpoints, Meep's are f(u_{j+1/2}); d_j and h_j are cell averages of s")
    return rows


def main():
    if not os.path.exists(MEEP_PY):
        print(f"BLOCKED: MEEP_PYTHON not found ({MEEP_PY})")
        T.row("3", "Meep available", "yes", "no", "—", False, note="BLOCKED")
        T.save(extra=dict(blocked=True))
        return 3
    nc.build()
    os.makedirs(FIG, exist_ok=True)
    conv = part_A()
    part_B()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for (st, pol), (ress, Rs, Rinf, Rt) in conv.items():
        ax.loglog([1 / r for r in ress], np.abs(np.array(Rs) - Rt), "o-", label=f"Meep {st} {pol} (vs TMM)")
    ax.set_xlabel("1/resolution [λ0]")
    ax.set_ylabel("|R − R_∞|")
    ax.set_title("Gate 3A: Meep convergence (thin Bloch cell)", fontsize=9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "meep_convergence.png"), dpi=130)
    plt.close(fig)
    T.save()
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
