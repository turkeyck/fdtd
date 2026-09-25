"""Theory values for the nonuniform gates, computed BEFORE any nonuniform simulation (SPEC_nonuniform §13.3).

    python3 tests/nu_predict.py            # (re)write results/nu_predictions.json (append-only for existing ids)
    python3 tests/nu_predict.py --check    # verify the file reproduces (no rewrite)

Also runs the preflight tests T-TMM1, T-TMM3, T-A1, T-A3 (exact discrete field through the forward-wave estimator)
and the Level 2 default-grid preflight (R8). Exit code 1 if a preflight fails.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import tmm  # noqa: E402
from fdtd_theory import Setup, discrete_fresnel, fresnel  # noqa: E402

REC = []


PFX = "v2/"   # D12: theory values for dt = T0/ceil(T0/dt_Courant) (the dt every nonuniform run uses)


def put(pid, quantity, value, formula, grid=None):
    g = nc.load_grid(grid) if grid else None
    REC.append(dict(id=PFX + pid, quantity=quantity, value=value, formula_ref=formula,
                    grid=grid, grid_hash=g["hash"] if g else None))


# ----------------------------------------------------------------------------- Level 1
def region(geo):
    """Solver node range for the vacuum scattering problem: uniform on both ends."""
    return geo.j0 + 2, geo.Ny - geo.npml_hi - 2


def graded_nodes(geo):
    js = set()
    for tag in ("down", "up"):
        if tag in geo.marks:
            c0, n = geo.marks[tag]
            if n:
                js |= set(range(c0, c0 + n + 1))
    return sorted(js)


def theta_series(geo, pol, dt, Kt, kt):
    """Exact discrete field on the graded grid -> forward-wave estimator -> (j, theta_meas, theta_pred, theta_cont)."""
    lo, hi = region(geo)
    res = tmm.tmm_discrete(geo.hy, geo.eps_cell, dt, Kt, pol, lo, hi)
    u, js = res["u"], res["nodes"]
    if pol == "p":                     # measured quantity is tangential E at primal nodes (as sampled by FDTD):
        jp = js[1:]                    # e_t(j) ∝ (h_s(j+1/2) - h_s(j-1/2)) / (eps_t d_j)
        u = (u[1:] - u[:-1]) / (geo.eps_t[jp] * geo.dy[jp])
        js = jp
    yy, d, ky = tmm.local_ky_forward(u, geo.y[js], dt, Kt)
    idx = js[1:-1]
    th_meas = np.arctan2(Kt, ky)
    th_pred = np.array([math.atan2(Kt, tmm.ky_local(dd, 1.0, dt, Kt)) for dd in d])
    th_cont = math.atan2(kt, tmm.ky_continuous(kt))
    return idx, yy, th_meas, th_pred, th_cont


def level1():
    th = nc.thresholds()
    kt = nc.kt_cont()
    names = list(th["gate1"]["grids"]["taper"]) + [th["gate1"]["grids"]["abrupt"], "L1_uniform_ctrl"]
    names += [f"L1_scan_r{r:g}_base{b}" for r in (1.05, 1.2, 1.5) for b in (20, 40, 80)]
    names += ["L1_scan_r4_base40", "L1_scan_r4_base80"]
    names += [f"L1_taper_r1.1_base20_N{N}" for N in (10, 20, 30, 60)]
    names += [f"L1_taper_r1.1_base{b}_xz" for b in (10, 20, 40)]          # amendment D9 (rev.)
    preflight = {}
    for name in names:
        geo = nc.Geo(nc.load_grid(name))
        dt = geo.dt_run()
        Kt = nc.Kt_of(geo)
        put(f"L1/{name}/dt", "time step dt [T0]", dt, "SPEC §7.5", name)
        put(f"L1/{name}/Kt", "|K~t| [1/λ0]", Kt, "derivation §5", name)
        lo, hi = region(geo)
        for pol in ("s", "p"):
            r = tmm.tmm_discrete(geo.hy, geo.eps_cell, dt, Kt, pol, lo, hi)
            put(f"L1/{name}/{pol}/R_disc", "exact discrete reflectance R_disc", r["R"], "derivation §5, tmm_discrete", name)
            put(f"L1/{name}/{pol}/T_disc", "exact discrete transmittance", r["T"], "derivation §5, tmm_discrete", name)
        # sub-region ky
        subs = {}
        for tag in ("tf_a", "hold", "tf_b"):
            if tag in geo.marks and geo.marks[tag][1] > 0:
                c0, n = geo.marks[tag]
                h = float(geo.hy[c0])
                ky = tmm.ky_local(h, 1.0, dt, Kt)
                subs[tag] = dict(first_cell=c0, n=n, h=h, ky=ky, theta=math.atan2(Kt, ky))
                put(f"L1/{name}/ky/{tag}", f"discrete ky in sub-region {tag} (h={h:.6g})", ky, "derivation §5", name)
        if "tf_a" in subs and "hold" in subs and subs["hold"]["h"] != subs["tf_a"]["h"]:
            ky1, ky2 = subs["tf_a"]["ky"], subs["hold"]["ky"]
            n1sq, n2sq = (kt ** 2 + ky1 ** 2) / nc.K0 ** 2, (kt ** 2 + ky2 ** 2) / nc.K0 ** 2
            L = geo.marks["hold"][1] * subs["hold"]["h"]
            for pol in ("s", "p"):
                est = tmm.tmm_continuous([(n2sq, L)], n1sq, n1sq, kt, pol=pol)["R"] if pol == "s" else \
                    tmm.tmm_continuous([(n2sq, L)], n1sq, n1sq, kt, pol=pol)["R"]
                put(f"L1/{name}/{pol}/R_fresnel_est", "numerical-index Fresnel estimate of the fine slab (INFO)",
                    est, "derivation §5 (n_num)", name)
                single = tmm.fresnel_numeric_estimate(ky1, ky2, pol="s")
                put(f"L1/{name}/{pol}/R_fresnel_single", "single-junction numerical Fresnel estimate (INFO)",
                    single, "derivation §5", name)
        put(f"L1/{name}/theta_cont", "continuous angle θ_cont [rad]", math.atan2(kt, tmm.ky_continuous(kt)),
            "derivation §5", name)
        if graded_nodes(geo) and name.startswith("L1_taper"):
            for pol in ("s", "p"):
                idx, yy, tm, tp, tc = theta_series(geo, pol, dt, Kt, kt)
                gset = set(graded_nodes(geo))
                sel = np.array([j in gset for j in idx])
                dev = np.abs(tm[sel] - tp[sel])
                off = np.abs(tp[sel] - tc)
                ratio = float(np.max(dev / off))
                preflight[(name, pol)] = dict(ratio=ratio, mean_err_cont=float(np.mean(np.abs(tm[sel] - tc))),
                                              n=int(sel.sum()))
    return preflight


# ----------------------------------------------------------------------------- Level 2
def level2():
    kt = nc.kt_cont()
    structs = {"F1": [(4.0, 0.08)], "F5": [(n * n, d) for n, d in [(2.0, 0.08), (1.46, 0.12), (2.0, 0.08),
                                                                     (1.46, 0.12), (2.0, 0.08)]]}
    out = {}
    for st, layers in structs.items():
        for pol in ("s", "p"):
            c = tmm.tmm_continuous(layers, 1.0, 1.46 ** 2, kt, pol=pol)
            put(f"L2/{st}/{pol}/R_TMM", "continuous TMM reflectance", c["R"], "tmm_continuous")
            put(f"L2/{st}/{pol}/T_TMM", "continuous TMM transmittance", c["T"], "tmm_continuous")
            out[(st, pol)] = c
    names = [f"L2_{s}_ppw{p}" for s in ("F1", "F5") for p in (20, 40, 80, 160)]
    names += [f"L2_{s}_ppw40_dxz{d}" for s in ("F1", "F5") for d in (40, 80)]
    names += [f"L2_{s}_bis{k}" for s in ("F1", "F5") for k in (1, 2, 4, 8)]      # amendment D10
    names += [f"L2_F1_ctrl_mh{m:g}" for m in (4.5, 9.5, 18.5, 37.5)]              # amendment D11
    names += [f"L2_{s}_bisxz{k}" for s in ("F1", "F5") for k in (1, 2)]           # amendment D14
    pre = {}
    for name in names:
        geo = nc.Geo(nc.load_grid(name))
        dt = geo.dt_run()
        if "_bis" in name:                                  # D10/D14: dt scales with h
            k = int(name.split("_bis")[1].lstrip("xz"))
            dt = nc.Geo(nc.load_grid(name.split("_bis")[0] + "_bis1")).dt_run() / k
        Kt = nc.Kt_of(geo)
        lo, hi = geo.j0 + 2, geo.Ny - geo.npml_hi - 2
        put(f"L2/{name}/dt", "time step dt [T0]", dt, "SPEC §7.5", name)
        for pol in ("s", "p"):
            r = tmm.tmm_discrete(geo.hy, geo.eps_cell, dt, Kt, pol, lo, hi,
                                 eps_t=geo.eps_t if "eps_t" in geo.g else None)
            put(f"L2/{name}/{pol}/R_disc", "exact discrete reflectance", r["R"], "tmm_discrete", name)
            put(f"L2/{name}/{pol}/T_disc", "exact discrete transmittance", r["T"], "tmm_discrete", name)
            st = "F5" if "F5" in name else "F1"
            pre[(name, pol)] = (r["R"] - out[(st, pol)]["R"], r["T"] - out[(st, pol)]["T"])
    return out, pre


# ----------------------------------------------------------------------------- preflight unit tests
def unit_tests():
    ok = True
    for pol in ("s", "p"):
        st = Setup(n_lambda=20, pol=pol)
        stm = Setup(n_lambda=20, pol=pol, n_ref=1.5)
        D = discrete_fresnel(st, stm, 0.5 * (1 + 2.25))
        Kt = math.hypot(tmm.ktilde(st.kx, st.Delta), tmm.ktilde(st.kz, st.Delta))
        h = np.full(200, st.Delta)
        e = np.ones(200)
        e[100:] = 2.25
        r = tmm.tmm_discrete(h, e, st.dt, Kt, pol, 50, 150)
        d1 = abs(r["R"] - D["R"])
        print(f"T-TMM1 ({pol}): |R_new - R_old_discrete| = {d1:.2e} (< 1e-12)")
        ok &= d1 < 1e-12
        kt = math.hypot(st.kx, st.kz)
        c1 = math.sqrt(1 - (kt / nc.K0) ** 2)
        c2 = math.sqrt(1 - (kt / (nc.K0 * 1.5)) ** 2)
        F = fresnel(1, 1.5, c1, c2)
        rc = tmm.tmm_continuous([], 1.0, 2.25, kt, pol=pol)
        d3 = abs(rc["R"] - F["R" + pol])
        print(f"T-TMM3 ({pol}): |R_TMM - Fresnel| = {d3:.2e} (< 1e-14)")
        ok &= d3 < 1e-14
    # T-A1: three-point estimator on a synthetic two-wave field
    h, ky = 0.05, 5.03
    y = np.arange(60) * h
    u = 1.0 * np.exp(1j * ky * y) + (0.3 - 0.2j) * np.exp(-1j * ky * y)
    kym, im = tmm.ky_three_point(u, h)
    print(f"T-A1: three-point ky error = {abs(kym - ky):.2e} (< 1e-13)")
    ok &= abs(kym - ky) < 1e-13
    return ok


def stageB():
    """RCWA reference for B3 (D15 geometry) and TMM for the D16 families (classical incidence, θ = 30°)."""
    import rcwa
    import grid_gen as gg
    ok = True
    kx, kz = nc.kvec(1, 1)
    kt = math.hypot(kx, kz)
    k = np.array([kx, math.sqrt(nc.K0 ** 2 - kt ** 2), kz])
    for pol in ("s", "p"):
        s = np.cross(k, [0, 1, 0]); s /= np.linalg.norm(s)
        pp = np.cross(s, k); pp /= np.linalg.norm(pp)
        e = s if pol == "s" else pp
        prev = None
        for M in (20, 40, 80, 160):
            g = rcwa.solve(M, 2.0, 0.3, 1.0, gg.N_SUB ** 2, gg.N_FILM ** 2, 1.0, 0.5, kx, kz, e)
            pr, pt = rcwa.propagating(g, "R"), rcwa.propagating(g, "T")
            v = np.concatenate([g["R_orders"][pr], g["T_orders"][pt]])
            if prev is not None:
                put(f"B/rcwa/{pol}/maxdiff_{2 * M + 1}", f"max efficiency change to {2 * M + 1} orders",
                    float(np.max(np.abs(v - prev))), "rcwa.py")
            prev = v
        put(f"B/rcwa/{pol}/orders_R", "propagating reflected orders", [int(q) for q in g["orders"][pr]], "rcwa.py")
        put(f"B/rcwa/{pol}/orders_T", "propagating transmitted orders", [int(q) for q in g["orders"][pt]], "rcwa.py")
        put(f"B/rcwa/{pol}/R_p", "RCWA reflected efficiencies (321 orders)", g["R_orders"][pr].tolist(), "rcwa.py")
        put(f"B/rcwa/{pol}/T_p", "RCWA transmitted efficiencies (321 orders)", g["T_orders"][pt].tolist(), "rcwa.py")
        last = REC_value(f"B/rcwa/{pol}/maxdiff_321")
        flag = last < 1e-5
        ok &= flag
        print(f"  B3-3 RCWA ({pol}): 161->321 orders max change {last:.2e} {'OK' if flag else 'FAIL'}; "
              f"sum R+T-1 = {g['R'] + g['T'] - 1:+.1e}")
    kt0 = math.pi                                     # m = 1, n = 0, Lx = 2
    for pol in ("s", "p"):
        c = tmm.tmm_continuous([(gg.N_FILM ** 2, gg.D_FILM)], 1.0, gg.N_SUB ** 2, kt0, pol=pol)
        put(f"B/film/{pol}/R_TMM", "TMM R, film on substrate, θ = 30° (D16 families)", c["R"], "tmm_continuous")
        put(f"B/film/{pol}/T_TMM", "TMM T, film on substrate, θ = 30°", c["T"], "tmm_continuous")
    return ok


def main():
    check_only = "--check" in sys.argv
    ok = unit_tests()
    pf1 = level1()
    l2, pf2 = level2()
    print("\nT-A3 preflight (exact discrete field -> forward-wave estimator, D3: max |θm-θp|/|θp-θc| <= 0.1):")
    errs = {}
    for (name, pol), v in pf1.items():
        flag = v["ratio"] <= 0.1
        ok &= flag
        print(f"  {name:<24} {pol}: max ratio {v['ratio']:.3e} over {v['n']} graded nodes, "
              f"mean |θm-θc| = {v['mean_err_cont']:.3e} rad -> {'OK' if flag else 'FAIL'}")
        if name in ("L1_taper_r1.1_base10_xz", "L1_taper_r1.1_base20_xz", "L1_taper_r1.1_base40_xz"):
            errs.setdefault(pol, []).append((1.0 / int(name.split("base")[1].split("_")[0]), v["mean_err_cont"]))
    for pol, e in errs.items():
        e.sort()
        o = nc.order_fit([q[0] for q in e], [q[1] for q in e])
        put(f"L1/theta_order_D9rev/{pol}", "predicted order of mean |θ_meas-θ_cont| (exact discrete field)", o,
            "derivation §5")
        print(f"  predicted 1-3b order ({pol}): {o:.3f} (window [1.8, 2.2])")
    print("\nLevel 2 preflight (R8): default grid D6 and series, exact discrete vs continuous TMM:")
    for (name, pol), (dR, dT) in pf2.items():
        mark = ""
        if name.endswith("_ppw40"):
            flag = abs(dR) < 1e-3 and abs(dT) < 1e-3
            ok &= flag
            mark = "  <- default grid (2-1) " + ("OK" if flag else "PREDICTED FAIL")
        print(f"  {name:<26} {pol}: R_disc-R_TMM = {dR:+.3e}, T_disc-T_TMM = {dT:+.3e}{mark}")
    # orders predicted for 2-3 (D10 bisection family, vs Richardson limit) and 2-4 (D11 control vs same limit)
    rinf = {}
    for st in ("F1", "F5"):
        for pol in ("s", "p"):
            hs = [1.0, 0.5, 0.25, 0.125]
            Rs = [REC_value(f"L2/L2_{st}_bis{k}/{pol}/R_disc") for k in (1, 2, 4, 8)]
            Rinf = Rs[-1] + (Rs[-1] - Rs[-2]) / 3.0
            rinf[(st, pol)] = Rinf
            put(f"L2/Rinf/{st}/{pol}", "Richardson limit of R_disc (bisection family, dx fixed)", Rinf, "tmm_discrete")
            o = nc.order_fit(hs[:3], np.array(Rs[:3]) - Rinf)
            put(f"L2/order_y/{st}/{pol}", "predicted order of R vs Δy (D10 family)", o, "tmm_discrete")
            flag = 1.8 <= o <= 2.2
            ok &= flag
            print(f"  predicted 2-3 order {st} {pol}: {o:.3f} {'OK' if flag else 'FAIL'}  "
                  f"(R_inf - R_TMM = {Rinf - l2[(st, pol)]['R']:+.2e})")
    for st in ("F1", "F5"):                     # D10 revision: 3-level estimate on k = 1, 2, 4
        for pol in ("s", "p"):
            R3 = [REC_value(f"L2/L2_{st}_bis{k}/{pol}/R_disc") for k in (1, 2, 4)]
            p3 = float(np.log2((R3[0] - R3[1]) / (R3[1] - R3[2])))
            rinf3 = R3[2] + (R3[2] - R3[1]) / (2 ** p3 - 1)
            put(f"L2/order_y3/{st}/{pol}", "predicted 3-level order (k=1,2,4)", p3, "tmm_discrete")
            put(f"L2/Rinf3/{st}/{pol}", "3-level Richardson limit", rinf3, "tmm_discrete")
            flag = 1.8 <= p3 <= 2.2
            ok &= flag
            print(f"  predicted 2-3 order (3-level) {st} {pol}: {p3:.3f} {'OK' if flag else 'FAIL'}")
            if st == "F1":
                mh3 = (4.5, 9.5, 18.5)
                C3 = [REC_value(f"L2/L2_F1_ctrl_mh{m:g}/{pol}/R_disc") for m in mh3]
                o3 = nc.order_fit([0.08 / m for m in mh3], np.array(C3) - rinf3)
                put(f"L2/order_control3/{pol}", "predicted control order (mh 4.5, 9.5, 18.5) vs 3-level R_inf", o3,
                    "tmm_discrete")
                flag = 0.8 <= o3 <= 1.2
                ok &= flag
                print(f"  predicted 2-4 order (control, 3 levels) {pol}: {o3:.3f} {'OK' if flag else 'FAIL'}")
    for st in ("F1", "F5"):                     # D14: 2-level combined Richardson (assumed order 2)
        for pol in ("s", "p"):
            R1, R2 = (REC_value(f"L2/L2_{st}_bisxz{k}/{pol}/R_disc") for k in (1, 2))
            ri = R2 + (R2 - R1) / 3.0
            put(f"L2/Rinf_xz2/{st}/{pol}", "2-level combined Richardson limit (D14)", ri, "tmm_discrete")
            flag = abs(ri - l2[(st, pol)]["R"]) < 1e-4
            ok &= flag
            print(f"  predicted 3A-3 limit {st} {pol}: R_inf - R_TMM = {ri - l2[(st, pol)]['R']:+.2e} {'OK' if flag else 'FAIL'}")
    for pol in ("s", "p"):
        mhs = (4.5, 9.5, 18.5, 37.5)
        hs = [0.08 / m for m in mhs]
        Rs = [REC_value(f"L2/L2_F1_ctrl_mh{m:g}/{pol}/R_disc") for m in mhs]
        o = nc.order_fit(hs, np.array(Rs) - rinf[("F1", pol)])
        put(f"L2/order_control/{pol}", "predicted order, unaligned uniform control (D11)", o, "tmm_discrete")
        flag = 0.8 <= o <= 1.2
        ok &= flag
        print(f"  predicted 2-4 order (control) {pol}: {o:.3f} {'OK' if flag else 'FAIL'}")
    for pol in ("s", "p"):
        o = REC_value(f"L1/theta_order_D9rev/{pol}")
        flag = 1.8 <= o <= 2.2
        ok &= flag
    print("\nStage B references:")
    ok &= stageB()
    for r in REC:
        r["created_utc"] = nc.utcnow()
    old = nc.load_predictions()
    if check_only:
        worst = 0.0
        for r in REC:
            if r["id"] in old and isinstance(r["value"], float):
                worst = max(worst, abs(r["value"] - old[r["id"]]["value"]) / max(abs(r["value"]), 1e-300))
        print(f"--check: max relative change vs stored predictions = {worst:.1e}")
        return 0 if ok else 1
    merged = dict(old)
    for r in REC:
        if r["id"] in old and old[r["id"]].get("grid_hash") == r["grid_hash"]:
            continue                                    # append-only: never overwrite an existing id
        if r["id"] in old:
            r["id"] = r["id"] + f"#{r['grid_hash'][:8]}"
        merged[r["id"]] = r
    os.makedirs(nc.RESULTS, exist_ok=True)
    with open(nc.PRED_FILE, "w", encoding="utf-8") as fh:
        json.dump(dict(doc="Frozen theory values (SPEC_nonuniform §13.3). Append-only.",
                       records=list(merged.values())), fh, indent=1, ensure_ascii=False)
    print(f"\n{len(merged)} predictions in {nc.PRED_FILE}; preflight {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def REC_value(pid):
    for r in REC:
        if r["id"] == PFX + pid:
            return r["value"]
    raise KeyError(pid)


if __name__ == "__main__":
    sys.exit(main())
