"""Physical parameters in, FDTD simulation out: one entry point for uniform and nonuniform grids.

    python3 simulate.py examples/film_nonuniform.json             # run every wavelength x polarization
    python3 simulate.py examples/film_nonuniform.json --dry-run   # grid, realized angle, cost; no run
    python3 simulate.py --example film > my.json                  # print a template (film|stack|uniform|spectrum|grating)

The configuration is a JSON file of physical quantities (wavelength, angles, layer stack, mesh density); lengths
are strings with a unit ("633nm", "0.08um") or bare numbers in nm. The script builds the grid (grid_gen.py),
chooses the time step, runs the C solver (fdtd3d_oblique), and measures R, T (planar stacks) or per-order
diffraction efficiencies (one lamellar grating layer), with references: continuous TMM and the exact discrete
(Yee-lattice) TMM for stacks, RCWA for a single grating layer on a substrate. See docs/USER_GUIDE.md.

Limits of the solver (checked here): incident medium vacuum; lossless, non-dispersive, non-magnetic media with
n >= 1; 0 < θ < 90° and sin θ < 0.95 (the s/p basis excludes normal incidence); periodic x/z without Bloch phase,
so kx = 2πm/Lx and kz = 2πn/Lz are quantized (the realized angle is reported); at most one grating layer, one
ridge per period, grating period = Lx.
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))
import analyze_nu as an  # noqa: E402
import grid_gen as gg  # noqa: E402
import nu_common as nc  # noqa: E402
import nu_runs as RN  # noqa: E402
import tmm  # noqa: E402
from oblique import C_LIGHT, parse_length, quantize  # noqa: E402

BIN = os.path.join(ROOT, "fdtd3d_oblique")
CELL_STEPS_PER_S = 6e7          # measured field-update throughput, 8 threads, idle machine; the DFT outputs add
                                # ~50% and at least ~1 ms per step (estimate only)

EXAMPLES = {
    "film": {
        "name": "film_nonuniform",
        "wavelength": "633nm", "angle_deg": 30, "azimuth_deg": 0, "polarizations": ["s", "p"],
        "layers": [{"name": "TiO2", "thickness": "80nm", "n": 2.0}],
        "substrate_n": 1.46,
        "mesh": {"type": "nonuniform", "ppw": 40, "transverse": "lambda/20", "r_max": 1.1},
        "run": {"periods": 90, "dft_periods": 30, "courant": 0.5},
    },
    "uniform": {
        "name": "film_uniform",
        "wavelength": "633nm", "angle_deg": 30, "azimuth_deg": 0, "polarizations": ["s", "p"],
        "layers": [{"name": "TiO2", "thickness": "80nm", "n": 2.0}],
        "substrate_n": 1.46,
        "mesh": {"type": "uniform", "ppw": 20},
        "run": {"periods": 90, "dft_periods": 30, "courant": 0.5},
    },
    "stack": {
        "name": "stack5",
        "wavelength": "1000nm", "angle_deg": 36.94, "azimuth_deg": 33.69, "polarizations": ["s", "p"],
        "layers": [{"thickness": "80nm", "n": 2.0}, {"thickness": "120nm", "n": 1.46}, {"thickness": "80nm", "n": 2.0},
                   {"thickness": "120nm", "n": 1.46}, {"thickness": "80nm", "n": 2.0}],
        "substrate_n": 1.46,
        "mesh": {"type": "nonuniform", "ppw": 40, "transverse": "lambda/20"},
    },
    "spectrum": {
        "name": "film_spectrum",
        "wavelength": ["500nm", "550nm", "600nm", "650nm", "700nm"], "angle_deg": 30, "polarizations": ["s"],
        "layers": [{"thickness": "80nm", "n": 2.0}],
        "substrate_n": 1.46,
        "mesh": {"type": "nonuniform", "ppw": 30, "transverse": "lambda/20"},
    },
    "grating": {
        "name": "lamellar_grating",
        "wavelength": "1000nm", "angle_deg": 36.94, "azimuth_deg": 33.69, "polarizations": ["s", "p"],
        "layers": [{"name": "grating", "thickness": "300nm",
                    "grating": {"period": "2000nm", "duty": 0.5, "n_ridge": 2.0, "n_groove": 1.0}}],
        "substrate_n": 1.46,
        "mesh": {"type": "nonuniform", "ppw": 20, "transverse": "lambda/20"},
        "run": {"periods": 200, "dft_periods": 60},
    },
}

DEFAULTS = {
    "azimuth_deg": 0.0, "polarizations": ["s", "p"], "substrate_n": 1.0, "layers": [],
    "mesh": {"type": "nonuniform", "ppw": 40, "transverse": None, "r_max": 1.1},
    "run": {"periods": 90, "dft_periods": 30, "courant": 0.5},
    "geometry": {"sf": 1.0, "tf": 3.0, "substrate": 3.0, "pml_cells": 20},
    "angle": {"tol_deg": 0.05, "max_period_cells": 400},
    "output": {"dir": None, "field_map": True},
}


# ----------------------------------------------------------------------------- configuration
def load_config(path):
    with open(path, encoding="utf-8") as fh:
        user = json.load(fh)
    cfg = json.loads(json.dumps(DEFAULTS))
    for k, v in user.items():
        if isinstance(v, dict) and isinstance(cfg.get(k), dict):
            cfg[k].update(v)
        else:
            cfg[k] = v
    cfg.setdefault("name", os.path.splitext(os.path.basename(path))[0])
    if cfg["output"]["dir"] is None:
        cfg["output"]["dir"] = os.path.join("runs", "user", cfg["name"])
    validate(cfg)
    return cfg


def validate(cfg):
    if "wavelength" not in cfg:
        raise ValueError("'wavelength' is required")
    if not 0.0 < float(cfg["angle_deg"]) < 90.0:
        raise ValueError("angle_deg must be in (0, 90): normal incidence is not supported (s/p basis)")
    for p in cfg["polarizations"]:
        if p not in ("s", "p"):
            raise ValueError(f"polarization '{p}' (use s and/or p)")
    if cfg["mesh"]["type"] not in ("uniform", "nonuniform"):
        raise ValueError("mesh.type must be 'uniform' or 'nonuniform'")
    if cfg["mesh"]["ppw"] < 10:
        raise ValueError("mesh.ppw < 10 is not supported")
    ns = [cfg["substrate_n"]]
    ng = 0
    for q, L in enumerate(cfg["layers"]):
        if "thickness" not in L:
            raise ValueError(f"layer {q}: 'thickness' missing")
        if "grating" in L:
            ng += 1
            G = L["grating"]
            for key in ("period", "duty", "n_ridge"):
                if key not in G:
                    raise ValueError(f"layer {q}: grating.{key} missing")
            if not 0.0 < G["duty"] < 1.0:
                raise ValueError(f"layer {q}: grating.duty must be in (0, 1)")
            ns += [G["n_ridge"], G.get("n_groove", 1.0)]
        else:
            if "n" not in L:
                raise ValueError(f"layer {q}: 'n' missing")
            ns.append(L["n"])
    if ng > 1:
        raise ValueError("at most one grating layer is supported")
    if min(ns) < 1.0:
        raise ValueError("refractive indices must be >= 1 (lossless, non-dispersive, non-magnetic media)")


def wavelengths(cfg):
    w = cfg["wavelength"]
    return [parse_length(v if isinstance(v, str) else str(v)) for v in (w if isinstance(w, list) else [w])]


def spacing(s, lam, default):
    """'lambda/N' -> 1/N; a length -> length/λ (in λ0 units); None -> default."""
    if s is None:
        return default
    s = str(s).strip()
    if s.startswith(("lambda/", "λ/")):
        return 1.0 / float(s.split("/")[1])
    return parse_length(s) / lam


# ----------------------------------------------------------------------------- grid construction
def layer_specs(cfg, lam):
    """Stack in λ0 units: dicts(d, eps (background), nmax (densest material), grating or None, name)."""
    out = []
    for q, L in enumerate(cfg["layers"]):
        d = parse_length(str(L["thickness"])) / lam
        if "grating" in L:
            G = L["grating"]
            nr, ngv = float(G["n_ridge"]), float(G.get("n_groove", 1.0))
            out.append(dict(d=d, eps=ngv ** 2, nmax=max(nr, ngv), name=L.get("name", f"layer{q}"),
                            grating=dict(period=parse_length(str(G["period"])) / lam, duty=float(G["duty"]),
                                         eps_ridge=nr ** 2, eps_groove=ngv ** 2)))
        else:
            out.append(dict(d=d, eps=float(L["n"]) ** 2, nmax=float(L["n"]), name=L.get("name", f"layer{q}"),
                            grating=None))
    return out


def build(cfg, lam, vacuum=False):
    """Grid dict + plan info for one wavelength. vacuum=True: same spacings, vacuum everywhere, no grating."""
    mesh, geo_c, ang = cfg["mesh"], cfg["geometry"], cfg["angle"]
    ppw, r_max = float(mesh["ppw"]), float(mesh.get("r_max", 1.1))
    npml = int(geo_c["pml_cells"])
    ns = float(cfg["substrate_n"])
    specs = layer_specs(cfg, lam)
    uniform = mesh["type"] == "uniform"
    nmax = max([1.0, ns] + [s["nmax"] for s in specs])
    D = 1.0 / (ppw * nmax)                                     # uniform mesh: one spacing everywhere
    info = dict(mesh=mesh["type"], layers=[])
    # ---- y blocks
    hv = D if uniform else 1.0 / ppw
    hs = D if uniform else 1.0 / (ns * ppw)
    blocks = [dict(kind="uniform", n=npml, h=hv, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=int(round(geo_c["sf"] / hv)), h=hv, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(geo_c["tf"] / hv)), h=hv, eps=1.0, tag="tf_vac")]
    if not uniform:                                            # vacuum is the coarsest medium: grade on its side
        blocks.append(dict(kind="grade", eps=1.0, tag="grade_vac"))
    for q, s in enumerate(specs):
        if uniform:
            k = int(round(s["d"] / D))
            if k < 1:
                raise ValueError(f"layer '{s['name']}' ({s['d'] * lam * 1e9:.3g} nm) is thinner than Δ/2 "
                                 f"({D * lam * 1e9 / 2:.3g} nm): raise mesh.ppw or use mesh.type = nonuniform")
            blocks.append(dict(kind="uniform", n=k, h=D, eps=s["eps"], tag=f"layer{q}"))
            info["layers"].append(dict(name=s["name"], requested_nm=s["d"] * lam * 1e9, realized_nm=k * D * lam * 1e9))
        else:
            blocks.append(dict(kind="layer", L=s["d"], eps=s["eps"], hmax=1.0 / (s["nmax"] * ppw), tag=f"layer{q}"))
            info["layers"].append(dict(name=s["name"], requested_nm=s["d"] * lam * 1e9, realized_nm=s["d"] * lam * 1e9))
    if specs and not uniform:
        last = specs[-1]
        hlast = last["d"] / math.ceil(last["d"] * last["nmax"] * ppw - 1e-12)
        if hs > hlast:        # coarser substrate: grade on the substrate side; a finer one is graded inside the layer
            blocks.append(dict(kind="grade", eps=ns ** 2, tag="grade_sub"))
    blocks += [dict(kind="uniform", n=int(round(geo_c["substrate"] / hs)), h=hs, eps=ns ** 2, tag="sub"),
               dict(kind="uniform", n=npml, h=hs, eps=ns ** 2, tag="pml_hi")]
    try:
        h, eps, marks = gg.build_y(blocks, r_max)
    except ValueError as e:
        raise ValueError(f"grid: {e} -- raise mesh.ppw or mesh.r_max") from None
    j0 = npml + blocks[1]["n"]
    # ---- transverse axes and angle
    th, ph = math.radians(cfg["angle_deg"]), math.radians(cfg["azimuth_deg"])
    sx, sz = math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph)
    ktol = math.cos(th) * math.radians(ang["tol_deg"]) / math.sqrt(2.0)
    dT = spacing(mesh.get("transverse"), lam, D if uniform else 1.0 / 20)
    grat = next((s["grating"] for s in specs if s["grating"]), None)
    gq = next((q for q, s in enumerate(specs) if s["grating"]), None)
    if grat is None:
        m, Nx = quantize(sx, 1.0 / dT, ang["max_period_cells"], ktol)
        x = gg.uniform_axis(Nx * dT, dT)
        i_hi = None
    else:
        Lam = grat["period"]
        m = int(round(sx * Lam))                                  # kx = 2πm/Λ: the period is the x cell
        if uniform:
            Nx = int(round(Lam / D))
            Nr = int(round(grat["duty"] * Nx))
            if not 0 < Nr < Nx:
                raise ValueError("grating ridge or groove narrower than one cell: raise mesh.ppw")
            x = gg.uniform_axis(Lam, Lam / Nx)
            i_hi = Nr
            info["grating_duty_realized"] = Nr / Nx
        else:
            x, starts = gg.periodic_axis_layers([(grat["duty"] * Lam, 1.0 / (math.sqrt(grat["eps_ridge"]) * ppw)),
                                                 ((1 - grat["duty"]) * Lam, 1.0 / (math.sqrt(grat["eps_groove"]) * ppw))],
                                                r_max)
            i_hi = starts[1]
            info["grating_duty_realized"] = grat["duty"]
    n, Nz = quantize(sz, 1.0 / dT, ang["max_period_cells"], ktol)
    z = gg.uniform_axis(Nz * dT, dT)
    if m == 0 and n == 0:
        raise ValueError("the realized tangential wavevector is zero (normal incidence is not supported); "
                         "change the angle or the grating period")
    Lx, Lz = x["L"], z["L"]
    kxr, kzr = m / Lx, n / Lz                                     # in units of k0
    kt = math.hypot(kxr, kzr)
    if kt >= 0.95:
        raise ValueError(f"realized sin θ = {kt:.4f} >= 0.95 (the solver needs a 5% propagation margin)")
    info.update(requested_angle_deg=cfg["angle_deg"], requested_azimuth_deg=cfg["azimuth_deg"],
                angle_deg=math.degrees(math.asin(kt)), azimuth_deg=math.degrees(math.atan2(kzr, kxr)),
                m=m, n=n, Lx_nm=Lx * lam * 1e9, Lz_nm=Lz * lam * 1e9)
    if vacuum:
        eps = np.ones_like(eps)
    name = f"{cfg['name']}_{lam * 1e9:.6g}nm" + ("_vac" if vacuum else "")
    g = gg.finish(name, h, eps, x, z, r_max if not uniform else 1.1, npml, npml, j0, 40,
                  None if uniform else ppw, "user_" + mesh["type"], dict(wavelength_nm=lam * 1e9), marks)
    if grat is not None and not vacuum:
        c0, nc_ = marks[f"layer{gq}"]
        g["grating"] = {"j_lo": c0, "j_hi": c0 + nc_, "i_lo": 0, "i_hi": i_hi, "eps_ridge": grat["eps_ridge"],
                        "eps_groove": grat["eps_groove"], "thickness": specs[gq]["d"]}
        g["hash"] = gg.canonical_hash(g)
    bad = gg.check(g)
    if bad:
        raise ValueError(f"grid check failed: {bad}")
    geo = nc.Geo(g)
    dt = geo.dt_run(cfg["run"]["courant"])
    Np = int(round(1.0 / dt))
    cells = geo.Nx * (geo.Ny + 1) * geo.Nz
    steps = cfg["run"]["periods"] * Np
    info.update(Nx=geo.Nx, Ny=geo.Ny, Nz=geo.Nz, cells=cells, steps=steps, steps_per_period=Np,
                dy_min_nm=float(geo.hy.min() * lam * 1e9), dy_max_nm=float(geo.hy.max() * lam * 1e9),
                dx_min_nm=float(geo.hx.min() * lam * 1e9), dz_nm=float(geo.hz[0] * lam * 1e9),
                r_max=float(np.max(np.maximum(geo.hy[1:] / geo.hy[:-1], geo.hy[:-1] / geo.hy[1:]))),
                dt_fs=dt * lam / C_LIGHT * 1e15, est_minutes=steps * max(1.5 * cells / CELL_STEPS_PER_S, 1e-3) / 60.0,
                inc="a" if geo.uniform_xz() else "p", grating=grat is not None,
                x_nonuniform=not bool(g["x"]["uniform"]))
    return g, geo, dt, info


# ----------------------------------------------------------------------------- running
def run_solver(outdir, g, **kw):
    """Run the C solver on grid dict g; reuse the run if the arguments, the grid and the binary are unchanged."""
    os.makedirs(outdir, exist_ok=True)
    gpath = os.path.join(outdir, "grid.json")
    with open(gpath, "w") as fh:
        json.dump(g, fh)
    argv = [f"{k}={','.join(map(str, v)) if isinstance(v, (list, tuple)) else v}" for k, v in kw.items()]
    key = dict(argv=argv, grid_hash=g["hash"])
    argfile, metafile = os.path.join(outdir, "args.json"), os.path.join(outdir, "meta.json")
    if os.path.exists(argfile) and os.path.exists(metafile) and os.path.getmtime(metafile) >= os.path.getmtime(BIN):
        with open(argfile) as fh:
            if json.load(fh) == key:
                return nc.load_meta(outdir)
    env = dict(os.environ)
    if len(g["x"]["h"]) * len(g["z"]["h"]) <= 64:
        env["OMP_NUM_THREADS"] = "1"
    res = subprocess.run([BIN, f"grid={gpath}"] + argv + [f"out={outdir}"], capture_output=True, text=True, env=env)
    if res.returncode != 0:
        raise RuntimeError(f"solver failed ({res.returncode}):\n{res.stderr}")
    with open(argfile, "w") as fh:
        json.dump(key, fh)
    return nc.load_meta(outdir)


def solver_args(cfg, info, dt, pol, geo, planes=None, field_map=False):
    Np = info["steps_per_period"]
    P, W = cfg["run"]["periods"], cfg["run"]["dft_periods"]
    kw = dict(pol=pol, inc=info["inc"], m=info["m"], n=info["n"], proj=1, dt=repr(dt), energy_every=0,
              nsteps=P * Np, dft0=(P - W) * Np, dft1=P * Np)
    if planes:
        kw["yplanes"] = planes
    if field_map:
        kw["zk"] = 0
    return kw


# ----------------------------------------------------------------------------- planar stack: R, T
def uniform_run_from(geo, j_start, direction):
    h0, e0 = geo.hy[j_start], geo.eps_cell[j_start]
    j = j_start
    while 0 <= j + direction < geo.Ny and geo.hy[j + direction] == h0 and geo.eps_cell[j + direction] == e0:
        j += direction
    return (j_start, j) if direction > 0 else (j, j_start)


def measure_planar(out, meta, geo):
    """R: forward/backward fit of the Floquet-projected tangential E in the uniform vacuum region after the TF/SF
    plane; T: conserved S_y in the uniform substrate region / exact incident flux (same as gate 2)."""
    used = nc.load_used(out)
    a = RN.load_proj(out, meta)
    c = "Ex" if np.nanmax(np.abs(a["Ex"])) >= np.nanmax(np.abs(a["Ez"])) else "Ez"
    _, c1 = uniform_run_from(geo, geo.j0, +1)
    js = np.arange(geo.j0 + 1, c1 + 2)
    u = a[c][js]
    ky, _ = tmm.ky_three_point(u[1:], float(geo.hy[geo.j0]))
    A = np.stack([np.exp(1j * ky * used["y"][js]), np.exp(-1j * ky * used["y"][js])], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    R = abs(sol[1] / sol[0]) ** 2
    s0, s1 = uniform_run_from(geo, geo.Ny - geo.npml_hi - 1, -1)
    yv, yd = used["y"], used["y_dual"]

    def S_at(j):
        wa = (yd[j] - yv[j]) / (yd[j] - yd[j - 1])
        wb = (yv[j] - yd[j - 1]) / (yd[j] - yd[j - 1])
        Hx = wa * a["Hx"][j - 1] + wb * a["Hx"][j]
        Hz = wa * a["Hz"][j - 1] + wb * a["Hz"][j]
        return 0.5 * float(np.real(a["Ez"][j] * np.conj(Hx) - a["Ex"][j] * np.conj(Hz)))
    Ssub = np.array([S_at(j) for j in range(s0 + 1, s1 + 1)])
    E0, H0 = np.array(meta["E0"]), np.array(meta["H0"])
    Sinc = 0.5 * math.cos(meta["ky"] * meta["Delta_a"] / 2) * (E0[2] * H0[0] - E0[0] * H0[2])
    return dict(R=float(R), T=float(Ssub.mean() / Sinc), flux_spread=float(np.ptp(Ssub) / abs(Sinc)))


def references_planar(cfg, lam, meta, geo, dt, info, pol):
    kx, kz = meta["kx"], meta["kz"]
    ns = float(cfg["substrate_n"])
    layers = [(s["eps"], L["realized_nm"] * 1e-9 / lam) for s, L in zip(layer_specs(cfg, lam), info["layers"])]
    cont = tmm.tmm_continuous(layers, 1.0, ns ** 2, math.hypot(kx, kz), pol=pol)
    Kt = math.hypot(tmm.ktilde(kx, geo.hx[0]), tmm.ktilde(kz, geo.hz[0]))
    disc = tmm.tmm_discrete(geo.hy, geo.eps_cell, dt, Kt, pol, geo.j0 + 2, geo.Ny - geo.npml_hi - 2)
    out = dict(R_tmm=cont["R"], T_tmm=cont["T"], R_disc=float(disc["R"]), T_disc=float(disc["T"]))
    if any(abs(L["realized_nm"] - L["requested_nm"]) > 1e-9 for L in info["layers"]):
        req = [(s["eps"], s["d"]) for s in layer_specs(cfg, lam)]    # uniform mesh: thicknesses were snapped
        c2 = tmm.tmm_continuous(req, 1.0, ns ** 2, math.hypot(kx, kz), pol=pol)
        out.update(R_tmm_requested=c2["R"], T_tmm_requested=c2["T"])
    return out


# ----------------------------------------------------------------------------- grating: per-order efficiencies
def grating_planes(geo):
    jR = geo.npml_lo + (geo.j0 - geo.npml_lo) // 2
    t0, tn = geo.marks["tf_vac"]
    s0, sn = geo.marks["sub"]
    return jR, t0 + tn // 2, s0 + sn // 2


def measure_grating(cfg, lam, out, meta, out_v, meta_v, geo, info, pol):
    used = nc.load_used(out)
    Lx = used["x"][-1]
    ns = float(cfg["substrate_n"])
    jR, jI, jT = grating_planes(geo)
    kx0, kz = meta["kx"], meta["kz"]
    P = int(math.ceil(Lx * (1 + ns))) + 2
    orders = np.arange(-P, P + 1)
    kxs = kx0 + 2 * math.pi * orders / Lx
    if meta_v is None:                                            # inc = a: exact discrete incident flux
        E0, H0 = np.array(meta["E0"]), np.array(meta["H0"])
        S_inc = S_inc_d = 0.5 * math.cos(meta["ky"] * meta["Delta_a"] / 2) * (E0[2] * H0[0] - E0[0] * H0[2])
    else:                                                         # inc = p: vacuum normalization run
        Dv, _ = nc.load_dft(out_v, f"y{jI}", meta_v)
        Sv, _ = an.order_fluxes(Dv, used, kxs, kz)
        S_inc, S_inc_d = Sv[P], an.plane_flux(Dv, used)
    DR, _ = nc.load_dft(out, f"y{jR}", meta)
    DT, _ = nc.load_dft(out, f"y{jT}", meta)
    SR, _ = an.order_fluxes(DR, used, kxs, kz)
    ST, _ = an.order_fluxes(DT, used, kxs, kz)
    k0 = 2 * math.pi
    propR = kxs ** 2 + kz ** 2 < k0 ** 2
    propT = kxs ** 2 + kz ** 2 < (k0 * ns) ** 2
    R = {int(p): float(-SR[q] / S_inc) for q, p in enumerate(orders) if propR[q]}
    T = {int(p): float(ST[q] / S_inc) for q, p in enumerate(orders) if propT[q]}
    energy = (an.plane_flux(DT, used) - an.plane_flux(DR, used)) / S_inc_d
    res = dict(R_orders=R, T_orders=T, R=sum(R.values()), T=sum(T.values()), energy_conserved_flux=float(energy))
    specs = layer_specs(cfg, lam)
    if len(specs) == 1:                                           # RCWA: one lamellar layer on the substrate
        import rcwa
        G = specs[0]["grating"]
        k = np.array([kx0, math.sqrt(k0 ** 2 - kx0 ** 2 - kz ** 2), kz])
        s = np.cross(k, [0, 1, 0])
        s /= np.linalg.norm(s)
        pp = np.cross(s, k)
        pp /= np.linalg.norm(pp)
        duty = info.get("grating_duty_realized", G["duty"])
        d_real = info["layers"][0]["realized_nm"] * 1e-9 / lam              # uniform mesh: snapped thickness
        rc = rcwa.solve(80, Lx, d_real, 1.0, ns ** 2, G["eps_ridge"], G["eps_groove"], duty, kx0, kz,
                        s if pol == "s" else pp)
        ro = {int(p): float(v) for p, v in zip(rc["orders"], rc["R_orders"])}
        to = {int(p): float(v) for p, v in zip(rc["orders"], rc["T_orders"])}
        res["rcwa_R_orders"] = {p: ro[p] for p in R}
        res["rcwa_T_orders"] = {p: to[p] for p in T}
        res["max_abs_diff_vs_rcwa"] = max([abs(R[p] - ro[p]) for p in R] + [abs(T[p] - to[p]) for p in T])
    return res


# ----------------------------------------------------------------------------- figures
def field_map(out, meta, geo, title, path, lam):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    used = nc.load_used(out)
    F, _ = nc.load_dft(out, "xy", meta)
    live = [c for c in ("Ex", "Ey", "Ez") if np.nanmax(np.abs(F[c])) > 1e-6]
    Lx = used["x"][-1]
    tiles = max(1, int(math.ceil(2.0 / Lx)))
    ymax = max(np.nanmax(np.abs(np.real(F[c]))) for c in live)
    fig, axes = plt.subplots(len(live), 1, figsize=(13, 1.2 + 2.6 * len(live)), sharex=True, squeeze=False)
    itf = [used["y"][j] * lam * 1e9 for j in range(1, geo.Ny) if geo.eps_cell[j] != geo.eps_cell[j - 1]]
    for ax, c in zip(axes[:, 0], live):
        xs = used["x"][:-1] if c in ("Ey", "Ez") else used["x_dual"]
        ys = used["y"] if c in ("Ex", "Ez") else used["y_dual"]
        data = np.real(F[c][:, :len(ys)])
        X = np.concatenate([xs + t * Lx for t in range(tiles)]) * lam * 1e9
        im = ax.pcolormesh(ys * lam * 1e9, X, np.tile(data, (tiles, 1)), shading="nearest", cmap="RdBu_r",
                           vmin=-ymax, vmax=ymax, rasterized=True)
        for yy in itf:
            ax.axvline(yy, color="k", lw=0.8)
        ax.axvline(used["y"][geo.j0] * lam * 1e9, color="lime", lw=1.0, ls="--")
        for lo, hi in ((0, used["y"][geo.npml_lo]), (used["y"][geo.Ny - geo.npml_hi], used["y"][-1])):
            ax.axvspan(lo * lam * 1e9, hi * lam * 1e9, color="gray", alpha=0.25, lw=0)
        ax.set_ylabel("x [nm]")
        ax.set_title(f"Re {c} (t = 0), DFT phasor; black: interfaces, green: TF/SF plane, gray: PML", fontsize=9,
                     loc="left")
        fig.colorbar(im, ax=ax, pad=0.01)
    axes[-1, 0].set_xlabel("y [nm] (propagation axis)")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def spectrum_plot(cfg, rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    for pol, col in (("s", "#1d4ed8"), ("p", "#b91c1c")):
        rr = sorted([r for r in rows if r["pol"] == pol], key=lambda r: r["wavelength_nm"])
        if not rr:
            continue
        w = [r["wavelength_nm"] for r in rr]
        ax.plot(w, [r["R"] for r in rr], "o", color=col, label=f"R FDTD ({pol})")
        ax.plot(w, [r["T"] for r in rr], "s", color=col, mfc="none", label=f"T FDTD ({pol})")
        if "R_tmm" in rr[0]:                                     # continuous TMM on a dense grid, realized angle
            kt = 2 * math.pi * math.sin(math.radians(rr[0]["info"]["angle_deg"]))
            wd = np.linspace(min(w), max(w), 300)
            ns = float(cfg["substrate_n"])
            ref = [tmm.tmm_continuous([(q["eps"], q["d"]) for q in layer_specs(cfg, x * 1e-9)], 1.0, ns ** 2, kt,
                                      pol=pol) for x in wd]
            ax.plot(wd, [c["R"] for c in ref], "-", color=col, lw=0.8, label=f"R, T TMM ({pol})")
            ax.plot(wd, [c["T"] for c in ref], "-", color=col, lw=0.8)
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("R, T")
    ax.set_title(f"{cfg['name']}: θ = {cfg['angle_deg']}°", fontsize=10)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- driver
def describe(lam, info):
    L = [f"λ = {lam * 1e9:g} nm  mesh = {info['mesh']}  injection = "
         f"{'aux line (exact)' if info['inc'] == 'a' else 'analytic plane wave (x nonuniform)'}",
         f"  angle: requested θ = {info['requested_angle_deg']:g}°, φ = {info['requested_azimuth_deg']:g}°; realized "
         f"θ = {info['angle_deg']:.4f}°, φ = {info['azimuth_deg']:.3f}° (m = {info['m']}, n = {info['n']}, "
         f"Lx = {info['Lx_nm']:.5g} nm, Lz = {info['Lz_nm']:.5g} nm)",
         f"  cells {info['Nx']} x {info['Ny']} x {info['Nz']} = {info['cells']:.3g}; Δy {info['dy_min_nm']:.4g}–"
         f"{info['dy_max_nm']:.4g} nm (max ratio {info['r_max']:.3f}); Δx min {info['dx_min_nm']:.4g} nm, "
         f"Δz {info['dz_nm']:.4g} nm",
         f"  Δt = {info['dt_fs']:.4g} fs ({info['steps_per_period']} steps/period), {info['steps']} steps; "
         f"estimated {info['est_minutes']:.1f} min per run" + (" (+ a vacuum normalization run)"
                                                              if info["grating"] and info["inc"] == "p" else "")]
    for q in info["layers"]:
        if abs(q["realized_nm"] - q["requested_nm"]) > 1e-9:
            L.append(f"  note: layer '{q['name']}' {q['requested_nm']:.4g} nm -> {q['realized_nm']:.4g} nm on the "
                     f"uniform mesh (interfaces must be on nodes; use mesh.type = nonuniform for the exact thickness)")
    if "grating_duty_realized" in info:
        L.append(f"  grating duty realized: {info['grating_duty_realized']:.4f}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", nargs="?", help="JSON file of physical parameters")
    ap.add_argument("--example", choices=sorted(EXAMPLES), help="print a template configuration and exit")
    ap.add_argument("--dry-run", action="store_true", help="build the grids and print the plan; do not run")
    ap.add_argument("--pol", nargs="+", choices=["s", "p"], help="override the polarizations of the file")
    a = ap.parse_args()
    if a.example:
        print(json.dumps(EXAMPLES[a.example], indent=2, ensure_ascii=False))
        return 0
    if not a.config:
        ap.error("a configuration file is required (or --example NAME)")
    cfg = load_config(a.config)
    if a.pol:
        cfg["polarizations"] = a.pol
    outdir = os.path.join(ROOT, cfg["output"]["dir"]) if not os.path.isabs(cfg["output"]["dir"]) else cfg["output"]["dir"]
    rows = []
    for lam in wavelengths(cfg):
        g, geo, dt, info = build(cfg, lam)
        print(describe(lam, info), flush=True)
        if a.dry_run:
            continue
        if not os.path.exists(BIN):
            subprocess.run(["make", "-s"], cwd=ROOT, check=True)
        tag = f"{lam * 1e9:.6g}nm"
        for pol in cfg["polarizations"]:
            fm = bool(cfg["output"]["field_map"])
            if info["grating"]:
                planes = list(grating_planes(geo))
                out = os.path.join(outdir, f"{tag}_{pol}")
                meta = run_solver(out, g, **solver_args(cfg, info, dt, pol, geo, planes, fm))
                out_v = meta_v = None
                if info["inc"] == "p":
                    gv, _, _, _ = build(cfg, lam, vacuum=True)
                    out_v = os.path.join(outdir, f"{tag}_{pol}_vac")
                    meta_v = run_solver(out_v, gv, **solver_args(cfg, info, dt, pol, geo, planes))
                r = measure_grating(cfg, lam, out, meta, out_v, meta_v, geo, info, pol)
                line = (f"  {pol}: ΣR = {r['R']:.6f}  ΣT = {r['T']:.6f}  energy (conserved flux) − 1 = "
                        f"{r['energy_conserved_flux'] - 1:+.1e}")
                if "max_abs_diff_vs_rcwa" in r:
                    line += f"  max |η − η_RCWA| = {r['max_abs_diff_vs_rcwa']:.2e}"
                print(line)
                for key, lab in (("R_orders", "R"), ("T_orders", "T")):
                    ref = r.get("rcwa_" + key, {})
                    print("     " + ", ".join(f"{lab}{p:+d} {v:.5f}" + (f" (RCWA {ref[p]:.5f})" if p in ref else "")
                                              for p, v in r[key].items()))
            else:
                out = os.path.join(outdir, f"{tag}_{pol}")
                meta = run_solver(out, g, **solver_args(cfg, info, dt, pol, geo, None, fm))
                r = measure_planar(out, meta, geo)
                r.update(references_planar(cfg, lam, meta, geo, dt, info, pol))
                print(f"  {pol}: R = {r['R']:.6f}  T = {r['T']:.6f}  R+T−1 = {r['R'] + r['T'] - 1:+.1e} | "
                      f"TMM R = {r['R_tmm']:.6f}  T = {r['T_tmm']:.6f} | exact discrete R = {r['R_disc']:.6f}  "
                      f"(FDTD − discrete {r['R'] - r['R_disc']:+.1e}, discretization {r['R_disc'] - r['R_tmm']:+.1e})")
                if "R_tmm_requested" in r:
                    print(f"     TMM above uses the snapped thicknesses; at the requested thicknesses TMM R = "
                          f"{r['R_tmm_requested']:.6f}  T = {r['T_tmm_requested']:.6f} (FDTD − TMM "
                          f"{r['R'] - r['R_tmm_requested']:+.1e})")
            if abs(r.get("energy_conserved_flux", r["R"] + r["T"]) - 1) > 1e-3:
                print("     warning: energy off by > 1e-3 -- not in steady state? raise run.periods "
                      "(resonant structures ring for a long time)")
            r.update(pol=pol, wavelength_nm=lam * 1e9, runtime_s=meta["runtime_s"], info=info, run_dir=out)
            rows.append(r)
            if fm:
                fig = os.path.join(out, "field_xy.png")
                field_map(out, meta, geo, f"{cfg['name']}: {pol}-pol, λ = {lam * 1e9:g} nm, θ = "
                                          f"{info['angle_deg']:.3f}°, φ = {info['azimuth_deg']:.2f}° (x–y slice, z = 0)",
                          fig, lam)
                print(f"     field map: {os.path.relpath(fig, ROOT)}  (run time {meta['runtime_s']:.0f} s)")
    if a.dry_run:
        return 0
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(config=cfg, results=rows), fh, indent=1, ensure_ascii=False, default=float)
    if len(wavelengths(cfg)) > 1:
        spectrum_plot(cfg, rows, os.path.join(outdir, "spectrum.png"))
        print(f"spectrum: {os.path.relpath(os.path.join(outdir, 'spectrum.png'), ROOT)}")
    print(f"summary: {os.path.relpath(os.path.join(outdir, 'summary.json'), ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
