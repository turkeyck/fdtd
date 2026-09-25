"""Tensor-product grid generator and checker for fdtd3d_oblique (SPEC_nonuniform.md §13.1).

    python3 grid_gen.py --all            # write every grid used by the validation gates into grids/
    python3 grid_gen.py --check FILE...  # validate grid files (exit 1 on any violation)

A y-grid is described by a list of blocks, bottom (near PML) to top (far PML):
  {"kind": "uniform", "n": cells, "h": spacing, "eps": eps, "tag": ...}        exact uniform cells
  {"kind": "layer", "L": thickness, "eps": eps, "hmax": max spacing}          fixed physical thickness; filled
        uniformly with L/ceil(L/hmax); if a neighbour is finer, a geometric grading (ratio <= r_max) is placed
        inside this layer next to that neighbour and the layer is rescaled (factor <= 1) to its exact thickness
  {"kind": "grade", "eps": eps}                                               geometric transition between the
        end spacings of the neighbouring blocks, equal ratio rho = (h_b/h_a)^(1/(m+1)) <= r_max
Material interfaces therefore always fall on primal nodes. Spacings are stored explicitly (never recovered
by subtracting node coordinates), so a uniform grid has h == Delta bit for bit.
"""
import argparse
import hashlib
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
GRIDS = os.path.join(ROOT, "grids")
FORMAT, VERSION = "fdtd-grid", 1
REL = 1e-14


def hmax_for(eps, ppw):
    return 1.0 / (math.sqrt(eps) * ppw)


# ----------------------------------------------------------------------------- 1D builders
def _geometric(ha, hb, r_max):
    """Intermediate cells strictly between spacings ha and hb with a constant ratio <= r_max."""
    if abs(hb / ha - 1.0) < 1e-12:
        return []
    m1 = math.ceil(abs(math.log(hb / ha)) / math.log(r_max) - 1e-12)
    rho = (hb / ha) ** (1.0 / m1)
    return [ha * rho ** q for q in range(1, m1)]


def _layer_cells(L, hmax, h_left, h_right, r_max):
    """Fill thickness L with spacings <= hmax, grading geometrically away from finer neighbours; exact total L.

    For n cells the profile is c_i = min(cap, hl*rho^(i+1), hr*rho^(n-i)) with a common ratio rho in [1, r_max];
    the smallest n for which the profile at rho = r_max (cap = L/ceil(L/hmax)) reaches L is used, and rho (then
    the cap) is found by bisection so that sum(c) = L exactly. Neighbour ratios are then <= rho <= r_max."""
    n0 = math.ceil(L / hmax - 1e-12)
    hu = L / n0
    hl = h_left if (h_left is not None and h_left < hu * (1 - 1e-12)) else math.inf
    hr = h_right if (h_right is not None and h_right < hu * (1 - 1e-12)) else math.inf
    if hl == math.inf and hr == math.inf:
        return [hu] * n0

    def prof(n, rho, cap):
        i = np.arange(n)
        return np.minimum(np.minimum(cap, hl * rho ** (i + 1.0)), hr * rho ** (n - i + 0.0))

    n = n0
    while prof(n, r_max, hu).sum() < L:
        n += 1
    lo_h = min(hl, hr)
    if n * min(lo_h, hu) > L * (1 + 1e-12):
        raise ValueError(f"layer L={L} too thin for neighbour spacing {lo_h} with r_max={r_max}")
    a, b = 1.0, r_max                                  # sum is increasing in rho
    for _ in range(200):
        m = 0.5 * (a + b)
        if prof(n, m, hu).sum() < L:
            a = m
        else:
            b = m
    cells = prof(n, b, hu)
    return list(cells * (L / cells.sum()))


def build_y(blocks, r_max):
    """Return (h, eps, marks): cell spacings, cell eps and {tag: (first_cell, n_cells)}."""
    # pass 1: nominal uniform spacing of every non-grade block
    nominal = []
    for b in blocks:
        if b["kind"] == "uniform":
            nominal.append(b["h"])
        elif b["kind"] == "layer":
            nominal.append(b["L"] / math.ceil(b["L"] / b["hmax"] - 1e-12))
        else:
            nominal.append(None)

    def end_h(i, side):
        """Spacing that block i presents at its 'left'/'right' end (layers: their nominal uniform spacing)."""
        return nominal[i]

    h, eps, marks = [], [], {}
    for i, b in enumerate(blocks):
        if b["kind"] == "uniform":
            cells = [b["h"]] * b["n"]
        elif b["kind"] == "layer":
            hl = h[-1] if h else None
            hr = end_h(i + 1, "left") if i + 1 < len(blocks) else None   # None next to a grade block
            cells = _layer_cells(b["L"], b["hmax"], hl, hr, r_max)
        else:
            ha = h[-1]
            hb = end_h(i + 1, "left")
            cells = _geometric(ha, hb, r_max)
        if "tag" in b:
            marks[b["tag"]] = (len(h), len(cells))
        h += cells
        eps += [b["eps"]] * len(cells)
    return np.array(h), np.array(eps, float), marks


def uniform_axis(L, delta):
    n = int(round(L / delta))
    if abs(n * delta - L) > 1e-12 * L:
        raise ValueError(f"L={L} is not a multiple of {delta}")
    return {"uniform": True, "h": [delta] * n, "L": L}


# ----------------------------------------------------------------------------- grid dict
def canonical_hash(g):
    core = {k: g[k] for k in ("x", "y", "z", "eps_y", "zones", "grating") if k in g}
    for k in ("eps_t", "eps_n"):
        if k in g:
            core[k] = g[k]
    s = json.dumps(core, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(s.encode()).hexdigest()


def make_grid(name, blocks, x, z, r_max, npml_lo, npml_hi, j0, na=40, ppw=None, mode="blocks", extra=None,
              eps_t=None, eps_n=None):
    h, eps, marks = build_y(blocks, r_max) if blocks is not None else (None, None, {})
    return finish(name, h, eps, x, z, r_max, npml_lo, npml_hi, j0, na, ppw, mode, extra, marks, eps_t, eps_n)


def finish(name, h, eps, x, z, r_max, npml_lo, npml_hi, j0, na, ppw, mode, extra, marks, eps_t=None, eps_n=None):
    h = [float(v) for v in h]
    eps = [float(v) for v in eps]
    interfaces = [j for j in range(1, len(eps)) if eps[j] != eps[j - 1]]
    g = {
        "format": FORMAT, "version": VERSION, "units": "lambda0", "name": name,
        "generator": {"mode": mode, "ppw": ppw, "r_max": r_max, **(extra or {})},
        "x": x, "z": z,
        "y": {"h": h, "y0": 0.0},
        "eps_y": eps,
        "interfaces_j": interfaces,
        "zones": {"npml_lo": npml_lo, "npml_hi": npml_hi, "j0": j0, "ja": j0 - na, "uniform_min_cells": 5},
        "marks": {k: list(v) for k, v in marks.items()},
    }
    if eps_t is not None:
        g["eps_t"] = [float(v) for v in eps_t]
    if eps_n is not None:
        g["eps_n"] = [float(v) for v in eps_n]
    g["hash"] = canonical_hash(g)
    return g


def periodic_axis_layers(layers, r_max=1.1):
    """Periodic axis from layers [(L, hmax), ...]; each layer grades from finer neighbours (wrapping around).
    Returns (axis dict, first cell index of every layer)."""
    nom = [L / math.ceil(L / hm - 1e-12) for L, hm in layers]
    cells, starts = [], []
    n = len(layers)
    for q, (L, hm) in enumerate(layers):
        hl, hr = nom[(q - 1) % n], nom[(q + 1) % n]
        starts.append(len(cells))
        cells += _layer_cells(L, hm, hl, hr, r_max)
    Ltot = sum(L for L, _ in layers)
    h = np.array(cells)
    uni = bool(np.all(np.abs(h / h[0] - 1) < REL))
    return {"uniform": uni, "h": [float(v) for v in cells], "L": Ltot}, starts


def stageB_family(k=1, film=True, npml=30, name=None):
    """Amendment D16: stage-B injection/order families. Classical incidence (m = 1, n = 0, Lx = 2): x graded
    (layer [0,1) at λ0/(28.6k), layer [1,2) at λ0/(20k)), z thin uniform (Lz = 0.1, 4k cells), y uniform aligned
    h = 0.02/k: PML | SF 0.5 | TF 1.0 | film 0.08 (n = 2) | substrate 1.0 (n = 1.46) | PML (vacuum if film=False)."""
    x, _ = periodic_axis_layers([(1.0, 0.035 / k), (1.0, 0.05 / k)])
    z = uniform_axis(0.1, 0.025 / k)
    h = 0.02 / k
    nf = int(round(D_FILM / h))
    blocks = [dict(kind="uniform", n=npml * k, h=h, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=int(round(0.5 / h)), h=h, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(1.0 / h)), h=h, eps=1.0, tag="tf_vac"),
              dict(kind="uniform", n=nf, h=h, eps=N_FILM ** 2 if film else 1.0, tag="film"),
              dict(kind="uniform", n=int(round(1.0 / h)), h=h, eps=N_SUB ** 2 if film else 1.0, tag="sub"),
              dict(kind="uniform", n=npml * k, h=h, eps=N_SUB ** 2 if film else 1.0, tag="pml_hi")]
    j0 = npml * k + int(round(0.5 / h))
    nm = name or f"B_{'film' if film else 'vac'}_k{k}"
    return make_grid(nm, blocks, x, z, 1.1, npml * k, npml * k, j0, ppw=None, mode="stageB", extra=dict(k=k, m=1, n=0))


def stageB_ops_grid():
    """B0-2: x AND z graded (periodic), y = the gate-1 taper grid with PEC walls (source-free operator tests)."""
    x, _ = periodic_axis_layers([(0.5, 0.035), (0.5, 0.05)])
    z, _ = periodic_axis_layers([(0.6, 0.03), (0.6, 0.05)])
    g = level1(20, 1.1)
    return finish("B_ops_xz_graded_pec", g["y"]["h"], g["eps_y"], x, z, 1.1, 0, 0, g["zones"]["j0"], 40, None,
                  "stageB_ops", dict(parent=g["name"]), {})


def stageB_grating(ppw=40, npml=20, name=None):
    """B3: lamellar grating along x (period Λ = Lx = 2), ridge n = 2 on x ∈ [0, 1) (x nodes on both edges),
    groove vacuum, thickness 0.3, conical incidence (A2: m = n = 1, Lz = 3). x: ridge λ0/(2 ppw), groove
    λ0/ppw graded (r_max 1.1); z uniform λ0/ppw; y as the gate-2 layout (vacuum | grating | substrate)."""
    x, starts = periodic_axis_layers([(1.0, hmax_for(N_FILM ** 2, ppw)), (1.0, hmax_for(1.0, ppw))])
    z = uniform_axis(3.0, 1.0 / ppw)
    hv = hmax_for(1.0, ppw)
    hs = hmax_for(N_SUB ** 2, ppw)
    blocks = [dict(kind="uniform", n=npml, h=hv, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=int(round(1.0 / hv)), h=hv, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(2.0 / hv)), h=hv, eps=1.0, tag="tf_vac"),
              dict(kind="grade", eps=1.0, tag="grade_vac"),
              dict(kind="layer", L=0.3, eps=1.0, hmax=hmax_for(N_FILM ** 2, ppw), tag="grating"),
              dict(kind="grade", eps=N_SUB ** 2, tag="grade_sub"),
              dict(kind="uniform", n=int(round(2.0 / hs)), h=hs, eps=N_SUB ** 2, tag="sub"),
              dict(kind="uniform", n=npml, h=hs, eps=N_SUB ** 2, tag="pml_hi")]
    j0 = npml + int(round(1.0 / hv))
    h, eps, marks = build_y(blocks, 1.1)
    c0, n = marks["grating"]
    g = finish(name or f"B_grating_ppw{ppw}", h, eps, x, z, 1.1, npml, npml, j0, 40, ppw, "grating",
               dict(period=2.0, duty=0.5, ridge_x=[0.0, 1.0]), marks)
    g["grating"] = {"j_lo": c0, "j_hi": c0 + n, "i_lo": starts[0], "i_hi": starts[1], "eps_ridge": N_FILM ** 2,
                    "eps_groove": 1.0, "thickness": 0.3}
    g["hash"] = canonical_hash(g)
    return g


def stageB_grating_vac(ppw=40):
    """Normalization grid for B3: identical spacings, vacuum everywhere, no grating."""
    g = stageB_grating(ppw)
    z = g["zones"]
    return finish(f"B_grating_ppw{ppw}_vac", g["y"]["h"], [1.0] * len(g["eps_y"]), g["x"], g["z"], 1.1, z["npml_lo"],
                  z["npml_hi"], z["j0"], z["j0"] - z["ja"], ppw, "grating_norm", dict(parent=g["name"]),
                  {k: tuple(v) for k, v in g["marks"].items()})


def nodes(h, y0=0.0):
    """Primal node coordinates from spacings (compensated summation)."""
    y = [y0]
    for q in range(len(h)):
        y.append(y0 + math.fsum(h[: q + 1]))
    return np.array(y)


# ----------------------------------------------------------------------------- checker
def check(g):
    """Return a list of violations (empty = usable)."""
    bad = []
    if g.get("format") != FORMAT or g.get("version") != VERSION:
        bad.append("format/version")
    if canonical_hash(g) != g.get("hash"):
        bad.append("hash mismatch (file edited by hand?)")
    r_max = g["generator"]["r_max"]
    for ax in ("x", "z"):
        a = g[ax]
        hx = np.array(a["h"])
        if abs(hx.sum() - a["L"]) > 1e-13 * max(a["L"], 1):
            bad.append(f"{ax}: sum(h) != L")
        uni = bool(np.all(np.abs(hx / hx[0] - 1) < REL))
        if uni != bool(a["uniform"]):
            bad.append(f"{ax}: 'uniform' flag inconsistent")
        rr = np.concatenate([hx[1:] / hx[:-1], [hx[0] / hx[-1]]])
        if np.any(rr > r_max * (1 + 1e-12)) or np.any(rr < 1 / (r_max * (1 + 1e-12))):
            bad.append(f"{ax}: adjacent ratio exceeds r_max (incl. periodic seam)")
    h = np.array(g["y"]["h"])
    e = np.array(g["eps_y"])
    Ny = len(h)
    if len(e) != Ny:
        bad.append("len(eps_y) != len(y.h)")
    if np.any(h <= 0):
        bad.append("non-positive spacing")
    rr = h[1:] / h[:-1]
    if np.any(rr > r_max * (1 + 1e-12)) or np.any(rr < 1 / (r_max * (1 + 1e-12))):
        k = int(np.argmax(np.maximum(rr, 1 / rr)))
        bad.append(f"y: adjacent ratio {max(rr[k], 1 / rr[k]):.4f} > r_max {r_max} at cell {k}")
    itf = [j for j in range(1, Ny) if e[j] != e[j - 1]]
    if itf != list(g["interfaces_j"]):
        bad.append("interfaces_j inconsistent with eps_y")
    z = g["zones"]
    u = z["uniform_min_cells"]

    def uniform(lo, hi, what, ref=None):
        lo, hi = max(lo, 0), min(hi, Ny)
        if hi <= lo:
            return
        seg = h[lo:hi]
        r0 = seg[0] if ref is None else ref
        if np.any(np.abs(seg / r0 - 1) > REL) or np.any(e[lo:hi] != e[lo]):
            bad.append(f"{what}: cells [{lo},{hi}) not uniform (max rel diff {np.max(np.abs(seg / r0 - 1)):.1e})")
    uniform(0, z["npml_lo"] + u, "near PML + margin")
    uniform(Ny - z["npml_hi"] - u, Ny, "far PML + margin")
    uniform(z["j0"] - u, z["j0"] + u, "TF/SF plane j0 +-5")
    ja = z["ja"]
    if ja - u < 0:
        uniform(0, ja + u, "aux source (extension below node 0)", ref=h[0])
    else:
        uniform(ja - u, ja + u, "aux source j_a +-5")
    ppw = g["generator"].get("ppw")
    if ppw and g["generator"]["mode"] != "mapped":
        hm = 1.0 / (np.sqrt(e) * ppw)
        if np.any(h > hm * (1 + 1e-12)):
            k = int(np.argmax(h / hm))
            bad.append(f"cell {k}: h={h[k]:.6g} > lambda0/(n ppw)={hm[k]:.6g}")
    for key, n_exp in (("eps_t", Ny + 1), ("eps_n", Ny)):
        if key in g and len(g[key]) != n_exp:
            bad.append(f"{key}: length {len(g[key])} != {n_exp}")
    return bad


# ----------------------------------------------------------------------------- the validation grids
DXZ = 0.05                                 # x/z spacing lambda0/20 (Stage A)
AX_X, AX_Z = (lambda d=DXZ: uniform_axis(2.0, d)), (lambda d=DXZ: uniform_axis(3.0, d))


def level1(base=20, r=1.1, npml_far=60, hold=0.5, name=None, ramp=True, dxz=DXZ):
    """near PML 20 | SF 20 | j0 | uniform 2 lambda0 | taper to base/4 | hold | taper back | uniform 3 | far PML."""
    h0 = 1.0 / base
    hf = h0 / 4
    blocks = [dict(kind="uniform", n=20, h=h0, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=20, h=h0, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(2.0 / h0)), h=h0, eps=1.0, tag="tf_a")]
    if ramp:
        blocks += [dict(kind="grade", eps=1.0, tag="down"),
                   dict(kind="uniform", n=int(round(hold / hf)), h=hf, eps=1.0, tag="hold"),
                   dict(kind="grade", eps=1.0, tag="up")]
    blocks += [dict(kind="uniform", n=int(round(3.0 / h0)), h=h0, eps=1.0, tag="tf_b"),
               dict(kind="uniform", n=npml_far, h=h0, eps=1.0, tag="pml_hi")]
    name = name or f"L1_taper_r{r:g}_base{base}"
    return make_grid(name, blocks, uniform_axis(2.0, dxz), uniform_axis(3.0, dxz), max(r, 1.0 + 1e-9) if ramp else 1.1,
                     20, npml_far, 40, ppw=base, extra=dict(base=base, r=r, hold=hold, dxz=dxz))


def pec_variant(g, name, Lx=None, Lz=None):
    """Same y grid with no PML (PEC walls at j = 0, Ny) for the source-free gate-0 tests; optionally a thin
    x/z cell (4 x 6 cells of the same spacing) -- the y operator, which carries the nonuniformity, is unchanged."""
    h = g["y"]["h"]
    x = uniform_axis(Lx, g["x"]["h"][0]) if Lx else g["x"]
    z = uniform_axis(Lz, g["z"]["h"][0]) if Lz else g["z"]
    return finish(name, h, g["eps_y"], x, z, g["generator"]["r_max"], 0, 0, g["zones"]["j0"],
                  g["zones"]["j0"] - g["zones"]["ja"], g["generator"]["ppw"], "pec", dict(parent=g["name"]), {})


def uniform_parity(nl=20, npml=20, sf=20, tf=200, eps2=None, y1=None, name=None):
    """Equally spaced grid identical to the legacy mesh=uniform layout (for gate 0-1a)."""
    D = 1.0 / nl  # lint-ok: uniform generator (legacy layout)
    Ny = 2 * npml + sf + tf
    eps = [1.0] * Ny
    j0 = npml + sf
    if eps2 is not None:
        for j in range(j0 + y1, Ny):
            eps[j] = eps2
    return finish(name or f"P_uniform_nl{nl}", [D] * Ny, eps, uniform_axis(2.0, D), uniform_axis(3.0, D), 1.1,
                  npml, npml, j0, 40, None, "uniform", dict(nl=nl, sf=sf, tf=tf), {})


N_FILM, D_FILM = 2.0, 0.08
N_SUB = 1.46
STACK5 = [(2.0, 0.08), (1.46, 0.12), (2.0, 0.08), (1.46, 0.12), (2.0, 0.08)]


def level2(struct="F1", ppw=40, dxz=DXZ, npml=20, r_max=1.1, name=None):
    """vacuum (near PML | SF 1 | TF 3 lambda0) | grade | layers | substrate 3 lambda0 | far PML (substrate)."""
    hv = hmax_for(1.0, ppw)
    hs = hmax_for(N_SUB ** 2, ppw)
    layers = [(N_FILM, D_FILM)] if struct == "F1" else STACK5
    blocks = [dict(kind="uniform", n=npml, h=hv, eps=1.0, tag="pml_lo"),
              dict(kind="uniform", n=int(round(1.0 / hv)), h=hv, eps=1.0, tag="sf"),
              dict(kind="uniform", n=int(round(3.0 / hv)), h=hv, eps=1.0, tag="tf_vac")]
    lay = [dict(kind="layer", L=d, eps=n * n, hmax=hmax_for(n * n, ppw), tag=f"layer{q}")
           for q, (n, d) in enumerate(layers)]
    blocks += [dict(kind="grade", eps=1.0, tag="grade_vac")] + lay
    blocks += [dict(kind="grade", eps=N_SUB ** 2, tag="grade_sub"),
               dict(kind="uniform", n=int(round(3.0 / hs)), h=hs, eps=N_SUB ** 2, tag="sub"),
               dict(kind="uniform", n=npml, h=hs, eps=N_SUB ** 2, tag="pml_hi")]
    j0 = npml + int(round(1.0 / hv))
    name = name or f"L2_{struct}_ppw{ppw}" + ("" if abs(dxz - DXZ) < 1e-15 else f"_dxz{int(round(1 / dxz))}")
    return make_grid(name, blocks, uniform_axis(2.0, dxz), uniform_axis(3.0, dxz), r_max, npml, npml, j0,
                     ppw=ppw, extra=dict(struct=struct, dxz=dxz))


def level2_unaligned(ppw=20, offset=0.3, npml=20, name=None):
    """Control 2-4: uniform grid at the vacuum spacing, film interfaces NOT on nodes (start = node + offset*h),
    point-sampled eps (tangential E at nodes, Ey at half nodes) = staircase."""
    hv = hmax_for(N_FILM ** 2, ppw)                  # uniform spacing that also resolves the film
    n_vac = npml + int(round(1.0 / hv)) + int(round(3.0 / hv))
    n_tot = n_vac + int(round(0.5 / hv)) + int(round(3.0 / hv)) + npml
    y = np.arange(n_tot + 1) * hv
    y_a = (n_vac + offset) * hv
    y_b = y_a + D_FILM

    def eps_at(yy):
        return np.where(yy < y_a, 1.0, np.where(yy < y_b, N_FILM ** 2, N_SUB ** 2))
    eps_t = eps_at(y)
    eps_n = eps_at(0.5 * (y[1:] + y[:-1]))
    j0 = npml + int(round(1.0 / hv))
    return finish(name or f"L2_F1_unaligned_ppw{ppw}", [hv] * n_tot, list(eps_n), AX_X(), AX_Z(), 1.1, npml, npml,
                  j0, 40, ppw, "unaligned", dict(offset=offset), {}, eps_t=eps_t, eps_n=eps_n)


def level2_bisect(struct="F1", k=1, npml=20):
    """Amendment D10: nested bisection of the PPW=20 grid (every cell split into k equal cells), x/z fixed;
    the gate runs these with dt = dt(ppw20 grid)/k so that dt scales with h."""
    g0 = level2(struct, 20, npml=npml)
    h = np.repeat(np.array(g0["y"]["h"]) / k, k)
    e = np.repeat(np.array(g0["eps_y"]), k)
    z = g0["zones"]
    return finish(f"L2_{struct}_bis{k}", h, e, g0["x"], g0["z"], 1.1, npml * k, npml * k, z["j0"] * k, 40 * k,
                  20 * k, "bisect", dict(struct=struct, k=k, parent=g0["name"], dt_rule="dt(parent)/k"), {})


def level2_bisect_xz(struct="F1", k=1, npml=20):
    """Amendment D14 (gate 3A-3): bisection of the PPW=20 grid in y AND x/z spacing λ0/(20k): every spacing
    halves with k (dt = dt(bis1)/k), so the 2-level Richardson limit approximates the continuous answer."""
    g = level2_bisect(struct, k, npml)
    g2 = finish(f"L2_{struct}_bisxz{k}", g["y"]["h"], g["eps_y"], uniform_axis(2.0, DXZ / k),
                uniform_axis(3.0, DXZ / k), 1.1, g["zones"]["npml_lo"], g["zones"]["npml_hi"], g["zones"]["j0"],
                g["zones"]["j0"] - g["zones"]["ja"], g["generator"]["ppw"], "bisect_xz",
                dict(struct=struct, k=k, dt_rule="dt(L2_*_bis1)/k"), {})
    return g2


def level2_control(mh=4.5, npml=20, name=None):
    """Amendment D11 (control 2-4): uniform grid, Delta = 0.08/mh (mh = M + 1/2), film interfaces NOT on nodes
    but at fixed cell fractions 0.25 (start) and 0.75 (end) at every resolution; eps point-sampled (staircase):
    eps_t = eps(y_j), eps_n = eps(y_{j+1/2}). The film thickness error is then O(Delta) with a fixed coefficient."""
    D = D_FILM / mh
    n_pml = npml
    n_sf = int(round(1.0 / D))
    n_tf = int(round(3.0 / D))
    n_a = n_pml + n_sf + n_tf
    n_tot = n_a + int(round(3.1 / D)) + n_pml
    y = np.arange(n_tot + 1) * D  # lint-ok: this control grid is uniform by definition
    y_a = (n_a + 0.25) * D  # lint-ok
    y_b = y_a + D_FILM

    def eps_at(yy):
        return np.where(yy < y_a, 1.0, np.where(yy < y_b, N_FILM ** 2, N_SUB ** 2))
    eps_t = eps_at(y)
    eps_n = eps_at(0.5 * (y[1:] + y[:-1]))
    return finish(name or f"L2_F1_ctrl_mh{mh:g}", [D] * n_tot, list(eps_n), AX_X(), AX_Z(), 1.1, n_pml, n_pml,
                  n_pml + n_sf, 40, None, "control", dict(mh=mh, fractions=[0.25, 0.75]), {}, eps_t=eps_t,
                  eps_n=eps_n)


# ---- mapped grid y = f(u) for the Meep transformation-optics comparison (gate 3B)
def _smooth_step(x):
    x = np.clip(x, 0.0, 1.0)
    a = np.where(x > 0, np.exp(-1.0 / np.where(x > 0, x, 1.0)), 0.0)
    b = np.where(x < 1, np.exp(-1.0 / np.where(x < 1, 1.0 - x, 1.0)), 0.0)
    return a / (a + b)


def bump(u, u1, u2, u3, u4):
    """C-infinity bump: 0 for u<u1 or u>u4, 1 on [u2,u3]."""
    return _smooth_step((u - u1) / (u2 - u1)) * (1.0 - _smooth_step((u - u3) / (u4 - u3)))


def s_of_u(u, pars):
    return 1.0 - pars["amp"] * bump(np.asarray(u, float), *pars["u14"])


def mapped_h(un, pars, npts=64):
    """Cell spacings h_j = int_{u_j}^{u_j+1} s(v) dv (Gauss-Legendre per cell; exactly du where s = 1)."""
    xg, wg = np.polynomial.legendre.leggauss(npts)
    h = []
    for a, b in zip(un[:-1], un[1:]):
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        s = s_of_u(mid + half * xg, pars)
        h.append(pars["du"] if np.all(s == 1.0) else half * float(np.dot(wg, s)))
    return np.array(h)


def level3_mapped(film_cells=6, npml=20, name=None):
    """u-grid uniform du = 0.05 (lambda0/20); s(u) = 1 - 0.75 bump, s = 1/4 around the film.
    Layout in u (cells): near PML 20 | SF 20 | j0 | vacuum ... | film (6 cells in the s=1/4 zone) | substrate | PML."""
    du = 0.05
    n_lo = npml + 20 + 40                           # PML, SF, 2 lambda0 uniform after j0
    n_ramp, n_flat = 40, 20                         # bump rise over 40 cells (2 lambda0), flat s=1/4 over 20 + film + 20
    u1 = n_lo * du
    u2 = u1 + n_ramp * du
    u3 = u2 + (n_flat + film_cells + n_flat) * du
    u4 = u3 + n_ramp * du
    n_total = int(round(u4 / du)) + 60 + npml       # 3 lambda0 after the bump, far PML
    pars = dict(du=du, amp=0.75, u14=(u1, u2, u3, u4))
    un = np.arange(n_total + 1) * du
    h = mapped_h(un, pars)
    j_fa = int(round(u2 / du)) + n_flat
    j_fb = j_fa + film_cells
    return j_fa, j_fb, h, pars, n_total


def level3_mapped_grids():
    out = []
    for case in ("vac", "film"):
        j_fa, j_fb, h, pars, n_total = level3_mapped()
        eps = np.ones(len(h))
        if case == "film":
            eps[j_fa:j_fb] = N_FILM ** 2
            eps[j_fb:] = N_SUB ** 2
        g = finish(f"L3B_mapped_{case}", h, eps, AX_X(), AX_Z(), 1.1, 20, 20, 40, 40, None, "mapped",
                   dict(du=pars["du"], amp=pars["amp"], u14=list(pars["u14"]), film_nodes=[j_fa, j_fb]), {})
        out.append(g)
    return out


def write(g):
    os.makedirs(GRIDS, exist_ok=True)
    path = os.path.join(GRIDS, g["name"] + ".json")
    with open(path, "w") as fh:
        json.dump(g, fh, indent=1)
    return path


def all_grids():
    gs = []
    for base in (20, 40, 80):
        gs.append(level1(base, 1.1))
        gs.append(level1(base, 1.1, npml_far=20, name=f"L1_taper_r1.1_base{base}_N20"))
    gs.append(level1(20, 4.0, name="L1_abrupt_r4_base20"))
    gs.append(level1(20, 4.0, npml_far=20, name="L1_abrupt_r4_base20_N20"))
    gs.append(level1(20, ramp=False, name="L1_uniform_ctrl"))
    for N in (10, 20, 30, 60):
        gs.append(level1(20, 1.1, npml_far=N, name=f"L1_taper_r1.1_base20_N{N}"))
    for r in (1.05, 1.2, 1.5):
        for base in (20, 40, 80):
            gs.append(level1(base, r, name=f"L1_scan_r{r:g}_base{base}"))
    for base in (40, 80):
        gs.append(level1(base, 4.0, name=f"L1_scan_r4_base{base}"))
    for src in (level1(20, 1.1), level1(20, 4.0, name="L1_abrupt_r4_base20")):     # gate 0 (PEC, no source)
        tag = src["name"]
        gs.append(pec_variant(src, tag + "_pec"))
        gs.append(pec_variant(src, tag + "_pec_thin", Lx=0.2, Lz=0.3))
    for base in (10, 20, 40):                           # amendment D9 (rev. 2026-09-25): x/z refined with the base
        gs.append(level1(base, 1.1, dxz=1.0 / base, name=f"L1_taper_r1.1_base{base}_xz"))
    for struct in ("F1", "F5"):
        for ppw in (20, 40, 80, 160):
            gs.append(level2(struct, ppw))
        for dxz in (0.025, 0.0125):
            gs.append(level2(struct, 40, dxz=dxz))
    for struct in ("F1", "F5"):                         # amendment D10
        for k in (1, 2, 4, 8):
            gs.append(level2_bisect(struct, k))
    for struct in ("F1", "F5"):                         # amendment D14
        for k in (1, 2):
            gs.append(level2_bisect_xz(struct, k))
    for k in (1, 2, 4):                                 # amendment D16 (stage B families)
        gs.append(stageB_family(k, film=True))
        gs.append(stageB_family(k, film=False))
    gs.append(stageB_ops_grid())
    gs.append(stageB_grating(40))
    gs.append(stageB_grating_vac(40))
    for mh in (4.5, 9.5, 18.5, 37.5):                   # amendment D11
        gs.append(level2_control(mh))
    gs += level3_mapped_grids()
    gs.append(uniform_parity(20))
    gs.append(uniform_parity(20, eps2=2.25, y1=60, name="P_uniform_nl20_med"))
    gs.append(uniform_parity(10, sf=10, tf=80, name="P_uniform_nl10"))
    return gs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--check", nargs="*")
    a = ap.parse_args()
    rc = 0
    if a.all:
        for g in all_grids():
            p = write(g)
            bad = check(g)
            Ny = len(g["y"]["h"])
            hy = np.array(g["y"]["h"])
            print(f"{os.path.basename(p):<34} Ny={Ny:<5} hmin={hy.min():.5f} hmax={hy.max():.5f} "
                  f"{'OK' if not bad else 'VIOLATIONS: ' + '; '.join(bad)}")
            rc |= bool(bad)
    for f in a.check or []:
        with open(f) as fh:
            g = json.load(fh)
        bad = check(g)
        print(f, "OK" if not bad else bad)
        rc |= bool(bad)
    return rc


if __name__ == "__main__":
    sys.exit(main())
