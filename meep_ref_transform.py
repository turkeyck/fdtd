"""Meep reference for gate 3 (B): the nonuniform y grid as transformation optics on a uniform u grid.

    MEEP_PYTHON meep_ref_transform.py --grid grids/L3B_mapped_film.json --pol s --dt 0.0102040816 \
        --variant cont|disc --out runs/nu/meep/B_film_s_cont

y = f(u), u uniform with Δu = du (= λ0/20, also the x/z spacing), s = dy/du. Equivalent material (Meep, u coordinates):
    eps' = eps·diag(s, 1/s, s),  mu' = diag(s, 1/s, s)          (material_function, eps_averaging = False)
Field map: tangential components unchanged; E'_u = s E_y, H'_u = s H_y.
variant cont: s(u) = 1 − amp·bump(u) sampled at each component's own position (the continuous metric).
variant disc: s at integer u nodes = d_j/Δu and at half nodes = h_j/Δu (the discrete metric of our grid) -- diagnostic
              3B-3: it removes the metric-sampling difference.
Interface nodes (film on primal nodes) get the arithmetic mean eps, as our tangential-E rule does for equal
neighbouring spacings. Meep Courant = dt_own / Δu so both codes step with the same dt.
Sample positions follow Meep's Yee rule (checked against every DFT array shape; y samples must land on our nodes).
Writes meep_B.npz: reduced complex profiles a_c(j) (Bloch projection over x, z) with our node index j.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

try:
    import meep as mp
except ImportError:
    sys.stderr.write("meep is not installed in this interpreter (use MEEP_PYTHON)\n")
    sys.exit(3)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid_gen  # noqa: E402

COMPS = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez, "Hx": mp.Hx, "Hy": mp.Hy, "Hz": mp.Hz}
KX, KZ = math.pi, 2 * math.pi / 3


def basis(pol):
    k0 = 2 * math.pi
    ky = math.sqrt(k0 ** 2 - KX ** 2 - KZ ** 2)
    k = np.array([KX, ky, KZ])
    s = np.cross(k, [0, 1, 0])
    s /= np.linalg.norm(s)
    p = np.cross(s, k)
    p /= np.linalg.norm(p)
    e = s if pol == "s" else p
    return np.array([e[0], 0.0, e[2]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--pol", default="s")
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--variant", default="cont", choices=("cont", "disc"))
    ap.add_argument("--decay", type=float, default=1e-10)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shapes-only", action="store_true", help="check the sample-position rule and exit")
    a = ap.parse_args()
    with open(a.grid) as fh:
        g = json.load(fh)
    gen = g["generator"]
    du = gen["du"]
    pars = dict(du=du, amp=gen["amp"], u14=tuple(gen["u14"]))
    h = np.array(g["y"]["h"])
    Ny = len(h)
    d = np.empty(Ny + 1)
    d[1:-1] = 0.5 * (h[:-1] + h[1:])
    d[0], d[-1] = h[0], h[-1]
    eps_n = np.array(g["eps_y"])
    eps_t = np.empty(Ny + 1)
    eps_t[1:-1] = np.where(eps_n[:-1] == eps_n[1:], eps_n[1:],
                           (h[:-1] * eps_n[:-1] + h[1:] * eps_n[1:]) / (h[:-1] + h[1:]))
    eps_t[0], eps_t[-1] = eps_n[0], eps_n[-1]
    npml, j0 = g["zones"]["npml_lo"], g["zones"]["j0"]
    res = int(round(1.0 / du))
    Lx = Lz = 4 * du
    Ly = Ny * du
    if Ny % 2:
        sys.exit("Ny must be even so that our nodes fall on Meep grid points")

    def u_of(p):
        return p.y + Ly / 2

    def node(u):
        x = u / du
        m = round(2 * x)
        return m, abs(2 * x - m) < 1e-6

    def s_at(u):
        if a.variant == "cont":
            return float(grid_gen.s_of_u(u, pars))
        m, ok = node(u)
        if ok and m % 2 == 0:
            return float(d[min(max(m // 2, 0), Ny)] / du)
        j = min(max(int(math.floor(u / du)), 0), Ny - 1)
        return float(h[j] / du)

    def eps_at(u):
        m, ok = node(u)
        if ok and m % 2 == 0:
            return float(eps_t[min(max(m // 2, 0), Ny)])
        j = min(max(int(math.floor(u / du)), 0), Ny - 1)
        return float(eps_n[j])

    def mat(p):
        u = u_of(p)
        s, e = s_at(u), eps_at(u)
        return mp.Medium(epsilon_diag=mp.Vector3(e * s, e / s, e * s), mu_diag=mp.Vector3(s, 1.0 / s, s))

    kp = mp.Vector3(KX / (2 * math.pi), 0, KZ / (2 * math.pi))
    courant = a.dt / du
    mp.verbosity(0)
    center = mp.Vector3(0, 0, 0)
    # ---- sample positions: Meep's Yee rule for a DFT volume [lo, hi] (verified in meep_ref.py against array shapes):
    # an integer-offset axis has round((hi-lo)/Δ)+1 samples from lo, a half-offset axis one more from lo−Δ/2.
    OFF = {"Ex": (0.5, 0, 0), "Ey": (0, 0.5, 0), "Ez": (0, 0, 0.5), "Hx": (0, 0.5, 0.5), "Hy": (0.5, 0, 0.5),
           "Hz": (0.5, 0.5, 0)}

    def axis(lo, hi, half):
        n = int(round((hi - lo) / du)) + 1
        return (lo - du / 2 + du * np.arange(n + 1)) if half else (lo + du * np.arange(n))

    def axis_strict(lo, hi, half):
        """Meep's Yee DFT samples for a volume whose limits are off the grid: every grid point inside plus the
        nearest one outside each end (floor/ceil); checked against every component's array shape below."""
        off = 0.5 if half else 0.0
        m0 = math.floor(lo / du - off)
        m1 = math.ceil(hi / du - off)
        return du * (np.arange(m0, m1 + 1) + off)
    # y limits a quarter cell off the grid, so that no sample lies on the boundary of the DFT volume
    lo_y = -Ly / 2 + npml * du + du / 4
    size = mp.Vector3(Lx, -2 * lo_y, Lz)
    pos = {}
    for name in COMPS:
        ox, oy, oz = OFF[name]
        pos[name] = (axis(-Lx / 2, Lx / 2, ox == 0.5), axis_strict(lo_y, -lo_y, oy == 0.5), axis(-Lz / 2, Lz / 2, oz == 0.5))
    # ---- the run
    J = basis(a.pol)

    def amp(p):
        return np.exp(2j * math.pi * (kp.x * p.x + kp.z * p.z))
    y_src = j0 * du - Ly / 2
    src = [mp.Source(mp.GaussianSource(1.0, fwidth=0.3), component=c, center=mp.Vector3(0, y_src, 0),
                     size=mp.Vector3(Lx, 0, Lz), amplitude=complex(v), amp_func=amp)
           for c, v in ((mp.Ex, J[0]), (mp.Ez, J[2])) if abs(v) > 1e-14]
    sim = mp.Simulation(cell_size=mp.Vector3(Lx, Ly, Lz), resolution=res, sources=src, material_function=mat,
                        boundary_layers=[mp.PML(npml * du, direction=mp.Y)], k_point=kp, Courant=courant,
                        eps_averaging=False)
    dft = sim.add_dft_fields(list(COMPS.values()), [1.0], center=center, size=size, yee_grid=True,
                             decimation_factor=1)
    comp = mp.Ez if abs(J[2]) >= abs(J[0]) else mp.Ex
    if a.shapes_only:
        sim.run(until=0.05)
        for name, c in COMPS.items():
            arr = sim.get_dft_array(dft, c, 0)
            xs, ys, zs = pos[name]
            print(name, arr.shape, (len(xs), len(ys), len(zs)), "OK" if arr.shape == (len(xs), len(ys), len(zs)) else "MISMATCH")
        return
    t0 = time.time()
    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, comp, mp.Vector3(0, (Ny - npml - 5) * du - Ly / 2, 0),
                                                            a.decay))
    runtime = time.time() - t0
    out, info = {}, {}
    for name, c in COMPS.items():
        arr = sim.get_dft_array(dft, c, 0)
        xs, ys, zs = pos[name]
        if arr.shape != (len(xs), len(ys), len(zs)):
            sys.exit(f"sample-position rule does not match the DFT array for {name}: {arr.shape} vs "
                     f"{(len(xs), len(ys), len(zs))}")
        # Bloch projection: drop the duplicated periodic end sample (x = ±Lx/2 or ±Lz/2 appears twice)
        nx = int(round(Lx / du))
        nz = int(round(Lz / du))
        ix = np.arange(len(xs))[:nx] if xs[0] >= -Lx / 2 - 1e-12 else np.arange(1, nx + 1)
        iz = np.arange(len(zs))[:nz] if zs[0] >= -Lz / 2 - 1e-12 else np.arange(1, nz + 1)
        sub = arr[np.ix_(ix, np.arange(len(ys)), iz)]
        ph = np.exp(-1j * (KX * xs[ix][:, None, None] + KZ * zs[iz][None, None, :]))
        red = (sub * ph).mean(axis=(0, 2))
        yv = ys + Ly / 2
        m = np.round(2 * yv / du).astype(int)
        err = float(np.max(np.abs(2 * yv / du - m)))
        half = bool(np.all(m % 2 == 1))
        full = bool(np.all(m % 2 == 0))
        j = (m - 1) // 2 if half else m // 2
        out[name] = red
        out[name + "_j"] = j
        info[name] = dict(position_error_cells=err / 2, dual=half, primal=full)
    os.makedirs(a.out, exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "meep_B.npz"), **out)
    meta = dict(grid=a.grid, grid_hash=g["hash"], pol=a.pol, variant=a.variant, dt=a.dt, courant=courant, du=du,
                resolution=res, Ly=Ly, meep_version=mp.__version__, runtime_s=runtime, meep_time=sim.meep_time(),
                probe=info, finite=bool(all(np.all(np.isfinite(out[k])) for k in COMPS)))
    with open(os.path.join(a.out, "meep_B.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps(dict(runtime_s=runtime, finite=meta["finite"], probe=info)))


if __name__ == "__main__":
    main()
