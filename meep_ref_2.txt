"""Meep reference simulations for Level 3 (run with the Meep conda env, e.g. ~/micromamba/envs/mp/bin/python).

Same cell as fdtd3d_oblique: Lx x Ly x Lz with Ly = (2*npml + sf + tf)*Delta, resolution = 1/Delta,
Courant = 0.5, PML only along y (mp.PML(d, direction=mp.Y)), Bloch k_point = (kx/2π, 0, kz/2π).

Meep conventions used here (checked, not assumed):
  * k_point is in units of 2π/a and the Bloch phase is exp(2πi k·r) -- Meep manual, Python User Interface,
    "k_point"; https://meep.readthedocs.io/en/latest/Python_User_Interface/ . The measured phase gradient of
    the DFT fields (kx, kz) is reported by compare.py as the empirical check.
  * Origin at the cell centre; grid points at multiples of Delta. With yee_grid=True each component lives on
    its own Yee sub-grid; for a volume [lo, hi] an integer-offset axis has round((hi-lo)/Δ)+1 samples starting
    at lo and a half-offset axis one more sample starting at lo-Δ/2 (verified against array shapes).
  * DFT decimation is disabled (decimation_factor=1); DFT objects are added only after the ramp so that the
    window [t0, t1] (integer periods) matches the C code's window.

Source: a planar current sheet at our TF/SF plane y0 with J = tangential part of E0 (s: ŝ; p: p̂_t ∥ k̂_t),
amp_func = exp(2πi k·r), time signal g(t) exp(-2πi f t) with the same erf ramp as the C code
(t0 = 20 T0, τ = 4 T0). A current sheet radiates to both ±y; only the +y side is used. Its in-plane
∇·J ≠ 0 (p-pol) leaves a static charge ∝ spectrum of g at ω = 0 relative to ω0, ~exp(-(ω0 τ)²/4) ≈ e^-158:
negligible with this ramp (it would not be with a sharp turn-on).

Outputs (npz in --out): DFT fields on a thin z-slab (x–y slice, same z as the C code's xy slice), on y-planes,
fluxes for R/T, and the metadata needed by compare.py.
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
    sys.stderr.write("meep is not installed in this interpreter; run with the Meep environment "
                     "(e.g. ~/micromamba/envs/mp/bin/python meep_ref.py ...)\n")
    sys.exit(3)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdtd_theory import Setup  # noqa: E402

COMPS = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez, "Hx": mp.Hx, "Hy": mp.Hy, "Hz": mp.Hz}
# Yee offsets (x, y, z) in cells -- identical convention to derivation.md §1
OFF = {"Ex": (0.5, 0, 0), "Ey": (0, 0.5, 0), "Ez": (0, 0, 0.5), "Hx": (0, 0.5, 0.5), "Hy": (0.5, 0, 0.5),
       "Hz": (0.5, 0.5, 0)}


def axis_coords(lo, hi, D, half):
    n = int(round((hi - lo) / D)) + 1
    return (lo - D / 2 + D * np.arange(n + 1)) if half else (lo + D * np.arange(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nl", type=int, default=20)
    ap.add_argument("--pol", default="s")
    ap.add_argument("--medium", action="store_true", help="n=1.5 half-space (Level 2 geometry)")
    ap.add_argument("--eps_averaging", type=int, default=1)
    ap.add_argument("--tf_l", type=int, default=10, help="TF length in lambda0 (10: Level 1 geometry, 8: Level 2)")
    ap.add_argument("--norm", default="", help="vacuum run directory whose incident R-plane data to subtract")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    nl, D = a.nl, 1.0 / a.nl
    npml, sf = 20, nl
    tf = a.tf_l * nl
    if a.medium and a.tf_l != 8:
        sys.exit("medium runs use the Level 2 geometry: --tf_l 8")
    Ny = 2 * npml + sf + tf
    Lx, Lz, Ly = 2.0, 3.0, Ny * D
    j0 = npml + sf
    j1 = j0 + 3 * nl
    st = Setup(n_lambda=nl, pol=a.pol)
    e, h = st.amplitudes()
    per = 1.0  # T0 = lambda0/c = 1
    t_on, t_win0, t_win1 = 20.0, 80.0, 110.0  # erf centre; DFT window in periods (same as the C runs)
    tau = 4.0

    # coordinate maps: ours = meep + L/2
    def to_meep_y(y):
        return y - Ly / 2

    J = np.array([e[0], 0.0, e[2]])  # tangential part of E0
    k_meep = mp.Vector3(st.kx / (2 * math.pi), 0, st.kz / (2 * math.pi))

    def amp(p):
        return np.exp(2j * math.pi * (k_meep.x * p.x + k_meep.z * p.z))

    def g(t):
        return 0.5 * math.erfc(-(t - t_on) / tau) * np.exp(-2j * math.pi * t)

    srcs = []
    for comp, val in ((mp.Ex, J[0]), (mp.Ez, J[2])):
        if abs(val) > 1e-14:
            srcs.append(mp.Source(mp.CustomSource(g, center_frequency=1.0, fwidth=0.5), component=comp,
                                  center=mp.Vector3(0, to_meep_y(j0 * D), 0), size=mp.Vector3(Lx, 0, Lz),
                                  amplitude=complex(val), amp_func=amp))
    geometry = []
    if a.medium:
        y1 = j1 * D
        geometry.append(mp.Block(center=mp.Vector3(0, to_meep_y((y1 + Ly) / 2), 0),
                                 size=mp.Vector3(mp.inf, Ly - y1, mp.inf), material=mp.Medium(index=1.5)))
    mp.verbosity(0)
    sim = mp.Simulation(cell_size=mp.Vector3(Lx, Ly, Lz), resolution=nl, sources=srcs, geometry=geometry,
                        boundary_layers=[mp.PML(npml * D, direction=mp.Y)], k_point=k_meep,
                        force_complex_fields=True, Courant=0.5, eps_averaging=bool(a.eps_averaging))
    t_start = time.time()
    sim.run(until=t_win0)
    tA = sim.meep_time()

    # DFT monitors (window starts now)
    zc = nl * D  # the C code's xy slice sits at k = zk = nl -> z = 1.0 (integer-z comps)
    slab = sim.add_dft_fields(list(COMPS.values()), [1.0], center=mp.Vector3(0, 0, zc + D / 4 - Lz / 2),
                              size=mp.Vector3(Lx, Ly - 2 * npml * D, D / 2), yee_grid=True, decimation_factor=1)
    ylist = sorted({npml + nl // 2, j0 + nl, j0 + nl + 1, j0 + 2 * nl, j0 + 5 * nl, j0 + 5 * nl + 1} |
                   ({j1 + nl, j1 + nl + 1, j1 + 2 * nl, j1 + 3 * nl} if a.medium else set()))
    planes = {}
    for j in ylist:
        # thin y-slab covering y = jΔ (integer comps) and y = (j+1/2)Δ (half comps)
        planes[j] = sim.add_dft_fields(list(COMPS.values()), [1.0],
                                       center=mp.Vector3(0, to_meep_y(j * D + D / 4), 0),
                                       size=mp.Vector3(Lx, D / 2, Lz), yee_grid=True, decimation_factor=1)
    flux = {}
    # R plane 1λ0 after the source (same in vacuum and medium runs so the normalization data match);
    # T plane 2λ0 inside the medium (vacuum run: same y, used as a second incident-flux plane)
    fl_y = {"R": (j0 + nl) * D + D / 4, "T": (j1 + 2 * nl) * D + D / 4}
    for key, y in fl_y.items():
        flux[key] = sim.add_flux(1.0, 0, 1, mp.FluxRegion(center=mp.Vector3(0, to_meep_y(y), 0),
                                                          size=mp.Vector3(Lx, 0, Lz)), decimation_factor=1)
    if a.norm:  # medium run: subtract the incident DFT fields on the R plane (standard Meep normalization)
        d = np.load(os.path.join(a.norm, "fluxdata_R.npz"))
        from meep.simulation import FluxData  # not re-exported at package level in Meep 1.34
        sim.load_minus_flux_data(flux["R"], FluxData(E=d["E"], H=d["H"]))
    sim.run(until=t_win1 - t_win0)
    tB = sim.meep_time()

    out = {}
    # slab (x-y slice): pick the z-layer that matches the C slice for each component
    for name, c in COMPS.items():
        arr = sim.get_dft_array(slab, c, 0)
        ox, oy, oz = OFF[name]
        xs = axis_coords(-Lx / 2, Lx / 2, D, ox == 0.5) + Lx / 2
        y_lo, y_hi = npml * D - Ly / 2, (Ny - npml) * D - Ly / 2
        ys = axis_coords(y_lo, y_hi, D, oy == 0.5) + Ly / 2
        z_lo, z_hi = zc - Lz / 2, zc + D / 2 - Lz / 2
        zs = axis_coords(z_lo, z_hi, D, oz == 0.5) + Lz / 2
        if arr.shape != (len(xs), len(ys), len(zs)):
            # Meep may drop the duplicated half-offset layer for zero-size / periodic dims; fall back to shape
            zs = zs[:arr.shape[2]]
        target = zc + oz * D
        kz_idx = int(np.argmin(np.abs(zs - target)))
        out[f"xy_{name}"] = arr[:, :, kz_idx]
        out[f"xy_{name}_x"] = xs[:arr.shape[0]]
        out[f"xy_{name}_y"] = ys[:arr.shape[1]]
        out[f"xy_{name}_z"] = zs[kz_idx]
    for j, obj in planes.items():
        for name, c in COMPS.items():
            arr = sim.get_dft_array(obj, c, 0)
            ox, oy, oz = OFF[name]
            y_lo = j * D - Ly / 2
            ys = axis_coords(y_lo, y_lo + D / 2, D, oy == 0.5) + Ly / 2
            target = (j + oy) * D
            jy = int(np.argmin(np.abs(ys[:arr.shape[1]] - target)))
            out[f"y{j}_{name}"] = arr[:, jy, :]
            out[f"y{j}_{name}_x"] = axis_coords(-Lx / 2, Lx / 2, D, ox == 0.5)[:arr.shape[0]] + Lx / 2
            out[f"y{j}_{name}_z"] = axis_coords(-Lz / 2, Lz / 2, D, oz == 0.5)[:arr.shape[2]] + Lz / 2
            out[f"y{j}_{name}_y"] = ys[jy]
    for key, f in flux.items():
        out[f"flux_{key}"] = np.array(mp.get_fluxes(f))
    os.makedirs(a.out, exist_ok=True)
    meta = dict(nl=nl, Delta=D, dt=0.5 * D, pol=a.pol, medium=a.medium, eps_averaging=a.eps_averaging,
                Lx=Lx, Ly=Ly, Lz=Lz, Ny=Ny, npml=npml, sf=sf, tf=tf, j0=j0, j1=j1 if a.medium else -1,
                k_point=[k_meep.x, k_meep.y, k_meep.z], J=J.tolist(), window=[tA, tB], planes=ylist,
                flux_y=fl_y, meep_version=mp.__version__, runtime_s=time.time() - t_start,
                window_periods=(tB - tA) / per)
    np.savez_compressed(os.path.join(a.out, "meep.npz"), **out)
    if not a.medium:  # incident DFT fields on the R plane, for the medium run's normalization
        fd = sim.get_flux_data(flux["R"])
        np.savez(os.path.join(a.out, "fluxdata_R.npz"), E=fd.E, H=fd.H)
    with open(os.path.join(a.out, "meep_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"meep {mp.__version__}: window [{tA:.3f}, {tB:.3f}] T0, runtime {meta['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
