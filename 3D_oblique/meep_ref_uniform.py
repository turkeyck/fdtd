"""Meep reference for gate 3 (A): converged R, T of the layered structures (run with the Meep interpreter).

    MEEP_PYTHON meep_ref_uniform.py --res 80 --pol s --struct F1 --out runs/nu/meep/A_F1_s_r80 [--cell full]

Thin 3D Bloch cell (4 x 4 pixels in x, z) with k_point = (kx/2π, 0, kz/2π) (Meep's k_point has no 2π factor;
Bloch phase exp(2πi k·r)). Source: planar current sheet J ∥ tangential E of the s/p wave, amp_func =
exp(2πi k·r), narrow-band GaussianSource(f0 = 1, fwidth = 0.1). A vacuum normalization run gives the incident
flux and the incident fields on the R plane (load_minus_flux_data). Only the f0 component is used (with a fixed
k_point, other frequencies are other angles). Subpixel averaging is on (Meep default) for second-order accuracy.
Writes meep_A.json with R, T, fluxes, resolution, runtime.
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

N_FILM, N_SUB = 2.0, 1.46
STACK = {"F1": [(2.0, 0.08)], "F5": [(2.0, 0.08), (1.46, 0.12), (2.0, 0.08), (1.46, 0.12), (2.0, 0.08)],
         "vac": []}
KX, KZ = math.pi, 2 * math.pi / 3              # m = n = 1, Lx = 2, Lz = 3 (SPEC A2)
PML, GAP_SRC_R, GAP_R_FILM, SUB, = 1.0, 0.5, 1.0, 1.5


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
    ap.add_argument("--res", type=int, required=True)
    ap.add_argument("--pol", default="s")
    ap.add_argument("--struct", default="F1")
    ap.add_argument("--cell", default="thin", choices=("thin", "full"))
    ap.add_argument("--norm", default="", help="normalization run directory (vacuum, same res/pol/cell)")
    ap.add_argument("--decay", type=float, default=1e-11)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = a.res
    layers = STACK[a.struct]
    t_stack = sum(d for _, d in layers)
    Ly = 2 * PML + GAP_SRC_R + GAP_R_FILM + t_stack + SUB + 0.5
    Ly = math.ceil(Ly * 40) / 40.0                      # a multiple of 1/40 (valid for res 40..320)
    if a.cell == "thin":
        Lx = Lz = 4.0 / res
    else:
        Lx, Lz = 2.0, 3.0
    y_src = -Ly / 2 + PML + 0.25
    y_R = y_src + GAP_SRC_R
    y_film = y_R + GAP_R_FILM
    y_T = y_film + t_stack + 0.75
    kp = mp.Vector3(KX / (2 * math.pi), 0, KZ / (2 * math.pi))

    def amp(p):
        return np.exp(2j * math.pi * (kp.x * p.x + kp.z * p.z))
    J = basis(a.pol)
    src = [mp.Source(mp.GaussianSource(1.0, fwidth=0.1), component=c, center=mp.Vector3(0, y_src, 0),
                     size=mp.Vector3(Lx, 0, Lz), amplitude=complex(v), amp_func=amp)
           for c, v in ((mp.Ex, J[0]), (mp.Ez, J[2])) if abs(v) > 1e-14]
    geometry = []
    if not a.norm == "" or a.struct != "vac":
        pass
    if a.struct != "vac" and a.norm:
        y = y_film
        for n, d in layers:
            geometry.append(mp.Block(center=mp.Vector3(0, y + d / 2, 0), size=mp.Vector3(mp.inf, d, mp.inf),
                                     material=mp.Medium(index=n)))
            y += d
        top = Ly / 2
        geometry.append(mp.Block(center=mp.Vector3(0, (y + top) / 2, 0), size=mp.Vector3(mp.inf, top - y, mp.inf),
                                 material=mp.Medium(index=N_SUB)))
    mp.verbosity(0)
    sim = mp.Simulation(cell_size=mp.Vector3(Lx, Ly, Lz), resolution=res, sources=src, geometry=geometry,
                        boundary_layers=[mp.PML(PML, direction=mp.Y)], k_point=kp, Courant=0.5)
    fR = sim.add_flux(1.0, 0, 1, mp.FluxRegion(center=mp.Vector3(0, y_R, 0), size=mp.Vector3(Lx, 0, Lz)))
    fT = sim.add_flux(1.0, 0, 1, mp.FluxRegion(center=mp.Vector3(0, y_T, 0), size=mp.Vector3(Lx, 0, Lz)))
    if a.norm:
        d = np.load(os.path.join(a.norm, "fluxdata_R.npz"))
        from meep.simulation import FluxData
        sim.load_minus_flux_data(fR, FluxData(E=d["E"], H=d["H"]))
    comp = mp.Ez if abs(J[2]) >= abs(J[0]) else mp.Ex
    t0 = time.time()
    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, comp, mp.Vector3(0, y_T, 0), a.decay))
    runtime = time.time() - t0
    os.makedirs(a.out, exist_ok=True)
    out = dict(res=res, pol=a.pol, struct=a.struct, cell=a.cell, Lx=Lx, Ly=Ly, Lz=Lz, k_point=[kp.x, kp.y, kp.z],
               flux_R=float(mp.get_fluxes(fR)[0]), flux_T=float(mp.get_fluxes(fT)[0]), meep_version=mp.__version__,
               runtime_s=runtime, meep_time=sim.meep_time(), y=dict(src=y_src, R=y_R, film=y_film, T=y_T))
    if not a.norm:
        fd = sim.get_flux_data(fR)
        np.savez(os.path.join(a.out, "fluxdata_R.npz"), E=fd.E, H=fd.H)
    else:
        with open(os.path.join(a.norm, "meep_A.json")) as fh:
            nrm = json.load(fh)
        out["flux_inc"] = nrm["flux_R"]
        out["R"] = -out["flux_R"] / nrm["flux_R"]
        out["T"] = out["flux_T"] / nrm["flux_R"]
    with open(os.path.join(a.out, "meep_A.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: out[k] for k in out if k in ("res", "R", "T", "runtime_s")}))


if __name__ == "__main__":
    main()
