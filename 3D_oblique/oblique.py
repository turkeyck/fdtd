"""Physical-units front end for fdtd3d_oblique: wavelength, incidence angle and grid in, E-field maps out.

    python3 oblique.py --wavelength 530nm --angle 30 --grid lambda/20
    python3 oblique.py --wavelength 1.55um --angle 45 --azimuth 20 --grid 25nm --pol s p

The solver works in normalized units (λ0 = 1, c = 1). In linear, non-dispersive media Maxwell's equations are
scale-invariant, so the wavelength only sets the unit conversion (Δ in nm, Δt in fs); the fields in units of λ0
are the same for every λ0. The angle is not free: with periodic x/z boundaries the tangential wavenumbers are
kx = 2πm/Lx, kz = 2πn/Lz, and Lx, Lz must be whole numbers of cells. The converter picks the (m, Lx), (n, Lz)
that reproduce the requested angle best and reports the angle actually simulated.

Angle convention: θ from +y (propagation axis), azimuth φ = atan2(kz, kx); sin θ = |k_t|/k0 (continuum
definition, the one Snell's law uses). The phase fronts on the Yee lattice travel at atan(|k_t|/ky_disc),
slightly smaller because ky_disc > ky_cont (numerical dispersion); both are reported.
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
import analyze as an  # noqa: E402
import fdtd_io  # noqa: E402
from fdtd_theory import Setup  # noqa: E402

C_LIGHT = 299792458.0
UNITS = {"nm": 1e-9, "um": 1e-6, "µm": 1e-6, "mm": 1e-3, "m": 1.0}


def parse_length(s):
    """'530nm', '0.53um', '5.3e-7m'; a bare number is nm. Returns metres."""
    m = re.fullmatch(r"\s*([0-9.eE+-]+)\s*([a-zµ]*)\s*", s)
    if not m:
        raise ValueError(f"cannot parse length '{s}'")
    unit = m.group(2) or "nm"
    if unit not in UNITS:
        raise ValueError(f"unknown unit '{unit}' (use nm, um, mm, m)")
    return float(m.group(1)) * UNITS[unit]


def parse_grid(s, lam):
    """'lambda/20', 'λ/20', '20' (cells per wavelength) or a length like '26.5nm'. Returns cells per λ0."""
    m = re.fullmatch(r"\s*(?:lambda|λ)\s*/\s*([0-9.]+)\s*", s)
    if m:
        n = float(m.group(1))
    elif re.fullmatch(r"\s*[0-9]+\s*", s):
        n = float(s)
    else:
        n = lam / parse_length(s)
    nl = int(round(n))
    if abs(n - nl) > 1e-6 * n:
        print(f"note: λ/Δ = {n:.6f} is not an integer; using Δ = λ/{nl} "
              f"({lam / nl * 1e9:.4f} nm instead of {lam / n * 1e9:.4f} nm)")
    if nl < 10:
        raise ValueError("fewer than 10 cells per wavelength is not supported")
    return nl


def quantize(frac, nl, max_cells, tol):
    """(m, N) with k/k0 = m*nl/N ≈ frac (N cells per period, m integer); frac = 0 -> (0, 2).

    Returns the smallest period N whose error |m*nl/N - frac| is within tol; if none fits within max_cells,
    the most accurate one."""
    if abs(frac) < 1e-12:
        return 0, 2
    cands = []
    for m in range(1, 64):
        N = int(round(m * nl / abs(frac)))
        if 2 <= N <= max_cells:
            cands.append((abs(m * nl / N - abs(frac)), N, m))
    if not cands:
        raise ValueError("no period fits within --max-cells; increase it")
    ok = [c for c in cands if c[0] <= tol]
    err, N, m = min(ok, key=lambda c: (c[1], c[0])) if ok else min(cands)
    return int(math.copysign(m, frac)), N


def plan(wavelength, angle_deg, nl, azimuth_deg=0.0, max_cells=400, sf_l=1.0, tf_l=6.0, npml=20, S=0.5,
         angle_tol=0.05):
    """Physical request -> solver parameters + everything needed to interpret the output."""
    th, ph = math.radians(angle_deg), math.radians(azimuth_deg)
    if not 0.0 < angle_deg < 90.0:
        raise ValueError("angle must be in (0, 90) degrees; normal incidence is not supported by the s/p basis")
    sx, sz = math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph)
    ktol = math.cos(th) * math.radians(angle_tol) / math.sqrt(2.0)  # per-axis k/k0 error ~ angle tolerance
    m, Nx = quantize(sx, nl, max_cells, ktol)
    n, Nz = quantize(sz, nl, max_cells, ktol)
    if m == 0 and n == 0:
        raise ValueError("angle too small for the allowed period; increase --max-cells")
    Lx, Lz = Nx / nl, Nz / nl
    kx, kz = m / Lx, n / Lz                                   # in units of k0 (k0 = 2π/λ0)
    kt = math.hypot(kx, kz)
    if kt >= 0.95:
        raise ValueError(f"realized sinθ = {kt:.4f} ≥ 0.95: the solver needs a 5% propagation margin")
    st = Setup(n_lambda=nl, S=S, Lx=Lx, Lz=Lz, m=m, n=n, pol="s")
    theta = math.degrees(math.asin(kt))
    theta_lat = math.degrees(math.atan2(kt * st.k0, st.ky))
    D = wavelength / nl
    return dict(
        wavelength_nm=wavelength * 1e9, frequency_THz=C_LIGHT / wavelength * 1e-12, period_fs=wavelength / C_LIGHT * 1e15,
        cells_per_wavelength=nl, Delta_nm=D * 1e9, dt_fs=S * D / C_LIGHT * 1e15, courant=S,
        requested=dict(angle_deg=angle_deg, azimuth_deg=azimuth_deg),
        realized=dict(angle_deg=theta, azimuth_deg=math.degrees(math.atan2(kz, kx)),
                      angle_error_deg=theta - angle_deg, lattice_phase_front_angle_deg=theta_lat),
        solver=dict(nl=nl, S=S, Lx=Lx, Lz=Lz, m=m, n=n, npml=npml, sf=int(round(sf_l * nl)), tf=int(round(tf_l * nl))),
        cells=dict(Nx=Nx, Nz=Nz, Ny=2 * npml + int(round(sf_l * nl)) + int(round(tf_l * nl))),
        size_nm=dict(Lx=Lx * wavelength * 1e9, Lz=Lz * wavelength * 1e9),
        ky_over_k0=dict(discrete=st.ky / st.k0, continuum=st.ky_cont / st.k0),
    )


def run(p, pol, outdir):
    s = p["solver"]
    per = int(round(s["nl"] / s["S"]))                       # steps per optical period
    meta = fdtd_io.run(outdir, pol=pol, inc="a", nsteps=70 * per, dft0=60 * per, dft1=70 * per, zk=0,
                       xi=0, energy_every=0, **s)
    return meta


def measured_angle(outdir, meta):
    """Phase-front direction fitted on the x–y slice (TF region), as an end-to-end check of the conversion."""
    st = an.setup_from_meta(meta)
    jlo, jhi = meta["j0"] + 2, meta["Ny"] - meta["npml"] - 1
    sets, _ = an.phase_sets(outdir, meta, st, ["xy"], jlo, jhi)
    k, rms, _, _ = an.wavefront_fit(sets)
    kt = math.hypot(k[0], st.kz)
    return math.degrees(math.atan2(kt, k[1])), rms


def plot(p, pol, outdir, meta, figpath, tiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    F, sl = fdtd_io.load_dft(outdir, "xy", meta)
    lam = p["wavelength_nm"]
    Ny = meta["Ny"]
    E = {c: np.real(F[c][:, :Ny]) for c in ("Ex", "Ey", "Ez")}           # instantaneous field at t = 0 (mod T)
    amp = np.sqrt(sum(np.abs(F[c][:, :Ny]) ** 2 for c in ("Ex", "Ey", "Ez")))
    ext = [0, Ny * meta["Delta"] * lam, 0, tiles * meta["Lx"] * lam]
    live = [c for c in ("Ex", "Ey", "Ez") if np.abs(F[c]).max() > 1e-6]
    zero = [c for c in ("Ex", "Ey", "Ez") if c not in live]
    panels = [(f"Re {c}  (t = 0)", np.tile(E[c], (tiles, 1)), "RdBu_r", (-1, 1)) for c in live]
    panels.append(("|E|  (phasor amplitude)", np.tile(amp, (tiles, 1)), "magma", (0, 1.05)))
    h = 13.0 * (ext[3] / ext[1]) * len(panels) + 1.6          # equal aspect: true angles on screen
    fig, axes = plt.subplots(len(panels), 1, figsize=(13, h), sharex=True)
    axes = np.atleast_1d(axes)
    yj0 = meta["j0"] * meta["Delta"] * lam
    ypml = (meta["npml"] * meta["Delta"] * lam, (Ny - meta["npml"]) * meta["Delta"] * lam)
    th = math.radians(p["realized"]["angle_deg"])
    for ax, (title, data, cmap, lim) in zip(axes, panels):
        im = ax.imshow(data, origin="lower", aspect="equal", extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1],
                       interpolation="nearest")
        for yy in ypml:
            ax.axvspan(0 if yy == ypml[0] else yy, yy if yy == ypml[0] else ext[1], color="gray", alpha=0.25, lw=0)
        ax.axvline(yj0, color="lime", lw=1.2, ls="--")
        ax.set_ylabel("x  [nm]")
        ax.set_title(title, fontsize=10, loc="left")
        fig.colorbar(im, ax=ax, pad=0.01, label="E / E0")
    ax = axes[0]
    ax.text(ypml[0] / 2, ext[3] * 0.92, "PML", ha="center", va="top", fontsize=8, color="k")
    ax.text(ypml[1] + (ext[1] - ypml[1]) / 2, ext[3] * 0.92, "PML", ha="center", va="top", fontsize=8, color="k")
    ax.text(yj0, ext[3] * 0.92, " TF/SF plane", color="lime", fontsize=8, va="top")
    L = 1.2 * lam
    x0, y0 = ext[3] * 0.25, yj0 + 0.8 * lam
    for ax in axes:
        ax.annotate("", xy=(y0 + L * math.cos(th), x0 + L * math.sin(th)), xytext=(y0, x0),
                    arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    axes[0].text(y0 + L * math.cos(th), x0 + L * math.sin(th), f" k  (θ = {p['realized']['angle_deg']:.2f}°)",
                 fontsize=9, va="bottom")
    axes[-1].set_xlabel("y  [nm]   (propagation axis)")
    r = p["realized"]
    fig.suptitle(f"{pol}-polarized plane wave: λ = {lam:g} nm, θ = {r['angle_deg']:.3f}° (requested "
                 f"{p['requested']['angle_deg']:g}°), φ = {r['azimuth_deg']:.1f}°, Δ = λ/{p['cells_per_wavelength']} = "
                 f"{p['Delta_nm']:.3g} nm, Δt = {p['dt_fs']:.4g} fs\n"
                 f"x–y slice at z = 0, x shown over {tiles} period(s) of Lx = {p['size_nm']['Lx']:.4g} nm (PBC), "
                 f"DFT at f = {p['frequency_THz']:.2f} THz over periods 60–70"
                 + (f"; {', '.join(zero)} ≡ 0 for this polarization" if zero else ""), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(figpath, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wavelength", default="530nm", help="vacuum wavelength, e.g. 530nm, 1.55um [530nm]")
    ap.add_argument("--angle", type=float, default=30.0, help="incidence angle θ from +y, degrees [30]")
    ap.add_argument("--azimuth", type=float, default=0.0, help="azimuth φ = atan2(kz, kx), degrees [0]")
    ap.add_argument("--grid", default="lambda/20", help="lambda/N, N, or a cell size like 26.5nm [lambda/20]")
    ap.add_argument("--pol", nargs="+", default=["s", "p"], choices=["s", "p"], help="polarization(s) [s p]")
    ap.add_argument("--tf", type=float, default=6.0, help="total-field length in wavelengths [6]")
    ap.add_argument("--sf", type=float, default=1.0, help="scattered-field length in wavelengths [1]")
    ap.add_argument("--max-cells", type=int, default=400, help="largest allowed period Nx or Nz in cells [400]")
    ap.add_argument("--angle-tol", type=float, default=0.05,
                    help="accept the smallest period whose angle error is within this many degrees [0.05]")
    ap.add_argument("--tiles", type=int, default=0, help="x periods shown in the plot [auto: ≥ 2 wavelengths]")
    ap.add_argument("--dry-run", action="store_true", help="only print the conversion")
    a = ap.parse_args()

    lam = parse_length(a.wavelength)
    nl = parse_grid(a.grid, lam)
    p = plan(lam, a.angle, nl, a.azimuth, a.max_cells, a.sf, a.tf, angle_tol=a.angle_tol)
    print(json.dumps(p, indent=2, ensure_ascii=False))
    if abs(p["realized"]["angle_error_deg"]) > a.angle_tol:
        print(f"note: the periodic grid realizes θ = {p['realized']['angle_deg']:.4f}° "
              f"(error {p['realized']['angle_error_deg']:+.4f}°); raise --max-cells or refine --grid to reduce it")
    if a.dry_run:
        return 0
    fdtd_io.build()
    tag = f"lam{p['wavelength_nm']:g}nm_th{a.angle:g}_phi{a.azimuth:g}_N{nl}"
    os.makedirs(os.path.join(ROOT, "figures", "oblique"), exist_ok=True)
    Ly = p["cells"]["Ny"] / nl
    tiles = a.tiles or max(1, round(0.4 * Ly / p["solver"]["Lx"]))  # enough x periods to see the fronts
    for pol in a.pol:
        outdir = os.path.join(ROOT, "runs", "oblique", f"{tag}_{pol}")
        meta = run(p, pol, outdir)
        th_meas, rms = measured_angle(outdir, meta)
        fig = os.path.join(ROOT, "figures", "oblique", f"{tag}_{pol}_E.png")
        plot(p, pol, outdir, meta, fig, tiles)
        p.setdefault("check", {})[pol] = dict(measured_phase_front_angle_deg=th_meas, fit_rms_rad=rms, figure=fig,
                                              run=outdir)
        print(f"{pol}: measured phase-front angle {th_meas:.6f}° (lattice theory "
              f"{p['realized']['lattice_phase_front_angle_deg']:.6f}°), fit RMS {rms:.1e} rad -> {fig}")
    with open(os.path.join(ROOT, "figures", "oblique", f"{tag}.json"), "w", encoding="utf-8") as fh:
        json.dump(p, fh, indent=2, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
