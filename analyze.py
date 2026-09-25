"""Post-processing and validation library for fdtd3d_oblique (numpy + matplotlib only).

Library functions are used by tests/level*.py; the CLI runs whole validation levels:
    python3 analyze.py --level 0|1|2|all
Conventions: DFT phasor F with f(t) = Re[F exp(-i w0 t)]; a forward plane wave has F = A exp(i k.r_c)
at each component's own Yee point r_c (derivation.md §1, §2.4).
"""
import argparse
import math
import os
import subprocess
import sys

import numpy as np

import fdtd_io
from fdtd_theory import ECOMP, HCOMP, Setup, group_velocity, interp_factors

ROOT = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------- basics
def setup_from_meta(meta, n_ref=1.0):
    return Setup(n_lambda=meta["nl"], S=meta["S"], Lx=meta["Lx"], Lz=meta["Lz"], m=meta["m"], n=meta["n"],
                 pol=meta["pol"], n_ref=n_ref)


def j_index(meta, sl):
    """Array of the y-index j of every sample in slice sl (broadcast to the slice shape)."""
    shape = fdtd_io.slice_shape(meta, sl)
    if sl["type"] == 0:
        return np.full(shape, sl["idx"])
    if sl["type"] == 1:
        return np.broadcast_to(np.arange(shape[1])[None, :], shape)
    return np.broadcast_to(np.arange(shape[0])[:, None], shape)


def select(meta, sl, comp, jlo, jhi):
    """Mask of valid samples of comp with jlo <= j <= jhi."""
    j = j_index(meta, sl)
    return fdtd_io.valid_j_mask(meta, sl, comp) & (j >= jlo) & (j <= jhi)


def plane_wave(meta, sl, comp, k):
    x, y, z = fdtd_io.slice_coords(meta, sl, comp)
    return np.exp(1j * (k[0] * x + k[1] * y + k[2] * z))


def mode_amplitude(F, meta, sl, comp, k, mask=None):
    """Least-squares amplitude A of F ~ A exp(i k.r) over the (masked) samples."""
    w = np.conj(plane_wave(meta, sl, comp, k)) * F
    return complex(w[mask].mean() if mask is not None else w.mean())


def k_forward(st):
    return np.array([st.kx, st.ky, st.kz])


def k_backward(st):
    return np.array([st.kx, -st.ky, st.kz])


# ----------------------------------------------------------------------------- echo
def echo_from_sf(out, meta, st, jplane):
    """Backward-wave amplitudes B_c on an SF y-plane (only the scattered/echo field lives there).

    Returns (B dict, relative residual of the single-plane-wave model, |B_E| / |E0|)."""
    F, sl = fdtd_io.load_dft(out, f"y{jplane}", meta)
    kb = k_backward(st)
    B, res, scale = {}, 0.0, 0.0
    for c in ECOMP + HCOMP:
        B[c] = mode_amplitude(F[c], meta, sl, c, kb)
        res = max(res, float(np.abs(F[c] - B[c] * plane_wave(meta, sl, c, kb)).max()))
        scale = max(scale, abs(B[c]))
    e0 = np.linalg.norm(st.amplitudes()[0])
    bE = math.sqrt(sum(abs(B[c]) ** 2 for c in ECOMP))
    return B, (res / scale if scale > 0 else 0.0), bE / e0


def subtract_echo(F, meta, sl, st, B):
    kb = k_backward(st)
    return {c: F[c] - B[c] * plane_wave(meta, sl, c, kb) for c in F}


# ----------------------------------------------------------------------------- wavefront fit
def unwrap2d(ph):
    return np.unwrap(np.unwrap(ph, axis=1), axis=0)


def wavefront_fit(sets):
    """3D least squares phi = k.r + phi0_s with one offset per (component, slice) set.

    sets: list of (x, y, z, phi) 1D arrays. Returns k (3,), residual RMS [rad], max |residual|."""
    ns = len(sets)
    rows = sum(len(s[3]) for s in sets)
    A = np.zeros((rows, 3 + ns))
    b = np.zeros(rows)
    r0 = 0
    for q, (x, y, z, ph) in enumerate(sets):
        n = len(ph)
        A[r0:r0 + n, 0], A[r0:r0 + n, 1], A[r0:r0 + n, 2] = x, y, z
        A[r0:r0 + n, 3 + q] = 1.0
        b[r0:r0 + n] = ph
        r0 += n
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = b - A @ sol
    return sol[:3], float(np.sqrt(np.mean(resid ** 2))), float(np.abs(resid).max()), resid


def phase_sets(out, meta, st, names, jlo, jhi, B=None, min_rel=0.05):
    """Unwrapped DFT phase of every significant component on each named slice, TF samples jlo..jhi."""
    e, h = st.amplitudes()
    amp = dict(zip(ECOMP, np.abs(e))) | dict(zip(HCOMP, np.abs(h)))
    sets, labels = [], []
    for name in names:
        F, sl = fdtd_io.load_dft(out, name, meta)
        if B is not None:
            F = subtract_echo(F, meta, sl, st, B)
        for c in ECOMP + HCOMP:
            grp = ECOMP if c in ECOMP else HCOMP
            if amp[c] < min_rel * max(amp[g] for g in grp):
                continue  # component identically ~0 for this polarization (e.g. Ey for s)
            mask = select(meta, sl, c, jlo, jhi)
            if not mask.any():
                continue
            x, y, z = fdtd_io.slice_coords(meta, sl, c)
            ph = np.angle(F[c])
            if sl["type"] in (1, 2):
                # restrict to the contiguous y-band first so unwrapping never crosses excluded rows
                jj = j_index(meta, sl)
                if sl["type"] == 1:
                    cols = np.where((jj[0] >= jlo) & (jj[0] <= jhi) & mask[0])[0]
                    sub = unwrap2d(ph[:, cols])
                    xs, ys, zs = x[:, cols], y[:, cols], z[:, cols]
                else:
                    rws = np.where((jj[:, 0] >= jlo) & (jj[:, 0] <= jhi) & mask[:, 0])[0]
                    sub = unwrap2d(ph[rws, :].T).T
                    xs, ys, zs = x[rws, :], y[rws, :], z[rws, :]
            else:
                sub = unwrap2d(ph)
                xs, ys, zs = x, y, z
            sets.append((xs.ravel(), ys.ravel(), zs.ravel(), sub.ravel()))
            labels.append((name, c))
    return sets, labels


# ----------------------------------------------------------------------------- divergence (plane pairs)
def div_from_planes(Fa, Fb, D):
    """div E at the integer nodes of plane b (= a+1) and div H at the cell centres between a and b.

    Fa, Fb: DFT dicts of consecutive y-planes j and j+1 (shape Nx x Nz, PBC in both axes)."""
    divE = ((Fb["Ex"] - np.roll(Fb["Ex"], 1, axis=0)) + (Fb["Ey"] - Fa["Ey"])
            + (Fb["Ez"] - np.roll(Fb["Ez"], 1, axis=1))) / D
    divH = ((np.roll(Fa["Hx"], -1, axis=0) - Fa["Hx"]) + (Fb["Hy"] - Fa["Hy"])
            + (np.roll(Fa["Hz"], -1, axis=1) - Fa["Hz"])) / D
    return divE, divH


# ----------------------------------------------------------------------------- flux and Poynting
def conserved_flux(F, D):
    """Discrete-conserved y-flux through the plane between integer y=j and y=j+1/2 (derivation.md §9).

    Phi = 1/2 Re sum_{i,k} [Ez(j) Hx*(j+1/2) - Ex(j) Hz*(j+1/2)] Delta^2."""
    return 0.5 * float(np.real(np.sum(F["Ez"] * np.conj(F["Hx"]) - F["Ex"] * np.conj(F["Hz"])))) * D * D


def to_cell_centres(Fa, Fb):
    """Average each component to the cell centres (i+1/2, j+1/2, k+1/2) between planes a (j) and b (j+1)."""
    ax = lambda f: 0.5 * (f + np.roll(f, -1, axis=0))  # noqa: E731
    az = lambda f: 0.5 * (f + np.roll(f, -1, axis=1))  # noqa: E731
    ay = lambda fa, fb: 0.5 * (fa + fb)  # noqa: E731
    return {
        "Ex": az(ay(Fa["Ex"], Fb["Ex"])),
        "Ey": ax(az(Fa["Ey"])),
        "Ez": ax(ay(Fa["Ez"], Fb["Ez"])),
        "Hx": ax(Fa["Hx"]),
        "Hy": ay(Fa["Hy"], Fb["Hy"]),
        "Hz": az(Fa["Hz"]),
    }


def poynting(Fc):
    E = np.stack([Fc["Ex"], Fc["Ey"], Fc["Ez"]])
    H = np.stack([Fc["Hx"], Fc["Hy"], Fc["Hz"]])
    return 0.5 * np.real(np.cross(E, np.conj(H), axis=0))


def interp_theory(st):
    """Exact co-located (cell-centre) amplitudes of the discrete plane wave and its S and |H'|/|E'|."""
    a = interp_factors(st)
    e, h = st.amplitudes()
    Ei = np.array([a["Ex"] * e[0], a["Ey"] * e[1], a["Ez"] * e[2]])
    Hi = np.array([a["Hx"] * h[0], a["Hy"] * h[1], a["Hz"] * h[2]])
    return Ei, Hi, 0.5 * np.cross(Ei, Hi), np.linalg.norm(Hi) / np.linalg.norm(Ei)


def angle_deg(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    c = np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1.0, 1.0)
    return math.degrees(math.acos(c)) if c < 1.0 else 0.0


def theory_flux(st, meta):
    e, h = st.amplitudes()
    return 0.5 * math.cos(st.ky * st.Delta / 2.0) * (e[2] * h[0] - e[0] * h[2]) * meta["Lx"] * meta["Lz"]


def vg_direction(st):
    return group_velocity(st)


# ----------------------------------------------------------------------------- CLI
LEVEL_SCRIPTS = {
    "0": ["tests/level0_exact_injection.py"],
    "1": ["tests/level0b_multistep.py", "tests/level_cpml_stability.py", "tests/level1_leakage.py",
          "tests/level1_all.py"],
    "2": ["tests/level2_fresnel.py"],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", default="all", choices=["0", "1", "2", "all"])
    args = ap.parse_args()
    levels = ["0", "1", "2"] if args.level == "all" else [args.level]
    status = 0
    for lv in levels:
        for script in LEVEL_SCRIPTS[lv]:
            print(f"\n=== {script} ===", flush=True)
            r = subprocess.run([sys.executable, os.path.join(ROOT, script)], cwd=ROOT)
            if r.returncode != 0:
                print(f"*** {script} FAILED (exit {r.returncode}); stopping.")
                return r.returncode
    return status


if __name__ == "__main__":
    sys.exit(main())
