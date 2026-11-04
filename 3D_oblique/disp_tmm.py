"""References for dispersive layered stacks (docs/DISPERSION_SUBPIXEL.md).

1. stack files:   write_stack(path, lambda0_um, layers) -- layers [(material name, y0, y1)] in lambda0 units,
                   relative to the TF/SF plane; materials from materials/models.json (Lorentz fits).
2. eps_fd:        effective permittivity of the solver's Lorentz ADE at angular frequency w (exact for the
                   discrete-time scheme: D^n = eps_inf E^n + sum P^n, P^n = chi(w) E^n).
3. lattice_rt:    exact Yee-lattice R, T at w for a run (stack_used.json + grid_used.json + meta.json): the
                   subpixel / staircase weights, eps_fd per material, the discrete K~t and w~ -- the reference that
                   the FDTD must reproduce up to CPML reflection and DFT truncation.
4. tmm_cont:      continuous transfer-matrix R, T (complex eps) -- the "physical" answer (Lumerical STACK equivalent)
                   with either the fitted model or the raw experimental n, k.
Units: lambda0 = 1, c = 1; f in c/lambda0 (f = lambda0 / lambda).
"""
import cmath
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "materials"))
import fit as mfit  # noqa: E402

MODELS = os.path.join(ROOT, "materials", "models.json")


def models():
    return json.load(open(MODELS, encoding="utf-8"))


# ----------------------------------------------------------------------------- stack file
def write_stack(path, lambda0_um, layers, pml_index_hi=None):
    """layers: [(name, y0, y1)] with y relative to the TF/SF plane (lambda0 units); y1 may be large (into the PML)."""
    md = models()
    names = []
    for nm, _, _ in layers:
        if nm not in names:
            names.append(nm)
    out = {"stack_format": 1, "lambda0_um": lambda0_um, "nmat": len(names),
           "names": names, "layer_mat": [names.index(nm) + 1 for nm, _, _ in layers],
           "layer_y0": [float(a) for _, a, _ in layers], "layer_y1": [float(b) for _, _, b in layers]}
    for q, nm in enumerate(names):
        e = md[nm]
        out[f"mat{q + 1}"] = {"eps_inf": e["eps_inf"], "de": [p["de"] for p in e["poles"]],
                              "nu0": [p["nu0"] for p in e["poles"]], "g": [p["g"] for p in e["poles"]]}
    if pml_index_hi is not None:
        out["pml_index_hi"] = pml_index_hi
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    return path


# ----------------------------------------------------------------------------- permittivities
def eps_cont_model(name, f, lambda0_um):
    """Continuous fitted model at f (c/lambda0)."""
    e = models()[name]
    nu = np.asarray(f, float) / lambda0_um          # 1/um
    return mfit.eps_model(nu, *mfit.model_of(e))


def eps_data(name, f, lambda0_um):
    """Raw experimental eps = (n + i k)^2 at f (interpolated table or formula)."""
    e = models()[name]
    _, n_of, _ = mfit.read_yml(os.path.join(ROOT, "materials", "data", e["source"]))
    lam_um = lambda0_um / np.asarray(f, float)
    n = n_of(lam_um)
    return n * n


def eps_fd(mat, w, dt):
    """Discrete-time ADE permittivity (solver units). mat: entry of stack_used.json 'materials'."""
    z = cmath.exp(-1j * w * dt)
    e = mat["eps_inf"]
    for w0, de, g in zip(mat["w0"], mat["de"], mat["g"]):
        a = 0.5 * g * dt
        e += de * w0 * w0 * dt * dt / ((1 + a) * z - 2 + w0 * w0 * dt * dt + (1 - a) / z)
    return e


def eps_cont_used(mat, w):
    e = mat["eps_inf"]
    for w0, de, g in zip(mat["w0"], mat["de"], mat["g"]):
        e += de * w0 * w0 / (w0 * w0 - w * w - 1j * g * w)
    return e


# ----------------------------------------------------------------------------- exact lattice (complex)
def _phi(c):
    """Root of cos(phi) = c for an outgoing / decaying wave e^{i phi m} (Im phi > 0, or Re phi > 0 if lossless)."""
    p = cmath.acos(c)
    if p.imag < -1e-15 or (abs(p.imag) <= 1e-15 and p.real < 0):
        p = -p
    return p


def _scatter(w, gm, gp, q, gL, wL, qL, gR, wR, qR):
    """(1/w_m)[(u_{m+1}-u_m)/gp_m - (u_m-u_{m-1})/gm_m] + q_m u_m = 0, complex coefficients; incident 1 from left.
    Returns (r, t, phiL, phiR)."""
    M = len(w) - 1
    cL = 1.0 - qL * wL * gL / 2.0
    cR = 1.0 - qR * wR * gR / 2.0
    phL, phR = _phi(cL), _phi(cR)
    a = 1.0 / (w * gm)
    c = 1.0 / (w * gp)
    b = (-(a + c) + q).astype(complex)
    a = a.astype(complex); c = c.astype(complex)
    rhs = np.zeros(M + 1, complex)
    b[0] += a[0] * cmath.exp(1j * phL)
    rhs[0] = a[0] * 2j * cmath.sin(phL)
    b[M] += c[M] * cmath.exp(1j * phR)
    # Thomas
    n = M + 1
    cp = np.empty(n, complex); rp = np.empty(n, complex)
    cp[0] = c[0] / b[0]; rp[0] = rhs[0] / b[0]
    for m in range(1, n):
        den = b[m] - a[m] * cp[m - 1]
        cp[m] = c[m] / den if m < n - 1 else 0.0
        rp[m] = (rhs[m] - a[m] * rp[m - 1]) / den
    u = np.empty(n, complex)
    u[-1] = rp[-1]
    for m in range(n - 2, -1, -1):
        u[m] = rp[m] - cp[m] * u[m + 1]
    return u[0] - 1.0, u[M], phL, phR


def lattice_rt(outdir, f):
    """Exact Yee-lattice R(f), T(f) of a stack run (T = nan when the last material is lossy)."""
    meta = json.load(open(os.path.join(outdir, "meta.json")))
    used = json.load(open(os.path.join(outdir, "grid_used.json")))
    st = json.load(open(os.path.join(outdir, "stack_used.json")))
    dt = meta["dt"]
    h = np.array(used["hy"])
    Ny = len(h)
    dprim = np.empty(Ny + 1)
    dprim[1:-1] = 0.5 * (h[:-1] + h[1:]); dprim[0], dprim[-1] = h[0], h[-1]
    fT, fN = np.array(st["fT"]), np.array(st["fN"])
    mats = st["materials"]
    Kt = math.hypot(meta["Ktilde"][0], meta["Ktilde"][2])
    j_lo, j_hi = meta["j0"], Ny - meta["npml_hi"] - 1
    pol = meta["pol"]
    R, T = [], []
    for fi in np.atleast_1d(f):
        w = 2 * math.pi * fi
        wt2 = (2.0 / dt * math.sin(w * dt / 2.0)) ** 2
        em = np.array([eps_fd(m, w, dt) for m in mats])
        et = fT @ em                                   # tangential: parallel (arithmetic) mixing
        en = 1.0 / (fN @ (1.0 / em))                    # normal: series (harmonic) mixing
        if wt2 <= Kt * Kt:
            R.append(np.nan); T.append(np.nan); continue
        if pol == "s":
            js = np.arange(j_lo, j_hi + 1)
            w_ = dprim[js]
            gm, gp = h[js - 1], h[js]
            q = wt2 * et[js] - Kt * Kt
            gL, gR = h[j_lo - 1], h[j_hi]
            eL, eR = et[j_lo - 1], et[j_hi + 1]
            r, t, phL, phR = _scatter(w_, gm, gp, q, gL, h[j_lo - 1], wt2 * eL - Kt * Kt, gR, h[j_hi], wt2 * eR - Kt * Kt)
            sR, sL = phR, phL
            fac = (cmath.sin(sR) / gR).real / (math.sin(sL.real) / gL)
        else:
            js = np.arange(j_lo, j_hi)                 # dual nodes j+1/2; links through primal nodes j, j+1
            w_ = h[js]
            gm = et[js] * dprim[js]
            gp = et[js + 1] * dprim[js + 1]
            q = wt2 - Kt * Kt / en[js]
            eL, eR = en[j_lo - 1], en[j_hi]
            gL, gR = gm[0], gp[-1]
            r, t, phL, phR = _scatter(w_, gm, gp, q, gL, h[j_lo - 1], wt2 - Kt * Kt / eL, gR, h[j_hi], wt2 - Kt * Kt / eR)
            fac = (cmath.sin(phR) / gR).real / (math.sin(phL.real) / gL)
        R.append(abs(r) ** 2)
        lossy = abs(complex(eR).imag) > 1e-6 * abs(eR)   # numerically lossless exit medium (e.g. SiO2 fit, Im ~1e-8)
        T.append(np.nan if lossy else abs(t) ** 2 * fac)
    return np.array(R), np.array(T)


# ----------------------------------------------------------------------------- continuous TMM (complex)
def tmm_cont(eps_layers, d_layers, eps_in, eps_out, kt, k0, pol):
    """Characteristic-matrix R, T for planar layers (complex eps), field e^{i(kt x + ky y - w t)}.
    Returns (R, T); T = nan for an absorbing exit medium."""
    def ky(e):
        v = cmath.sqrt(e * k0 * k0 - kt * kt)
        return v if v.imag >= 0 else -v

    def adm(e):  # tangential H / tangential E admittance (s: ky/(w mu); p: w eps/ky), common factors dropped
        return ky(e) if pol == "s" else e * k0 * k0 / ky(e)
    M = np.eye(2, dtype=complex)
    for e, d in zip(eps_layers, d_layers):
        k, Y = ky(e), adm(e)
        c, s = cmath.cos(k * d), cmath.sin(k * d)
        M = M @ np.array([[c, -1j * s / Y], [-1j * Y * s, c]])
    Y0, Ys = adm(eps_in), adm(eps_out)
    B, C = M @ np.array([1.0, Ys])
    r = (Y0 * B - C) / (Y0 * B + C)
    t = 2 * Y0 / (Y0 * B + C)
    R = abs(r) ** 2
    T = (Ys.real / Y0.real) * abs(t) ** 2 if abs(complex(eps_out).imag) < 1e-12 else np.nan
    return R, T


def stack_cont(layers, sub, f, lambda0_um, kt, pol, source="model"):
    """Continuous R(f) for [(name, d)] on substrate 'sub' (vacuum ambient), eps from the fit or the raw data."""
    R, T = [], []
    for fi in np.atleast_1d(f):
        epsf = (lambda nm: complex(eps_cont_model(nm, [fi], lambda0_um)[0])) if source == "model" else \
               (lambda nm: complex(eps_data(nm, [fi], lambda0_um)[0]))
        r, t = tmm_cont([epsf(nm) for nm, _ in layers], [d for _, d in layers], 1.0, epsf(sub), kt,
                        2 * math.pi * fi, pol)
        R.append(r); T.append(t)
    return np.array(R), np.array(T)
