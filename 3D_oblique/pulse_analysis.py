"""Broadband (src=g, nfreq>0) post-processing: R(f), T(f) at fixed (kx, kz) from spec_proj.bin / spec_aux.bin.

Spectra layout (fdtd3d_oblique.c write_spec): complex128 [comp][freq][j = 0..Ny], comps Ex Ey Ez Hx Hy Hz;
E rows sit on primal y nodes y_j, Hx/Hz rows on dual nodes y_{j+1/2}. Each row is the x-z Floquet projection of the
field on exp(i(kx x + kz z)), so a layered structure leaves exactly one transverse mode per row.

Flux through the dual face between y_j and y_{j+1/2}: Phi_j = 1/2 Re[Ez_j conj(Hx_{j+1/2}) - Ex_j conj(Hz_{j+1/2})]
(the discrete-conserved y flux of a single Floquet mode; constant in j in a lossless layer).
T = Phi(transmitted rows) / Phi(aux line).
R = |B|^2 / |A|^2: forward (A) / backward (B) split of the main-grid TOTAL field on the uniform vacuum rows of the TF
region (between j0 and the first graded cell or interface), two consecutive primal samples per estimate,
ky from the local Yee dispersion. This is the reflection of the whole y grid seen from j0, which is what
tmm.tmm_discrete computes. (On a nonuniform y grid the vacuum aux line carries its own grading echo r_g, so the SF-region
flux R_sf = -Phi(SF)/Phi(aux) is |r - r_g|^2 instead; reported as INFO only -- results/nu_failures.json, gate P3-1.)
Coordinates come from grid_used.json (never idx*Delta).
"""
import json
import math
import os

import numpy as np

import tmm


def load(outdir):
    meta = json.load(open(os.path.join(outdir, "meta.json")))
    if meta.get("stack", {}).get("file"):
        meta["stack"]["_dir"] = outdir
    used = json.load(open(os.path.join(outdir, "grid_used.json")))
    sp = meta["spec"]
    nf, SY = sp["nfreq"], meta["Ny"] + 1
    f = np.array([sp["fmin"]] if nf == 1 else np.linspace(sp["fmin"], sp["fmax"], nf))

    def rd(name):
        p = os.path.join(outdir, name)
        return np.fromfile(p, dtype=np.complex128).reshape(6, nf, SY) if os.path.exists(p) else None
    return dict(meta=meta, used=used, f=f, proj=rd("spec_proj.bin"), aux=rd("spec_aux.bin"))


def flux_rows(S):
    """Phi_j for every frequency and row j = 0..Ny-1 (shape nf x Ny)."""
    Ex, Ez, Hx, Hz = S[0], S[2], S[3], S[5]
    return 0.5 * np.real(Ez[:, :-1] * np.conj(Hx[:, :-1]) - Ex[:, :-1] * np.conj(Hz[:, :-1]))


def regions(meta, used):
    """Row ranges (inclusive) away from the TF/SF plane, the interface layers and the PML (2-cell margins)."""
    Ny, j0, lo, hi = meta["Ny"], meta["j0"], meta["npml_lo"], meta["npml_hi"]
    eps = np.array(used["eps_n"])
    jump = [j for j in range(1, Ny) if eps[j] != eps[j - 1]]
    first_if = jump[0] if jump else Ny - hi
    last_if = jump[-1] if jump else None
    stk = meta.get("stack", {})
    if stk.get("file"):                          # dispersive stack: layers start at j_first
        first_if = stk["j_first"]
        last_if = None
        su = json.load(open(os.path.join(stk["_dir"], "stack_used.json"))) if "_dir" in stk else None
        if su is not None:                       # T rows only in a (numerically) lossless exit medium
            import disp_tmm
            ylast = max(l[1] for l in su["layers"])
            mexit = [l[0] for l in su["layers"] if l[1] == ylast][0]
            w = 2 * math.pi * meta["pulse"]["fc"]
            if abs(disp_tmm.eps_fd(su["materials"][mexit], w, meta["dt"]).imag) < 1e-5:
                y = np.array(used["y"])
                last_if = int(np.searchsorted(y, ylast)) + 1
    sf = (lo + 2, j0 - 2)
    inc = (j0 + 2, first_if - 2)
    hy = np.array(used["hy"])
    jg = next((j for j in range(j0, Ny) if abs(hy[j] / hy[j0] - 1.0) > 1e-12), Ny)
    vac = (j0 + 2, min(first_if, jg) - 3)        # uniform vacuum TF rows: both samples j, j+1 inside
    tr = (last_if + 2, Ny - hi - 2) if last_if is not None else None
    return dict(sf=sf, inc=inc, tr=tr, vac=vac)


def split_fb(S, f, used, meta, rows):
    """|A|^2, |B|^2 (summed over Ex, Ez) of A e^{iky y} + B e^{-iky y} for each frequency, averaged over the
    two-sample estimates at rows (j, j+1), j in rows; also the relative spread of B/A over rows."""
    y, hy = np.array(used["y"]), np.array(used["hy"])
    Kt = math.hypot(meta["Ktilde"][0], meta["Ktilde"][2])
    R, spread = np.full(len(f), np.nan), np.full(len(f), np.nan)
    for q, fi in enumerate(f):
        w = 2 * math.pi * fi
        vals = []
        for j in range(rows[0], rows[1] + 1):
            h = y[j + 1] - y[j]
            try:
                k = tmm.ky_local(h, 1.0, meta["dt"], Kt, w)
            except ValueError:
                break
            c, sK = 2 * math.cos(k * h / 2), 2j * math.sin(k * h / 2)
            a2 = b2 = 0.0
            for comp in (0, 2):
                u0, u1 = S[comp, q, j], S[comp, q, j + 1]
                A = 0.5 * ((u0 + u1) / c + (u1 - u0) / sK)
                B = 0.5 * ((u0 + u1) / c - (u1 - u0) / sK)
                a2 += abs(A) ** 2
                b2 += abs(B) ** 2
            vals.append(b2 / a2)
        if vals:
            R[q], spread[q] = np.mean(vals), np.ptp(vals)
    return R, spread


def measure(outdir):
    d = load(outdir)
    m, u = d["meta"], d["used"]
    rg = regions(m, u)
    Pin = flux_rows(d["aux"])
    Pm = flux_rows(d["proj"])
    a, b = rg["inc"]
    pin = Pin[:, a:b + 1].mean(axis=1)
    pin_var = np.abs(Pin[:, a:b + 1] - pin[:, None]).max(axis=1) / np.abs(pin)
    a, b = rg["sf"]
    R_sf = (-Pm[:, a:b + 1] / pin[:, None]).mean(axis=1)
    R, R_spread = split_fb(d["proj"], d["f"], u, m, rg["vac"])
    out = dict(f=d["f"], R=R, R_spread=R_spread, R_sf=R_sf, pin=pin, pin_var=pin_var,
               spec_inc=np.abs(d["aux"][0, :, m["j0"]]) ** 2 + np.abs(d["aux"][2, :, m["j0"]]) ** 2, regions=rg)
    if rg["tr"] is not None:
        a, b = rg["tr"]
        T_rows = Pm[:, a:b + 1] / pin[:, None]
        out["T"] = T_rows.mean(axis=1)
        out["T_spread"] = np.ptp(T_rows, axis=1)
    return out


def theory(outdir, f=None):
    """Exact Yee-lattice R(f), T(f) (tmm.tmm_discrete; disp_tmm.lattice_rt for dispersive stacks) for the run's
    y grid, dt and fixed tangential K~."""
    d = load(outdir)
    m, u = d["meta"], d["used"]
    f = d["f"] if f is None else np.asarray(f)
    if m.get("stack", {}).get("file"):
        import disp_tmm
        R, T = disp_tmm.lattice_rt(outdir, f)
        return dict(f=f, R=R, T=T)
    Kt = math.hypot(m["Ktilde"][0], m["Ktilde"][2])
    h = np.array(u["hy"])
    eps = np.array(u["eps_n"])
    j_lo, j_hi = m["j0"], m["Ny"] - m["npml_hi"] - 1
    R, T = [], []
    for fi in f:
        w = 2 * math.pi * fi
        wt = 2.0 / m["dt"] * math.sin(w * m["dt"] / 2.0)
        if wt * wt * eps[j_lo - 1] <= Kt * Kt:   # incident wave evanescent at this frequency
            R.append(np.nan); T.append(np.nan); continue
        r = tmm.tmm_discrete(h, eps, m["dt"], Kt, pol=m["pol"], j_lo=j_lo, j_hi=j_hi, omega=w,
                             eps_t=np.array(u["eps_t"]))
        R.append(r["R"]); T.append(r["T"])
    return dict(f=f, R=np.array(R), T=np.array(T))


def angle_deg(meta, f):
    """Continuous incidence angle at frequency f for the fixed kx, kz: sin(theta) = |k_t| / (2 pi f)."""
    kt = math.hypot(meta["kx"], meta["kz"])
    s = kt / (2 * math.pi * np.asarray(f))
    return np.degrees(np.arcsin(np.clip(s, 0, 1)))


if __name__ == "__main__":
    import sys
    me = measure(sys.argv[1])
    th = theory(sys.argv[1])
    meta = load(sys.argv[1])["meta"]
    print(" f       theta    R_meas        R_theory      |dR|      T_meas        T_theory      |dT|     R+T-1")
    for q, fi in enumerate(me["f"]):
        T = me.get("T", np.full_like(me["R"], np.nan))[q]
        print(f"{fi:6.3f} {float(angle_deg(meta, fi)):7.2f}  {me['R'][q]:.6e}  {th['R'][q]:.6e}  "
              f"{abs(me['R'][q] - th['R'][q]):.1e}  {T:.6e}  {th['T'][q]:.6e}  {abs(T - th['T'][q]):.1e}  {me['R'][q] + T - 1:+.1e}")
