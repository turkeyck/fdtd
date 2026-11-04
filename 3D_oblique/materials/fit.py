"""Material dispersion: experimental data (refractiveindex.info, CC0) -> multi-pole Lorentz model for the FDTD ADE.

Model (time convention e^{-i w t}, Im eps > 0 = loss), frequencies as wavenumbers nu = 1/lambda in 1/um:

    eps(nu) = eps_inf + sum_k  de_k nu_k^2 / (nu_k^2 - nu^2 - i g_k nu)

Fitted on complex eps over the band (default 400-800 nm) with Levenberg-Marquardt (numpy only). Positivity of
de_k, nu_k, g_k is enforced by fitting their logarithms, so the model is passive and causal (Kramers-Kronig
consistent), the property Lumerical's Multi-Coefficient Model (MCM) also enforces. Error metric = the one MCM
reports: RMS of |eps_model - eps_data| / |eps_data| over the sampled band.

    python3 materials/fit.py            # fit all, write materials/models.json and figures/disp/materials_fit.png
"""
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ROOT = os.path.dirname(HERE)

# name -> (data file, number of poles, initial poles [(de, nu0 [1/um], g [1/um])], eps_inf0)
SPECS = {
    "SiO2": ("SiO2_nk_Malitson.yml", 1, [(1.1, 9.0, 0.0)], 1.0),
    "Si3N4": ("Si3N4_nk_Luke.yml", 1, [(3.0, 7.0, 0.0)], 1.0),
    "Si": ("Si_nk_Green-2008.yml", 3, [(4.0, 2.75, 0.03), (0.5, 2.7, 0.3), (7.0, 5.4, 0.4)], 1.0),
    "GaAs": ("GaAs_nk_Aspnes.yml", 4, [(0.4, 2.47, 0.23), (0.55, 2.35, 0.17), (5.0, 2.8, 0.6), (1.0, 3.4, 0.5)], 1.0),
}
BAND_NM = (400.0, 800.0)


# ----------------------------------------------------------------------------- data (refractiveindex.info yml)
def read_yml(path):
    """Returns (refs text, evaluator lam_um -> complex n) for 'tabulated nk' or 'formula 1' entries."""
    txt = open(path, encoding="utf-8").read()
    refs = re.search(r"REFERENCES: \|\n(.*?)\n[A-Z]+:", txt, re.S)
    refs = refs.group(1).strip() if refs else ""
    if "type: tabulated nk" in txt:
        block = txt.split("data: |", 1)[1]
        rows = []
        for line in block.splitlines()[1:]:
            p = line.split()
            if len(p) != 3:
                break
            rows.append([float(v) for v in p])
        a = np.array(rows)

        def n_of(lam):
            lam = np.asarray(lam, float)
            if lam.min() < a[0, 0] - 1e-12 or lam.max() > a[-1, 0] + 1e-12:
                raise ValueError(f"{os.path.basename(path)}: outside tabulated range")
            return np.interp(lam, a[:, 0], a[:, 1]) + 1j * np.interp(lam, a[:, 0], a[:, 2])
        return refs, n_of, a
    if "type: formula 1" in txt:
        c = [float(v) for v in re.search(r"coefficients: (.*)", txt).group(1).split()]

        def n_of(lam):
            lam = np.asarray(lam, float)
            n2 = 1.0 + c[0] + sum(c[i] * lam ** 2 / (lam ** 2 - c[i + 1] ** 2) for i in range(1, len(c), 2))
            return np.sqrt(n2).astype(complex)
        return refs, n_of, None
    raise ValueError(f"unsupported entry in {path}")


# ----------------------------------------------------------------------------- model
def eps_model(nu, einf, poles):
    nu = np.asarray(nu, float)
    e = np.full(nu.shape, einf, complex)
    for de, n0, g in poles:
        e += de * n0 ** 2 / (n0 ** 2 - nu ** 2 - 1j * g * nu)
    return e


def _lin(nu, eps, w, n0g):
    """Variable projection: for fixed (nu0_k, g_k) the model is linear in (eps_inf, de_k); weighted LS on the
    real-valued stacked system. Returns (eps_inf, de list, residual vector)."""
    cols = [np.ones_like(nu, dtype=complex)]
    for n0, g in n0g:
        cols.append(n0 ** 2 / (n0 ** 2 - nu ** 2 - 1j * g * nu))
    A = np.stack(cols, axis=1) * w[:, None]
    b = eps * w
    Ar = np.concatenate([A.real, A.imag])
    br = np.concatenate([b.real, b.imag])
    c, *_ = np.linalg.lstsq(Ar, br, rcond=None)
    r = Ar @ c - br
    pen = max(0.0, 1.0 - c[0]) + sum(max(0.0, -ck) for ck in c[1:])   # eps_inf >= 1, de >= 0
    return c[0], list(c[1:]), np.concatenate([r, [1e3 * pen]])


def fit(nu, eps, einf0, poles0, iters=300):
    """Levenberg-Marquardt over log(nu0_k), log(g_k) with the linear parameters projected out (relative weights)."""
    w = 1.0 / np.abs(eps)
    x = np.array([v for _, n0, g in poles0 for v in (math.log(n0), math.log(max(g, 1e-4)))])

    def unp(x):
        return [(math.exp(x[2 * k]), math.exp(x[2 * k + 1])) for k in range(len(x) // 2)]

    def res(x):
        return _lin(nu, eps, w, unp(x))[2]
    lam = 1e-2
    r = res(x)
    for _ in range(iters):
        J = np.empty((len(r), len(x)))
        for q in range(len(x)):
            xp = x.copy(); xp[q] += 1e-7
            J[:, q] = (res(xp) - r) / 1e-7
        A = J.T @ J
        g = J.T @ r
        improved = False
        while lam < 1e12:
            dx = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-14), -g)
            sc = np.abs(dx).max()
            if sc > 0.5:
                dx *= 0.5 / sc
            xn = np.clip(x + dx, -12.0, 6.0)
            r2 = res(xn)
            if r2 @ r2 < r @ r:
                x, r, lam, improved = xn, r2, max(lam / 3, 1e-12), True
                break
            lam *= 4
        if not improved or np.abs(dx).max() < 1e-10:
            break
    n0g = unp(x)
    einf, de, _ = _lin(nu, eps, w, n0g)
    return einf, [(d, n0, g) for d, (n0, g) in zip(de, n0g)]


def fit_multistart(nu, eps, starts):
    best = None
    for einf0, p0 in starts:
        einf, poles = fit(nu, eps, einf0, p0)
        if einf < 1.0 - 1e-6 or any(d < -1e-9 for d, _, _ in poles):
            continue
        einf = max(einf, 1.0)                              # eps_inf >= 1 keeps the vacuum Courant bound valid
        poles = [(d, n0, g) for d, n0, g in poles if d > 1e-6]
        e = np.abs(eps_model(nu, einf, poles) - eps) / np.abs(eps)
        err = float(np.sqrt(np.mean(e ** 2)))
        if best is None or err < best[0]:
            best = (err, einf, poles)
    if best is None:
        raise RuntimeError("no passive fit found")
    return best[1], best[2]


def fit_all(band_nm=BAND_NM, npts=81):
    lam_um = np.linspace(band_nm[0], band_nm[1], npts) / 1000.0
    nu = 1.0 / lam_um
    out = {}
    for name, (fn, npole, p0, e0) in SPECS.items():
        refs, n_of, tab = read_yml(os.path.join(DATA, fn))
        n = n_of(lam_um)
        eps = n * n
        rng = np.random.default_rng(1)
        starts = [(e0, p0[:npole])]
        for _ in range(30):   # random multi-start around the physical guesses (fixed seed: reproducible)
            starts.append((e0, [(de, n0 * math.exp(rng.normal(0, 0.3)), g * math.exp(rng.normal(0, 0.8)) + 1e-3)
                                for de, n0, g in p0[:npole]]))
        einf, poles = fit_multistart(nu, eps, starts)
        em = eps_model(nu, einf, poles)
        nm = np.sqrt(em)
        rel = np.abs(em - eps) / np.abs(eps)
        out[name] = dict(
            source=fn, references=refs, band_nm=list(band_nm), units="nu = 1/lambda in 1/um; eps = eps_inf + "
            "sum de nu0^2/(nu0^2 - nu^2 - i g nu), e^{-i w t}",
            eps_inf=einf, poles=[dict(de=de, nu0=n0, g=g) for de, n0, g in poles],
            fit_rms_rel_eps=float(np.sqrt(np.mean(rel ** 2))), fit_max_rel_eps=float(rel.max()),
            fit_max_abs_dn=float(np.abs(nm.real - n.real).max()), fit_max_abs_dk=float(np.abs(nm.imag - n.imag).max()))
    return out, lam_um


def model_of(entry):
    return entry["eps_inf"], [(p["de"], p["nu0"], p["g"]) for p in entry["poles"]]


def main():
    out, lam = fit_all()
    with open(os.path.join(HERE, "models.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = os.path.join(ROOT, "figures", "disp")
    os.makedirs(figdir, exist_ok=True)
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    lam_f = np.linspace(380, 820, 400) / 1000.0
    for ax, (name, e) in zip(axs.flat, out.items()):
        _, n_of, tab = read_yml(os.path.join(DATA, e["source"]))
        lo = max(lam_f[0], tab[0, 0]) if tab is not None else lam_f[0]
        hi = min(lam_f[-1], tab[-1, 0]) if tab is not None else lam_f[-1]
        ll = lam_f[(lam_f >= lo) & (lam_f <= hi)]
        n = n_of(ll)
        nm = np.sqrt(eps_model(1 / ll, *model_of(e)))
        ax.plot(ll * 1000, n.real, "k.", ms=3, label="n data")
        ax.plot(ll * 1000, nm.real, "C0-", label="n model")
        ax2 = ax.twinx()
        ax2.plot(ll * 1000, n.imag, "k+", ms=3, label="k data")
        ax2.plot(ll * 1000, nm.imag, "C3-", label="k model")
        ax.axvspan(*BAND_NM, color="0.92", zorder=0)
        ax.set_title(f"{name}: {len(e['poles'])} Lorentz poles, RMS rel ε err {e['fit_rms_rel_eps']:.2e}", fontsize=9)
        ax.set_xlabel("λ (nm)"); ax.set_ylabel("n"); ax2.set_ylabel("k")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "materials_fit.png"), dpi=130)
    for name, e in out.items():
        print(f"{name:6s} poles={len(e['poles'])} eps_inf={e['eps_inf']:.4f} RMS rel eps={e['fit_rms_rel_eps']:.2e} "
              f"max rel={e['fit_max_rel_eps']:.2e} max|dn|={e['fit_max_abs_dn']:.2e} max|dk|={e['fit_max_abs_dk']:.2e}")
        for p in e["poles"]:
            print(f"        de={p['de']:.5g} nu0={p['nu0']:.5g}/um ({1000 / p['nu0']:.1f} nm) g={p['g']:.4g}/um")


if __name__ == "__main__":
    main()
