"""Analysis library for the nonuniform-grid gates (numpy only). Coordinates come from grid_used.json (never idx*Delta).

Field arrays have the solver dump layout F[c][i, j, k] (Nx, Ny+1, Nz); component positions as in
docs/derivation_nonuniform.md §0. x and z are periodic, y has PEC walls at j = 0, Ny (tangential E = 0).
"""
# ==== nonuniform (whole module is scanned by tests/nu_lint_coords.py)
import math

import numpy as np

import tmm

COMPS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")


# ----------------------------------------------------------------------------- weights and operators (§2)
def spacings(used):
    return dict(hx=used["hx"], dx=used["dx"], hy=used["hy"], dy=used["dy"], hz=used["hz"], dz=used["dz"])


def weights(used):
    """W per component, broadcastable to (Nx, Ny+1, Nz); zero on unused half-grid rows (j = Ny)."""
    s = spacings(used)
    Ny = len(s["hy"])
    hy = np.concatenate([s["hy"], [0.0]])      # half-y rows: last row unused
    X = lambda a: a[:, None, None]  # noqa: E731
    Y = lambda a: a[None, :, None]  # noqa: E731
    Z = lambda a: a[None, None, :]  # noqa: E731
    W = {
        "Ex": X(s["hx"]) * Y(s["dy"]) * Z(s["dz"]),
        "Ey": X(s["dx"]) * Y(hy) * Z(s["dz"]),
        "Ez": X(s["dx"]) * Y(s["dy"]) * Z(s["hz"]),
        "Hx": X(s["dx"]) * Y(hy) * Z(s["hz"]),
        "Hy": X(s["hx"]) * Y(s["dy"]) * Z(s["hz"]),
        "Hz": X(s["hx"]) * Y(hy) * Z(s["dz"]),
    }
    assert W["Ey"].shape[1] == Ny + 1
    return W


def fwd(a, ax):  # a[m+1] - a[m], periodic
    return np.roll(a, -1, axis=ax) - a


def bwd(a, ax):  # a[m] - a[m-1], periodic
    return a - np.roll(a, 1, axis=ax)


def curl_E(F, used):
    """(curl E) at the H positions; mu dH/dt = -curl E."""
    s = spacings(used)
    X = lambda a: a[:, None, None]  # noqa: E731
    Z = lambda a: a[None, None, :]  # noqa: E731
    Ny = len(s["hy"])
    hy = np.concatenate([s["hy"], [1.0]])[None, :, None]
    dyEz = np.zeros_like(F["Ez"])
    dyEz[:, :-1, :] = F["Ez"][:, 1:, :] - F["Ez"][:, :-1, :]
    dyEx = np.zeros_like(F["Ex"])
    dyEx[:, :-1, :] = F["Ex"][:, 1:, :] - F["Ex"][:, :-1, :]
    Hx = dyEz / hy - fwd(F["Ey"], 2) / Z(s["hz"])
    Hy = fwd(F["Ex"], 2) / Z(s["hz"]) - fwd(F["Ez"], 0) / X(s["hx"])
    Hz = fwd(F["Ey"], 0) / X(s["hx"]) - dyEx / hy
    Hx[:, Ny, :] = 0.0
    Hz[:, Ny, :] = 0.0
    return dict(Hx=Hx, Hy=Hy, Hz=Hz)


def curl_H(F, used):
    """(curl H) at the E positions, projected on the PEC space (tangential E rows j = 0, Ny set to 0)."""
    s = spacings(used)
    X = lambda a: a[:, None, None]  # noqa: E731
    Z = lambda a: a[None, None, :]  # noqa: E731
    Y = lambda a: a[None, :, None]  # noqa: E731
    Ny = len(s["hy"])
    dyHz = np.zeros_like(F["Hz"])
    dyHz[:, 1:, :] = F["Hz"][:, 1:, :] - F["Hz"][:, :-1, :]
    dyHx = np.zeros_like(F["Hx"])
    dyHx[:, 1:, :] = F["Hx"][:, 1:, :] - F["Hx"][:, :-1, :]
    Ex = dyHz / Y(s["dy"]) - bwd(F["Hy"], 2) / Z(s["dz"])
    Ey = bwd(F["Hx"], 2) / Z(s["dz"]) - bwd(F["Hz"], 0) / X(s["dx"])
    Ez = bwd(F["Hy"], 0) / X(s["dx"]) - dyHx / Y(s["dy"])
    for a in (Ex, Ez):
        a[:, 0, :] = 0.0
        a[:, Ny, :] = 0.0
    Ey[:, Ny, :] = 0.0
    return dict(Ex=Ex, Ey=Ey, Ez=Ez)


def inner(A, B, W, comps, eps=None):
    s = 0.0
    for c in comps:
        w = W[c] if eps is None or c not in eps else W[c] * eps[c]
        s += float(np.sum(w * A[c] * B[c]))
    return s


def eps_fields(used):
    et, en = used["eps_t"], np.concatenate([used["eps_n"], [used["eps_n"][-1]]])
    return {"Ex": et[None, :, None], "Ez": et[None, :, None], "Ey": en[None, :, None]}


def random_state(used, seed=1):
    rng = np.random.default_rng(seed)
    Nx, Nz, Ny = len(used["hx"]), len(used["hz"]), len(used["hy"])
    F = {c: rng.uniform(-1, 1, (Nx, Ny + 1, Nz)) for c in COMPS}
    for c in ("Ex", "Ez"):
        F[c][:, 0, :] = 0.0
        F[c][:, Ny, :] = 0.0
    for c in ("Ey", "Hx", "Hz"):
        F[c][:, Ny, :] = 0.0
    return F


# ----------------------------------------------------------------------------- divergence (§2)
def div_E(F, used, uniform_delta=None):
    """div E at primal nodes (Nx, Ny+1, Nz); rows j = 0 and j = Ny are not meaningful (set to 0).
    uniform_delta: the (wrong) uniform operator dividing every difference by one spacing (control)."""
    s = spacings(used)
    ddx = bwd(F["Ex"], 0)
    ddz = bwd(F["Ez"], 2)
    ddy = np.zeros_like(F["Ey"])
    ddy[:, 1:, :] = F["Ey"][:, 1:, :] - F["Ey"][:, :-1, :]
    if uniform_delta is not None:
        d = (ddx + ddy + ddz) / uniform_delta
    else:
        d = ddx / s["dx"][:, None, None] + ddy / s["dy"][None, :, None] + ddz / s["dz"][None, None, :]
    d[:, 0, :] = 0.0
    d[:, -1, :] = 0.0
    return d


def div_H(F, used, uniform_delta=None):
    """div H at cell centres (Nx, Ny, Nz)."""
    s = spacings(used)
    a = fwd(F["Hx"], 0)[:, :-1, :]
    c = fwd(F["Hz"], 2)[:, :-1, :]
    b = F["Hy"][:, 1:, :] - F["Hy"][:, :-1, :]
    if uniform_delta is not None:
        return (a + b + c) / uniform_delta
    return a / s["hx"][:, None, None] + b / s["hy"][None, :, None] + c / s["hz"][None, None, :]


# ----------------------------------------------------------------------------- phasor reductions
def transverse_basis(Kt_vec):
    """t = K~t/|K~t| and s = t x y-hat (unit vectors in the x-z plane, as (x, z) pairs)."""
    kx, kz = Kt_vec
    n = math.hypot(kx, kz)
    t = np.array([kx, kz]) / n
    s = np.array([-t[1], t[0]])   # K~ x y-hat has (x, z) components proportional to (-Kz, Kx)
    return t, s


def e_ts(Ex, Ez):
    """Given reduced phasors e_x(y), e_z(y) (transverse phase removed), return dict with e_t, e_s."""
    return Ex, Ez


def reduce_xz(F2d, phase):
    """Least-squares amplitude a in F(x, z) ~ a exp(i phi(x, z)) over a y-plane (Nx, Nz) and the residual."""
    w = F2d * np.conj(phase)
    a = w.mean()
    return a, F2d - a * phase


def two_wave_fit(u, y, ky):
    """u_j = a exp(i ky y_j) + b exp(-i ky y_j): least squares (a, b); returns a, b, relative residual."""
    A = np.stack([np.exp(1j * ky * y), np.exp(-1j * ky * y)], axis=1)
    sol, *_ = np.linalg.lstsq(A, u, rcond=None)
    res = np.linalg.norm(A @ sol - u) / np.linalg.norm(u)
    return complex(sol[0]), complex(sol[1]), float(res)


def ky_three_point(u, h):
    return tmm.ky_three_point(u, h)


def local_ky_forward(u, y, dt, Kt):
    return tmm.local_ky_forward(u, y, dt, Kt)


def plane_phase_fit(F2d, x, z):
    """Fit unwrap-free phase: F ~ A exp(i(kx x + kz z + c)); returns (kx, kz, rms phase residual [rad]).
    Uses the phase of F / F[0, 0] relative to the best plane through the complex ratios of neighbours."""
    ph = np.angle(F2d / F2d[0, 0])
    # unwrap along x then z
    ph = np.unwrap(np.unwrap(ph, axis=0), axis=1)
    X, Z = np.meshgrid(x, z, indexing="ij")
    A = np.stack([X.ravel(), Z.ravel(), np.ones(X.size)], axis=1)
    sol, *_ = np.linalg.lstsq(A, ph.ravel(), rcond=None)
    r = ph.ravel() - A @ sol
    return float(sol[0]), float(sol[1]), float(np.sqrt(np.mean(r ** 2)))


# ----------------------------------------------------------------------------- flux
def sy_plane_weighted(D, used, j, comps_E=("Ex", "Ez")):
    """Time-averaged S_y on the primal plane y_j from phasors of the y-plane DFTs of rows j-1/2 and j+1/2 (H)
    and j (E). H is interpolated to y_j with distance weights from the true node coordinates:
    H(y_j) = [H_{j-1/2} (y_{j+1/2} - y_j) + H_{j+1/2} (y_j - y_{j-1/2})] / (y_{j+1/2} - y_{j-1/2}).
    D: dict with keys (j) -> phasor dict for plane j (E valid) and (j-1), (j) for H rows. Returns mean S_y."""
    yd = used["y_dual"]
    y = used["y"]
    wa = (yd[j] - y[j]) / (yd[j] - yd[j - 1])
    wb = (y[j] - yd[j - 1]) / (yd[j] - yd[j - 1])
    Fm, F0 = D[j - 1], D[j]
    Hx = wa * Fm["Hx"] + wb * F0["Hx"]
    Hz = wa * Fm["Hz"] + wb * F0["Hz"]
    return 0.5 * float(np.real(np.mean(F0["Ez"] * np.conj(Hx) - F0["Ex"] * np.conj(Hz))))


def conserved_flux_rows(Ex, Ez, Hx_up, Hz_up):
    """Discrete-conserved flux through the dual plane between E row j and H row j+1/2 (per unit area)."""
    return 0.5 * float(np.real(np.mean(Ez * np.conj(Hx_up) - Ex * np.conj(Hz_up))))


# ----------------------------------------------------------------------------- Floquet decomposition (stage B)
def floquet(F2d, xs, wx, zs, wz, kxs, kz, Lx, Lz):
    """Coefficients c_p of F(x, z) = sum_p c_p exp(i(kx_p x + kz z)) on nonuniform nodes (quadrature weights wx, wz,
    i.e. the dual spacing for primal samples and the primal spacing for dual samples)."""
    ph_z = np.exp(-1j * kz * zs) * wz
    g = F2d @ ph_z                                     # (Nx,)
    return np.array([np.sum(g * wx * np.exp(-1j * kx * xs)) for kx in kxs]) / (Lx * Lz)


def plane_positions(used, comp):
    """(x positions, x weights, z positions, z weights) of a component's samples on a y-plane slice."""
    px = comp in ("Ey", "Ez", "Hx")                    # primal x
    pz = comp in ("Ex", "Ey", "Hz")                    # primal z
    xs = used["x"][:-1] if px else used["x_dual"]
    wx = used["dx"] if px else used["hx"]
    zs = used["z"][:-1] if pz else used["z_dual"]
    wz = used["dz"] if pz else used["hz"]
    return xs, wx, zs, wz


def plane_flux(D, used):
    """z/x-averaged S_y of a y-plane DFT slice from the real-space samples (E at y_j, H at y_{j+1/2}, quadrature
    weights): the discretely conserved flux, independent of any Floquet partition."""
    Lx, Lz = used["x"][-1], used["z"][-1]
    s = 0.0
    for e, h, sg in (("Ez", "Hx", 1), ("Ex", "Hz", -1)):
        xs, wx, zs, wz = plane_positions(used, e)
        s += sg * 0.5 * np.real(np.sum(D[e] * np.conj(D[h]) * wx[:, None] * wz[None, :]))
    return s / (Lx * Lz)


def floquet_gram_offdiag(used, orders):
    """max |G_pq|, p != q, of exp(i 2π p x / Lx) under the x quadrature weights (primal and dual samples); 0 on a
    uniform grid."""
    Lx = used["x"][-1]
    P = np.asarray(orders)
    m = 0.0
    for xs, w in ((np.asarray(used["x"][:-1]), np.asarray(used["dx"])), (np.asarray(used["x_dual"]), np.asarray(used["hx"]))):
        G = np.exp(2j * np.pi * (P[None, :] - P[:, None])[:, :, None] * xs[None, None, :] / Lx) @ w / Lx
        m = max(m, float(np.max(np.abs(G - np.diag(np.diag(G))))))
    return m


def order_fluxes(D, used, kxs, kz):
    """Per-order z-averaged S_y of a y-plane DFT slice (E at y_j, H at y_{j+1/2}: the discrete-conserved pairing)."""
    Lx, Lz = used["x"][-1], used["z"][-1]
    c = {}
    for comp in ("Ex", "Ez", "Hx", "Hz"):
        xs, wx, zs, wz = plane_positions(used, comp)
        c[comp] = floquet(D[comp], xs, wx, zs, wz, kxs, kz, Lx, Lz)
    return 0.5 * np.real(c["Ez"] * np.conj(c["Hx"]) - c["Ex"] * np.conj(c["Hz"])), c
