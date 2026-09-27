"""Reference solutions for layered media (numpy only).

tmm_continuous : exact continuous-Maxwell transfer-matrix solution (the physical reference, "R_TMM").
tmm_discrete   : exact solution of the *discrete* Yee equations reduced to 1D for a given (K~x, K~z)
                 on an arbitrary y-node array ("R_disc", docs/derivation_nonuniform.md §5).

Units: c = eps0 = mu0 = 1, lambda0 = 1, k0 = omega0 = 2*pi.
"""
import math

import numpy as np


# ----------------------------------------------------------------------------- continuous TMM
def tmm_continuous(layers, eps_in, eps_out, kt, k0=2 * math.pi, pol="s"):
    """layers: list of (eps, thickness) from the incidence side. Returns dict(r, t, R, T).

    s: E tangential (admittance Y = q); p: H tangential (admittance Y = eps/q), q = sqrt(eps k0^2 - kt^2).
    """
    def q_of(eps):
        return np.sqrt(complex(eps * k0 * k0 - kt * kt))

    def Y(eps):
        q = q_of(eps)
        return q if pol == "s" else eps / q

    M = np.eye(2, dtype=complex)
    for eps, d in layers:
        q, y = q_of(eps), Y(eps)
        c, s = np.cos(q * d), np.sin(q * d)
        M = M @ np.array([[c, -1j * s / y], [-1j * y * s, c]])
    y0, ys = Y(eps_in), Y(eps_out)
    den = y0 * M[0, 0] + y0 * ys * M[0, 1] + M[1, 0] + ys * M[1, 1]
    r = (y0 * M[0, 0] + y0 * ys * M[0, 1] - M[1, 0] - ys * M[1, 1]) / den
    t = 2 * y0 / den
    R = abs(r) ** 2
    T = float((ys.real / y0.real) * abs(t) ** 2)
    return dict(r=complex(r), t=complex(t), R=float(R), T=T)


def fresnel_numeric_estimate(ky1, ky2, eps1=1.0, eps2=1.0, pol="s"):
    """Fresnel reflectance with 'numerical indices' (phase-velocity ky); INFO-only estimate (§5)."""
    if pol == "s":
        r = (ky1 - ky2) / (ky1 + ky2)
    else:
        r = (ky1 / eps1 - ky2 / eps2) / (ky1 / eps1 + ky2 / eps2)
    return float(abs(r) ** 2)


# ----------------------------------------------------------------------------- discrete theory
def omega_tilde(dt, omega=2 * math.pi):
    return 2.0 / dt * math.sin(omega * dt / 2.0)


def ktilde(k, h):
    return 2.0 / h * math.sin(k * h / 2.0)


def ky_local(h, eps, dt, Kt, omega=2 * math.pi, mu=1.0):
    """Discrete ky in a uniform sub-region of spacing h: (2/h)^2 sin^2(ky h/2) = w~^2 eps mu - Kt^2."""
    q = (h / 2.0) ** 2 * (omega_tilde(dt, omega) ** 2 * eps * mu - Kt * Kt)
    if not (0.0 < q < 1.0):
        raise ValueError(f"no propagating discrete ky for h={h}, eps={eps} (sin^2 = {q})")
    return 2.0 / h * math.asin(math.sqrt(q))


def ky_continuous(kt, eps=1.0, k0=2 * math.pi):
    return math.sqrt(eps * k0 * k0 - kt * kt)


def theta_local(Kt, ky):
    return math.atan2(abs(Kt), ky)


def eps_nodes(h, eps_cell):
    """eps_t at primal nodes (h-weighted average of the two adjacent cells; end nodes take their one cell)."""
    h = np.asarray(h, float)
    e = np.asarray(eps_cell, float)
    et = np.empty(len(h) + 1)
    et[1:-1] = (h[:-1] * e[:-1] + h[1:] * e[1:]) / (h[:-1] + h[1:])
    et[0], et[-1] = e[0], e[-1]
    return et


def _thomas(a, b, c, r):
    """Solve a[m] u[m-1] + b[m] u[m] + c[m] u[m+1] = r[m] (a[0], c[-1] ignored)."""
    n = len(b)
    cp = np.empty(n, complex)
    rp = np.empty(n, complex)
    cp[0] = c[0] / b[0]
    rp[0] = r[0] / b[0]
    for m in range(1, n):
        den = b[m] - a[m] * cp[m - 1]
        cp[m] = c[m] / den if m < n - 1 else 0.0
        rp[m] = (r[m] - a[m] * rp[m - 1]) / den
    u = np.empty(n, complex)
    u[-1] = rp[-1]
    for m in range(n - 2, -1, -1):
        u[m] = rp[m] - cp[m] * u[m + 1]
    return u


def three_point_scatter(w, g, q, gL, gR, wL, wR, qL, qR):
    """Scattering solution of (1/w_m)[(u_{m+1}-u_m)/g_{m+1/2} - (u_m-u_{m-1})/g_{m-1/2}] + q_m u_m = 0,
    m = 0..M, with uniform half-lines on both sides (left: w=wL, g=gL, q=qL; right: wR, gR, qR).
    g has length M (links between interior nodes); the outer links use gL, gR.
    Incident amplitude 1 from the left, reference phase at node 0. Returns (u, r, t_mod, R, T)."""
    w, g, q = (np.asarray(v, float) for v in (w, g, q))
    M = len(w) - 1
    cL = 1.0 - qL * wL * gL / 2.0
    cR = 1.0 - qR * wR * gR / 2.0
    if not (-1 < cL < 1 and -1 < cR < 1):
        raise ValueError("end regions are not propagating")
    phiL, phiR = math.acos(cL), math.acos(cR)
    gm = np.concatenate([[gL], g])            # g_{m-1/2}
    gp = np.concatenate([g, [gR]])            # g_{m+1/2}
    a = 1.0 / (w * gm)
    c = 1.0 / (w * gp)
    b = (-(a + c) + q).astype(complex)
    a = a.astype(complex)
    c = c.astype(complex)
    rhs = np.zeros(M + 1, complex)
    eL = np.exp(1j * phiL)
    b[0] += a[0] * eL                          # u_{-1} = u_0 e^{i phiL} - 2i sin(phiL)
    rhs[0] = a[0] * 2j * math.sin(phiL)
    b[M] += c[M] * np.exp(1j * phiR)           # u_{M+1} = u_M e^{i phiR}
    u = _thomas(a, b, c, rhs)
    r = u[0] - 1.0
    R = abs(r) ** 2
    T = abs(u[M]) ** 2 * (math.sin(phiR) / gR) / (math.sin(phiL) / gL)
    return u, complex(r), float(abs(u[M])), float(R), float(T)


def tmm_discrete(h, eps_cell, dt, Kt, pol="s", j_lo=None, j_hi=None, omega=2 * math.pi, eps_t=None):
    """Exact Yee-lattice reflection/transmission for the y-grid (primal spacings h[j] between y_j and y_{j+1},
    eps_cell[j] in that cell = eps_n at the dual node). eps_t (primal nodes) defaults to the h-weighted average.
    The solver region is nodes j_lo..j_hi; the grid must be uniform (spacing and eps) to the left of j_lo and to
    the right of j_hi (those half-lines are closed with exact discrete modes).

    s: unknown e_s at primal nodes; p: unknown h_s at dual nodes (§5 table).
    Returns dict(R, T, r, u, nodes) where nodes are the y-indices (s: j, p: j+1/2 as j) of u."""
    h = np.asarray(h, float)
    e = np.asarray(eps_cell, float)
    N = len(h)
    j_lo = 1 if j_lo is None else j_lo
    j_hi = N - 1 if j_hi is None else j_hi
    w2 = omega_tilde(dt, omega) ** 2
    et = eps_nodes(h, e) if eps_t is None else np.asarray(eps_t, float)
    if pol == "s":
        js = np.arange(j_lo, j_hi + 1)
        d = 0.5 * (h[js - 1] + h[js])
        w = d
        g = h[j_lo:j_hi]
        q = w2 * et[js] - Kt * Kt
        hL, hR = h[j_lo - 1], h[j_hi]
        eL, eR = e[j_lo - 1], e[j_hi]
        res = three_point_scatter(w, g, q, gL=hL, gR=hR, wL=hL, wR=hR, qL=w2 * eL - Kt * Kt, qR=w2 * eR - Kt * Kt)
    else:
        js = np.arange(j_lo, j_hi)            # dual nodes j+1/2
        w = h[js]
        dnext = 0.5 * (h[js[:-1]] + h[js[:-1] + 1])
        g = et[js[:-1] + 1] * dnext
        q = w2 - Kt * Kt / e[js]
        hL, hR = h[j_lo - 1], h[j_hi]
        eL, eR = e[j_lo - 1], e[j_hi]
        res = three_point_scatter(w, g, q, gL=eL * hL, gR=eR * hR, wL=hL, wR=hR,
                                  qL=w2 - Kt * Kt / eL, qR=w2 - Kt * Kt / eR)
    u, r, tm, R, T = res
    return dict(R=R, T=T, r=r, u=u, nodes=js)


# ----------------------------------------------------------------------------- estimators
def ky_three_point(u, h):
    """Least-squares cos(ky h) from u_{j+1} + u_{j-1} = 2 cos(ky h) u_j on a uniform node set.
    Returns (ky, relative imaginary residual)."""
    u = np.asarray(u, complex)
    num = np.vdot(u[1:-1], u[2:] + u[:-2])
    den = 2 * np.vdot(u[1:-1], u[1:-1])
    c = num / den
    return math.acos(float(c.real)) / h, float(abs(c.imag) / abs(c))


def forward_amplitude(u, y, dt, Kt, eps=1.0, omega=2 * math.pi):
    """Local forward-wave amplitude A_{j+1/2} at the midpoints of consecutive samples u_j, u_{j+1} located at y_j
    (derivation §5): A = 1/2[(u_j+u_{j+1})/(2 cos(k h/2)) + (u_{j+1}-u_j)/(i K~y h)], k = ky_local(h)."""
    u = np.asarray(u, complex)
    y = np.asarray(y, float)
    h = np.diff(y)
    Ky = math.sqrt(omega_tilde(dt, omega) ** 2 * eps - Kt * Kt)
    k = np.array([ky_local(hh, eps, dt, Kt, omega) for hh in h])
    A = 0.5 * ((u[:-1] + u[1:]) / (2 * np.cos(k * h / 2)) + (u[1:] - u[:-1]) / (1j * Ky * h))
    return A, 0.5 * (y[:-1] + y[1:])


def local_ky_forward(u, y, dt, Kt, eps=1.0):
    """ky_meas at the interior samples j (between midpoints j-1/2, j+1/2): arg(A_{j+1/2}/A_{j-1/2}) / d_j.
    Returns (y_j, d_j, ky_meas)."""
    A, ym = forward_amplitude(u, y, dt, Kt, eps)
    d = np.diff(ym)
    ky = np.angle(A[1:] / A[:-1]) / d
    return np.asarray(y)[1:-1], d, ky
