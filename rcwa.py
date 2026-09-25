"""Conical RCWA for a 1D lamellar grating (numpy only) -- reference for stage B (SPEC_nonuniform §20.5).

Geometry in the solver's frame: grating period along x (period Lambda), invariant along z, layers stacked along y
(superstrate y < 0, grating layer 0 <= y <= d, substrate y > d). Internally the RCWA frame is (xr, yr, zr) =
(x, −z, y) (right-handed, zr = layer normal). Time convention exp(−iωt); fields ∝ exp(i(kx_p x + kz z)).

Fourier factorization (Li 1996/97): Ex (normal to the ridge walls) uses the inverse rule [eps Ex] = A^{-1} Ex with
A = Toeplitz(1/eps); Ey_r and Ez_r (tangential to the walls) use Laurent's rule E = Toeplitz(eps).
First-order system for psi = [Ex, Ey, Hx, Hy]_r (H scaled by Z0, lengths by 1/k0):
    psi_E' = i P psi_H,  psi_H' = i Q psi_E,
    P = [[Kx E^-1 Ky, I − Kx E^-1 Kx], [Ky E^-1 Ky − I, −Ky E^-1 Kx]],
    Q = [[−Kx Ky, Kx² − E], [A^-1 − Ky², Ky Kx]].
Modes: P Q W = W q², V = Q W q^-1 (Im q >= 0). Boundary matching is solved directly with exponentials referenced
to the near face of the layer (unconditionally stable for one layer). Efficiencies come from the z-flux of each
order, so they do not depend on the (degenerate) eigen-basis in the homogeneous half-spaces.
"""
import math

import numpy as np


def toeplitz_fourier(eps_ridge, eps_groove, duty, x0_frac, M):
    """Toeplitz matrices of eps and 1/eps for a lamellar profile: ridge on x/Lambda in [x0, x0 + duty)."""
    n = np.arange(-2 * M, 2 * M + 1)

    def coeffs(a_ridge, a_groove):
        c = np.empty(len(n), complex)
        for q, m in enumerate(n):
            if m == 0:
                c[q] = a_ridge * duty + a_groove * (1 - duty)
            else:
                # (1/Λ)∫ over the ridge of exp(−2πi m x/Λ) dx, ridge on [x0, x0+duty)Λ
                v = (np.exp(-2j * math.pi * m * (x0_frac + duty)) - np.exp(-2j * math.pi * m * x0_frac)) / (-2j * math.pi * m)
                c[q] = (a_ridge - a_groove) * v
        return c
    ce = coeffs(eps_ridge, eps_groove)
    ca = coeffs(1.0 / eps_ridge, 1.0 / eps_groove)
    N = 2 * M + 1
    idx = np.arange(N)
    D = idx[:, None] - idx[None, :] + 2 * M
    return ce[D], ca[D]


def solve(M, Lambda, d, eps_sup, eps_sub, eps_ridge, eps_groove, duty, kx0, kz0, e_inc, k0=2 * math.pi, x0_frac=0.0):
    """Returns dict(R_orders, T_orders, orders, R, T). e_inc: incident E unit vector in the solver frame (x, y, z)."""
    N = 2 * M + 1
    p = np.arange(-M, M + 1)
    kxp = (kx0 + 2 * math.pi * p / Lambda) / k0
    ky = -kz0 / k0                      # RCWA yr = −z
    Kx = np.diag(kxp.astype(complex))
    Ky = np.diag(np.full(N, ky, complex))
    I = np.eye(N, dtype=complex)

    def PQ(E, Ainv):
        Ei = np.linalg.inv(E)
        P = np.block([[Kx @ Ei @ Ky, I - Kx @ Ei @ Kx], [Ky @ Ei @ Ky - I, -Ky @ Ei @ Kx]])
        Q = np.block([[-Kx @ Ky, Kx @ Kx - E], [Ainv - Ky @ Ky, Ky @ Kx]])
        return P, Q

    def homogeneous(eps):
        qz = np.sqrt(eps - kxp ** 2 - ky ** 2 + 0j)
        qz = np.where(qz.imag < 0, -qz, qz)
        q = np.concatenate([qz, qz])
        P, Q = PQ(eps * I, eps * I)
        W = np.eye(2 * N, dtype=complex)
        V = Q @ W / q[None, :]
        return W, V, q

    W1, V1, q1 = homogeneous(eps_sup)
    W3, V3, q3 = homogeneous(eps_sub)
    E, A = toeplitz_fourier(eps_ridge, eps_groove, duty, x0_frac, M)
    P, Q = PQ(E, np.linalg.inv(A))
    lam, W = np.linalg.eig(P @ Q)
    q = np.sqrt(lam + 0j)
    q = np.where(q.imag < 0, -q, q)
    V = Q @ W / q[None, :]
    X = np.diag(np.exp(1j * q * k0 * d))
    # incident tangential fields (RCWA frame: Ex_r = Ex, Ey_r = −Ez), order 0 only
    einc = np.zeros(2 * N, complex)
    einc[M] = e_inc[0]
    einc[N + M] = -e_inc[2]
    c_inc = np.linalg.solve(W1, einc)
    n2 = 2 * N
    Z = np.zeros((n2, n2), complex)
    # unknowns [r, a, b, t]
    A_ = np.block([[-W1, W, W @ X, Z],
                   [V1, V, -V @ X, Z],
                   [Z, W @ X, W, -W3],
                   [Z, V @ X, -V, -V3]])
    rhs = np.concatenate([W1 @ c_inc, V1 @ c_inc, np.zeros(n2), np.zeros(n2)])
    sol = np.linalg.solve(A_, rhs)
    r, t = sol[:n2], sol[3 * n2:]

    def flux(Wm, Vm, c, sign):
        e = Wm @ c
        h = sign * (Vm @ c)
        ex, ey = e[:N], e[N:]
        hx, hy = h[:N], h[N:]
        return 0.5 * np.real(ex * np.conj(hy) - ey * np.conj(hx))   # z-flux per order (H scaled by Z0)
    S_inc = flux(W1, V1, c_inc, +1).sum()
    R_p = -flux(W1, V1, r, -1) / S_inc
    T_p = flux(W3, V3, t, +1) / S_inc
    return dict(orders=p, R_orders=np.real(R_p), T_orders=np.real(T_p), R=float(np.sum(R_p)), T=float(np.sum(T_p)),
                q_sup=q1[:N], q_sub=q3[:N])


def propagating(res, which="R", tol=1e-12):
    q = res["q_sup"] if which == "R" else res["q_sub"]
    return np.abs(q.imag) < tol
