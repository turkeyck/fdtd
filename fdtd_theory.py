"""Discrete-plane-wave theory for the 3D Yee FDTD oblique-incidence project.

Units are normalized: c = eps0 = mu0 = 1 (so eta0 = 1) and lambda0 = 1.
Formulas and sign conventions follow derivation.md.
"""
import math
import numpy as np

C0 = 1.0
EPS0 = 1.0
MU0 = 1.0
ETA0 = math.sqrt(MU0 / EPS0)

# Yee offsets in units of Delta: (x, y, z) and time offset in units of dt.
YEE = {
    "Ex": ((0.5, 0.0, 0.0), 0.0),
    "Ey": ((0.0, 0.5, 0.0), 0.0),
    "Ez": ((0.0, 0.0, 0.5), 0.0),
    "Hx": ((0.0, 0.5, 0.5), 0.5),
    "Hy": ((0.5, 0.0, 0.5), 0.5),
    "Hz": ((0.5, 0.5, 0.0), 0.5),
}
ECOMP = ("Ex", "Ey", "Ez")
HCOMP = ("Hx", "Hy", "Hz")
COMPS = ECOMP + HCOMP


class PropagationError(ValueError):
    pass


class Setup:
    """Grid + mode parameters. n_ref is the refractive index of the medium the wave lives in."""

    def __init__(self, n_lambda=20, S=0.5, Lx=2.0, Lz=3.0, m=1, n=1, pol="s",
                 lam0=1.0, n_ref=1.0, margin=0.05):
        self.n_lambda = n_lambda
        self.lam0 = lam0
        self.Delta = lam0 / n_lambda
        self.S = S
        self.dt = S * self.Delta / C0
        self.Lx, self.Lz = Lx, Lz
        self.Nx = int(round(Lx / self.Delta))
        self.Nz = int(round(Lz / self.Delta))
        if abs(self.Nx * self.Delta - Lx) > 1e-12 * Lx or abs(self.Nz * self.Delta - Lz) > 1e-12 * Lz:
            raise ValueError("Lx and Lz must be integer multiples of Delta")
        if S >= 1.0 / math.sqrt(3.0):
            raise ValueError("Courant number must be < 1/sqrt(3)")
        self.m, self.n, self.pol = m, n, pol
        self.n_ref = n_ref
        self.omega = 2.0 * math.pi * C0 / lam0
        self.k0 = self.omega / C0
        self.kx = 2.0 * math.pi * m / Lx
        self.kz = 2.0 * math.pi * n / Lz
        kt = math.hypot(self.kx, self.kz)
        if kt == 0.0:
            raise PropagationError("normal incidence (m=n=0): s/p basis K~ x y-hat is degenerate; unsupported")
        if kt >= (1.0 - margin) * n_ref * self.k0:
            raise PropagationError(f"|k_t|/(n k0) = {kt / (n_ref * self.k0):.4f} exceeds 1-margin={1 - margin}")
        self.ky = ky_discrete(self.omega, self.kx, self.kz, self.Delta, self.dt, n_ref)
        self.ky_cont = ky_continuous(self.omega, self.kx, self.kz, n_ref)

    # ---- discrete quantities -------------------------------------------------
    @property
    def k(self):
        return np.array([self.kx, self.ky, self.kz])

    @property
    def Kt(self):
        return ktilde(self.k, self.Delta)

    @property
    def wt(self):
        return omega_tilde(self.omega, self.dt)

    def amplitudes(self, E0=1.0, continuous=None):
        """Return real amplitude vectors (E0vec, H0vec).

        continuous=None      : exact discrete plane wave (K~ from discrete ky).
        continuous="ky"      : control C1 -- continuous ky, but the same K~/omega~ construction.
        continuous="full"    : control C2 -- textbook continuous wave, E0 perp k, H0 = k x E0/(mu0 omega).
        """
        if continuous is None:
            K = self.Kt
        elif continuous == "ky":
            K = ktilde(np.array([self.kx, self.ky_cont, self.kz]), self.Delta)
        elif continuous == "full":
            K = np.array([self.kx, self.ky_cont, self.kz])
        else:
            raise ValueError(continuous)
        s_hat, p_hat = pol_basis(K)
        e = E0 * (s_hat if self.pol == "s" else p_hat)
        if continuous == "full":
            h = np.cross(K, e) / (MU0 * self.omega)
        else:
            h = np.cross(K, e) / (MU0 * self.wt)
        return e, h

    def phase_k(self, continuous=None):
        if continuous is None:
            return self.k
        return np.array([self.kx, self.ky_cont, self.kz])


def ky_discrete(omega, kx, kz, Delta, dt, n_ref=1.0):
    """Closed-form ky from the 3D Yee dispersion relation (medium index n_ref)."""
    c = C0 / n_ref
    q = (Delta * math.sin(omega * dt / 2.0) / (c * dt)) ** 2 \
        - math.sin(kx * Delta / 2.0) ** 2 - math.sin(kz * Delta / 2.0) ** 2
    if not (0.0 < q <= 1.0):
        raise PropagationError(f"sin^2(ky Delta/2) = {q:.6g} is outside (0, 1]")
    return 2.0 / Delta * math.asin(math.sqrt(q))


def ky_continuous(omega, kx, kz, n_ref=1.0):
    return math.sqrt((n_ref * omega / C0) ** 2 - kx ** 2 - kz ** 2)


def ktilde(k, Delta):
    return 2.0 / Delta * np.sin(np.asarray(k, dtype=float) * Delta / 2.0)


def omega_tilde(omega, dt):
    return 2.0 / dt * math.sin(omega * dt / 2.0)


def pol_basis(K):
    """s_hat ∝ K x y_hat, p_hat ∝ s_hat x K (both unit vectors)."""
    K = np.asarray(K, dtype=float)
    s = np.cross(K, [0.0, 1.0, 0.0])
    s /= np.linalg.norm(s)
    p = np.cross(s, K)
    p /= np.linalg.norm(p)
    return s, p


def yee_coords(comp, I, J, K, Delta):
    (ox, oy, oz), _ = YEE[comp]
    return (I + ox) * Delta, (J + oy) * Delta, (K + oz) * Delta


def analytic(setup, comp, I, J, K, nstep, E0=1.0, continuous=None, phi0=0.0):
    """Real field of the plane wave sampled at comp's own Yee point and time.

    I, J, K may be broadcastable integer arrays; nstep is the integer time index n
    (E sampled at n*dt, H at (n+1/2)*dt).
    """
    e, h = setup.amplitudes(E0, continuous)
    vec = e if comp in ECOMP else h
    a = vec["xyz".index(comp[1])]
    x, y, z = yee_coords(comp, I, J, K, setup.Delta)
    kx, ky, kz = setup.phase_k(continuous)
    t = (nstep + YEE[comp][1]) * setup.dt
    return a * np.cos(kx * x + ky * y + kz * z - setup.omega * t + phi0)


def phasor_amplitude(setup, comp, E0=1.0):
    """Complex phasor amplitude A such that field = Re[A exp(i(k.r - w t))] at comp's Yee point."""
    e, h = setup.amplitudes(E0)
    vec = e if comp in ECOMP else h
    return complex(vec["xyz".index(comp[1])])


def interp_factors(setup):
    """Amplitude factors produced by averaging each component to the cell centre (i+1/2, j+1/2, k+1/2).

    Averaging two samples a half-cell either side of the target multiplies a plane-wave
    phasor by cos(k_d Delta/2) for each averaged direction d.
    """
    cx, cy, cz = (math.cos(kd * setup.Delta / 2.0) for kd in setup.k)
    return {
        "Ex": cy * cz, "Ey": cx * cz, "Ez": cx * cy,
        "Hx": cx, "Hy": cy, "Hz": cz,
    }


def group_velocity(setup):
    """Yee group velocity v_g = grad_k omega on the discrete dispersion shell (medium n_ref)."""
    c = C0 / setup.n_ref
    D, dt, w = setup.Delta, setup.dt, setup.omega
    return c ** 2 * dt * np.sin(setup.k * D) / (D * math.sin(w * dt))


def discrete_fresnel(st_vac, st_med, eps_if):
    """Exact R, T of the Yee lattice for a planar interface on an integer-y plane (derivation.md §8.1).

    Vacuum (y < y1) holds incident + reflected discrete plane waves, the medium (y > y1) the transmitted one;
    the unknowns r, t follow from continuity of tangential E on the interface plane and Ampère's law at that
    plane with eps_if (the value used for Ex, Ez there). Fluxes use the discrete-conserved Phi."""
    D, wt = st_vac.Delta, st_vac.wt
    kv = np.array([st_vac.kx, st_vac.ky, st_vac.kz])
    km = np.array([st_med.kx, st_med.ky, st_med.kz])
    kr = kv * np.array([1.0, -1.0, 1.0])

    def wave(k):
        K = ktilde(k, D)
        s, p = pol_basis(K)
        e = s if st_vac.pol == "s" else p
        return e, np.cross(K, e) / (MU0 * wt)
    ei, hi = wave(kv)
    er, hr = wave(kr)
    et, ht = wave(km)
    K = ktilde(kv, D)
    tan = np.array([1.0, 0.0, 1.0])
    u = ei * tan if st_vac.pol == "s" else np.array([K[0], 0.0, K[2]])
    u = u / np.linalg.norm(u)
    at = lambda k, y: np.exp(1j * k[1] * y)  # noqa: E731

    def residuals(r, t):
        Ev = (ei + r * er) * tan                                  # tangential E on the interface plane (y = 0)
        Em = t * et * tan
        Hdn = hi * at(kv, -D / 2) + r * hr * at(kr, -D / 2)      # vacuum H at y = -D/2
        Hup = t * ht * at(km, D / 2)                             # medium H at y = +D/2
        Hy0 = (K[2] * Ev[0] - K[0] * Ev[2]) / (MU0 * wt)         # Faraday at y = 0 (only tangential E enters)
        rx = 1j * wt * EPS0 * eps_if * Ev[0] + (Hup[2] - Hdn[2]) / D - 1j * K[2] * Hy0
        rz = 1j * wt * EPS0 * eps_if * Ev[2] + 1j * K[0] * Hy0 - (Hup[0] - Hdn[0]) / D
        return np.array([u @ (Ev - Em), u[0] * rx + u[2] * rz])
    f0 = residuals(0.0, 0.0)
    A = np.stack([residuals(1.0, 0.0) - f0, residuals(0.0, 1.0) - f0], axis=1)
    r, t = np.linalg.solve(A, -f0)

    def phi(e, h, k):
        return 0.5 * math.cos(k[1] * D / 2) * np.real(e[2] * np.conj(h[0]) - e[0] * np.conj(h[2]))
    pin = phi(ei, hi, kv)
    return dict(r=complex(r), t=complex(t), R=float(-abs(r) ** 2 * phi(er, hr, kr) / pin),
                T=float(abs(t) ** 2 * phi(et, ht, km) / pin))


def fresnel(n1, n2, cos1, cos2):
    rs = (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    rp = (n2 * cos1 - n1 * cos2) / (n2 * cos1 + n1 * cos2)
    Rs, Rp = rs ** 2, rp ** 2
    Ts = 1.0 - Rs
    Tp = 1.0 - Rp
    return {"rs": rs, "rp": rp, "Rs": Rs, "Rp": Rp, "Ts": Ts, "Tp": Tp}
