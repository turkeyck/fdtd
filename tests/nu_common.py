"""Shared helpers for the nonuniform-grid gates (tests/nu_*.py): grids, time step, wave numbers, thresholds,
result tables. Coordinates always come from the grid (never idx*Delta)."""
import datetime
import hashlib
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import grid_gen  # noqa: E402
import tmm  # noqa: E402

GRIDS = os.path.join(ROOT, "grids")
RESULTS = os.path.join(ROOT, "results")
FIGS = os.path.join(ROOT, "figures", "nu")
RUNS = os.path.join(ROOT, "runs", "nu")
THRESH_FILE = os.path.join(ROOT, "tests", "thresholds_nu.json")
PRED_FILE = os.path.join(RESULTS, "nu_predictions.json")
K0 = 2 * math.pi
S_DEFAULT = 0.5


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def thresholds():
    with open(THRESH_FILE, encoding="utf-8") as fh:
        return json.load(fh)


def grid_path(name):
    return os.path.join(GRIDS, name + ".json")


def load_grid(name):
    with open(grid_path(name)) as fh:
        g = json.load(fh)
    bad = grid_gen.check(g)
    if bad:
        raise ValueError(f"grid {name}: {bad}")
    return g


class Geo:
    """Arrays of a grid: primal/dual node coordinates and spacings per axis, eps, zones."""

    def __init__(self, g):
        self.g = g
        self.name = g["name"]
        self.hy = np.array(g["y"]["h"])
        self.y = grid_gen.nodes(list(self.hy))
        self.yh = 0.5 * (self.y[:-1] + self.y[1:])
        self.dy = np.empty(len(self.y))
        self.dy[1:-1] = 0.5 * (self.hy[:-1] + self.hy[1:])
        self.dy[0], self.dy[-1] = self.hy[0], self.hy[-1]
        self.hx = np.array(g["x"]["h"])
        self.hz = np.array(g["z"]["h"])
        self.Ny = len(self.hy)
        self.Nx, self.Nz = len(self.hx), len(self.hz)
        self.eps_cell = np.array(g.get("eps_n", g["eps_y"]))
        self.eps_t = np.array(g["eps_t"]) if "eps_t" in g else tmm.eps_nodes(self.hy, self.eps_cell)
        z = g["zones"]
        self.j0, self.ja, self.npml_lo, self.npml_hi = z["j0"], z["ja"], z["npml_lo"], z["npml_hi"]
        self.marks = {k: tuple(v) for k, v in g.get("marks", {}).items()}
        self.Lx, self.Lz = g["x"]["L"], g["z"]["L"]

    def dt(self, S=S_DEFAULT):
        return dt_courant(self.hx.min(), self.hy.min(), self.hz.min(), S)

    def dt_run(self, S=S_DEFAULT):
        """dt actually used for every nonuniform run (implementation note D12): an integer number of steps per
        period, dt = T0 / ceil(T0 / dt_Courant) <= dt_Courant, so a DFT window of whole periods is exact at f0."""
        return 1.0 / math.ceil(1.0 / self.dt(S) - 1e-9)

    def uniform_xz(self):
        return bool(self.g["x"]["uniform"]) and bool(self.g["z"]["uniform"])


def dt_courant(dx, dy, dz, S=S_DEFAULT):
    """dt = S sqrt(3)/sqrt(1/dx^2 + 1/dy^2 + 1/dz^2); exactly S*dx when the three minima are equal
    (same rule as the C code, so the uniform path reproduces dt = S*Delta bit for bit)."""
    if dx == dy == dz:
        return S * dx
    return S * math.sqrt(3.0) / math.sqrt(1.0 / dx ** 2 + 1.0 / dy ** 2 + 1.0 / dz ** 2)


def kvec(m=1, n=1, Lx=2.0, Lz=3.0):
    return 2 * math.pi * m / Lx, 2 * math.pi * n / Lz


def Kt_of(geo, m=1, n=1):
    kx, kz = kvec(m, n, geo.Lx, geo.Lz)
    dx, dz = geo.hx[0], geo.hz[0]
    return math.hypot(tmm.ktilde(kx, dx), tmm.ktilde(kz, dz))


def kt_cont(m=1, n=1, Lx=2.0, Lz=3.0):
    return math.hypot(*kvec(m, n, Lx, Lz))


# ----------------------------------------------------------------------------- predictions store
def load_predictions():
    if not os.path.exists(PRED_FILE):
        return {}
    with open(PRED_FILE, encoding="utf-8") as fh:
        return {r["id"]: r for r in json.load(fh)["records"]}


def pred(pid, prefix="v2/"):
    p = load_predictions()
    pid = prefix + pid
    if pid not in p:
        raise KeyError(f"prediction {pid} missing: run tests/nu_predict.py before the gate")
    return p[pid]["value"]


# ----------------------------------------------------------------------------- result tables
class Table:
    """Rows: item | quantity | theory | measured | threshold | PASS/FAIL/INFO | grid info | note."""

    def __init__(self, gate):
        self.gate = gate
        self.rows = []

    def row(self, item, quantity, theory, measured, threshold, passed, grid="", note=""):
        self.rows.append(dict(item=item, quantity=quantity, theory=theory, measured=measured, threshold=threshold,
                              passed=None if passed is None else bool(passed), grid=grid, note=note))
        flag = "INFO" if passed is None else ("PASS" if passed else "FAIL")
        print(f"[{item}] {quantity}: theory={theory} measured={measured} thr={threshold} -> {flag} {note}",
              flush=True)
        return passed

    def ok(self):
        g = [r for r in self.rows if r["passed"] is not None]
        return bool(g) and all(r["passed"] for r in g)

    def save(self, extra=None, name=None):
        os.makedirs(RESULTS, exist_ok=True)
        path = os.path.join(RESULTS, (name or f"nu_{self.gate}") + ".json")
        out = dict(gate=self.gate, created_utc=utcnow(), table=self.rows, passed=self.ok(), **(extra or {}))
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, ensure_ascii=False, default=_js)
        g = [r for r in self.rows if r["passed"] is not None]
        print(f"\n{self.gate}: {sum(r['passed'] for r in g)}/{len(g)} gates PASS ->", "PASS" if self.ok() else "FAIL")
        return path


def _js(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, complex):
        return [o.real, o.imag]
    return str(o)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def grid_info(geo, dt=None):
    hy = geo.hy
    rr = hy[1:] / hy[:-1]
    r = float(max(rr.max(), (1 / rr).max())) if len(rr) else 1.0
    s = f"Δy_min=λ0/{1 / hy.min():.4g}, Δy_max=λ0/{1 / hy.max():.4g}, r_max={r:.4f}, Δx=λ0/{1 / geo.hx[0]:.4g}"
    if dt is not None:
        s += f", Δt={dt:.6g} T0"
    return s


def order_fit(h, err):
    """log-log slope of err vs h (least squares)."""
    return float(np.polyfit(np.log(np.asarray(h, float)), np.log(np.abs(np.asarray(err, float))), 1)[0])


# ----------------------------------------------------------------------------- running the C solver
BIN = os.path.join(ROOT, "fdtd3d_oblique")
BIN_REF = os.path.join(ROOT, "fdtd3d_oblique_ref")
COMPS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")
FRESH = os.environ.get("FDTD_FRESH", "0") == "1"


def build():
    import subprocess
    subprocess.run(["make", "-s"], cwd=ROOT, check=True)
    subprocess.run(["make", "-s", "ref"], cwd=ROOT, check=True)


def run(outdir, grid=None, binary=BIN, expect_fail=False, **kw):
    """Run the solver (key=value args; grid = grid name from grids/). Reuses a previous run only if the argument
    list, the grid hash and the binary are unchanged. Returns (meta or None, returncode, stderr)."""
    import subprocess
    os.makedirs(outdir, exist_ok=True)
    if grid is not None:
        g = load_grid(grid)
        kw = dict(grid=grid_path(grid), **kw)
        ghash = g["hash"]
    else:
        ghash = None
    argv = [f"{k}={','.join(map(str, v)) if isinstance(v, (list, tuple)) else v}" for k, v in kw.items()]
    argfile, metafile = os.path.join(outdir, "args.json"), os.path.join(outdir, "meta.json")
    key = dict(argv=argv, grid_hash=ghash, binary=os.path.basename(binary))
    if not FRESH and not expect_fail and os.path.exists(argfile) and os.path.exists(metafile):
        with open(argfile) as fh:
            same = json.load(fh) == key
        if same and os.path.getmtime(metafile) >= os.path.getmtime(binary):
            return load_meta(outdir), 0, ""
    res = subprocess.run([binary] + argv + [f"out={outdir}"], capture_output=True, text=True)
    if res.returncode != 0 and not expect_fail:
        raise RuntimeError(f"solver failed ({res.returncode}):\n{res.stderr}")
    if res.returncode == 0:
        with open(argfile, "w") as fh:
            json.dump(key, fh)
        return load_meta(outdir), 0, res.stderr
    return None, res.returncode, res.stderr


def load_meta(outdir):
    with open(os.path.join(outdir, "meta.json")) as fh:
        return json.load(fh)


def load_used(outdir):
    """grid_used.json -> dict of numpy arrays (the only coordinate source for nonuniform analysis)."""
    with open(os.path.join(outdir, "grid_used.json")) as fh:
        u = json.load(fh)
    return {k: (np.array(v) if isinstance(v, list) else v) for k, v in u.items()}


def load_log(outdir):
    return np.genfromtxt(os.path.join(outdir, "log.csv"), delimiter=",", names=True)


def load_fields(outdir, name, meta=None):
    meta = meta or load_meta(outdir)
    shape = (meta["Nx"], meta["Ny"] + 1, meta["Nz"])
    raw = np.fromfile(os.path.join(outdir, name), dtype=np.float64).reshape((6,) + shape)
    return {c: raw[q] for q, c in enumerate(COMPS)}


def load_auxref_dump(outdir, step, meta=None):
    meta = meta or load_meta(outdir)
    raw = np.fromfile(os.path.join(outdir, f"auxref_n{step}.bin"), dtype=np.float64)
    z = (raw[0::2] + 1j * raw[1::2]).reshape(6, meta["Ny"] + 1)
    return {c: z[q] for q, c in enumerate(COMPS)}


def load_auxref_dft(outdir, meta=None):
    meta = meta or load_meta(outdir)
    raw = np.fromfile(os.path.join(outdir, "dft_auxref.bin"), dtype=np.float64)
    z = (raw[0::2] + 1j * raw[1::2]).reshape(6, meta["Ny"] + 1)
    return {c: z[q] for q, c in enumerate(COMPS)}


def load_dft(outdir, name, meta=None):
    meta = meta or load_meta(outdir)
    sl = next(s for s in meta["slices"] if s["name"] == name)
    Nx, Ny, Nz = meta["Nx"], meta["Ny"], meta["Nz"]
    shape = {0: (Nx, Nz), 1: (Nx, Ny + 1), 2: (Ny + 1, Nz)}[sl["type"]]
    raw = np.fromfile(os.path.join(outdir, f"dft_{name}.bin"), dtype=np.float64)
    z = (raw[0::2] + 1j * raw[1::2]).reshape((6,) + shape)
    return {c: z[q] for q, c in enumerate(COMPS)}, sl


def xz_phase(used, comp, kx, kz):
    """exp(i(kx x_c + kz z_c)) on the component's own (x, z) sample positions, shape (Nx, Nz)."""
    x = used["x"][:-1] if comp in ("Ey", "Ez", "Hx") else used["x_dual"]
    z = used["z"][:-1] if comp in ("Ex", "Ey", "Hz") else used["z_dual"]
    return np.exp(1j * (kx * x[:, None] + kz * z[None, :]))


def y_of(used, comp):
    """y coordinates of a component's samples (primal for Ex, Ez, Hy; dual for Ey, Hx, Hz)."""
    return used["y"] if comp in ("Ex", "Ez", "Hy") else used["y_dual"]
