"""Run definitions shared by the nonuniform gates (one place, so every gate reuses the same cached runs).

All runs use dt = T0 / ceil(T0 / dt_Courant) (implementation note D12) and the erf ramp (t0 = 20 T0, tau = 4 T0).
"""
import math
import os

import nu_common as nc

PERIODS, DFT_PERIODS = 90, 30


def steps_per_period(geo):
    return int(round(1.0 / geo.dt_run()))


def level1(grid, pol, tag="", full=False, stability=False, **extra):
    """Level-1 vacuum run on `grid`: aux line + aux_ref, projected DFT of every y row, DFT window 60..90 T0.
    full=True adds the full-volume planarity DFT, divergence log, a final field dump and an x-y snapshot."""
    geo = nc.Geo(nc.load_grid(grid))
    dt = geo.dt_run()
    Np = steps_per_period(geo)
    out = os.path.join(nc.RUNS, "gate1", f"{grid}_{pol}{tag}")
    kw = dict(pol=pol, inc="a", auxref=1, proj=1, dt=repr(dt), energy_every=0, ref_every=10 * Np)
    if stability:
        n = max(20000, 160 * Np)
        kw.update(nsteps=n, off_t=40.0, energy_every=max(Np // 4, 1), dft0=0, dft1=0, proj=0, ref_every=0)
    else:
        kw.update(nsteps=PERIODS * Np, dft0=(PERIODS - DFT_PERIODS) * Np, dft1=PERIODS * Np)
    if full:
        kw.update(planar=1, div_every=Np, dump_at=[PERIODS * Np], snap=PERIODS * Np, snapcomp="Ez",
                  zk=geo.Nz // 3, xi=geo.Nx // 2)
    kw.update(extra)
    meta, _, _ = nc.run(out, grid=grid, **kw)
    return out, meta, geo


def load_proj(outdir, meta=None):
    meta = meta or nc.load_meta(outdir)
    import numpy as np
    raw = np.fromfile(os.path.join(outdir, "dft_proj.bin"), dtype=np.float64)
    z = (raw[0::2] + 1j * raw[1::2]).reshape(6, meta["Ny"] + 1)
    return {c: z[q] for q, c in enumerate(nc.COMPS)}


def load_planar(outdir, meta=None):
    meta = meta or nc.load_meta(outdir)
    import numpy as np
    v = np.fromfile(os.path.join(outdir, "dft_planar.bin"), dtype=np.float64).reshape(6, meta["Ny"] + 1)
    return {c: v[q] for q, c in enumerate(nc.COMPS)}


def front_window_steps(meta, used):
    """Steps before the incident front can reach the far PML (speed c, 10% margin), from the aux source."""
    y_far = used["y"][meta["Ny"] - meta["npml_hi"]]
    y_src = used["y"][0] + meta["ja"] * used["hy"][0] if meta["ja"] < 0 else used["y"][meta["ja"]]
    return int(math.floor(0.9 * (y_far - y_src) / meta["dt"]))


def level2(grid, pol, dt=None, tag="", **extra):
    """Level-2 (layered media) run: same source, window and dt rule as level 1; aux_ref checked every 30 T0."""
    geo = nc.Geo(nc.load_grid(grid))
    if dt is None:
        dt = geo.dt_run()
        if "_bis" in grid:                           # D10 / D14: dt scales with h
            k = int(grid.split("_bis")[1].lstrip("xz"))
            dt = nc.Geo(nc.load_grid(grid.split("_bis")[0] + "_bis1")).dt_run() / k
    Np = int(round(1.0 / dt))
    out = os.path.join(nc.RUNS, "gate2", f"{grid}_{pol}{tag}")
    kw = dict(pol=pol, inc="a", auxref=1, proj=1, dt=repr(dt), energy_every=0, ref_every=30 * Np,
              nsteps=PERIODS * Np, dft0=(PERIODS - DFT_PERIODS) * Np, dft1=PERIODS * Np)
    kw.update(extra)
    meta, _, _ = nc.run(out, grid=grid, **kw)
    return out, meta, geo


def stageB(grid, pol, inc="p", tag="", periods=70, dft_periods=20, m=1, n=0, **extra):
    """Stage-B run (x and/or z nonuniform): modal injection inc=m (D22), analytic injection inc=p or current sheet
    inc=j, weighted Floquet
    projection of every y row, DFT window of the last `dft_periods` periods."""
    geo = nc.Geo(nc.load_grid(grid))
    dt = geo.dt_run()
    Np = int(round(1.0 / dt))
    out = os.path.join(nc.RUNS, "stageB", f"{grid}_{pol}_{inc}{tag}")
    kw = dict(pol=pol, inc=inc, m=m, n=n, proj=1, dt=repr(dt), energy_every=0, nsteps=periods * Np,
              dft0=(periods - dft_periods) * Np, dft1=periods * Np)
    kw.update(extra)
    meta, _, _ = nc.run(out, grid=grid, **kw)
    return out, meta, geo
