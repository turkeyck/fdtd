"""Run fdtd3d_oblique and read its outputs (meta.json, dft_*.bin, log.csv, snapshots, dumps)."""
import json
import os
import subprocess

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "fdtd3d_oblique")
COMPS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")


def build():
    subprocess.run(["make", "-s"], cwd=ROOT, check=True)


FRESH = os.environ.get("FDTD_FRESH", "0") == "1"


def run(outdir, quiet=True, cache=True, **kw):
    """Run the solver with key=value args; returns the parsed meta.json.

    A previous run in outdir is reused only if it was made with identical arguments by a binary no newer
    than its meta.json; set FDTD_FRESH=1 to force re-running everything (used by the regression driver).
    """
    os.makedirs(outdir, exist_ok=True)
    argv = [f"{k}={','.join(map(str, v)) if isinstance(v, (list, tuple)) else v}" for k, v in kw.items()]
    argfile, metafile = os.path.join(outdir, "args.json"), os.path.join(outdir, "meta.json")
    if cache and not FRESH and os.path.exists(argfile) and os.path.exists(metafile):
        with open(argfile) as fh:
            same = json.load(fh) == argv
        if same and os.path.getmtime(metafile) >= os.path.getmtime(BIN):
            return load_meta(outdir)
    res = subprocess.run([BIN] + argv + [f"out={outdir}"], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"solver failed ({res.returncode}):\n{res.stderr}")
    with open(argfile, "w") as fh:
        json.dump(argv, fh)
    if not quiet:
        print(res.stderr.strip().splitlines()[-1])
    return load_meta(outdir)


def load_meta(outdir):
    with open(os.path.join(outdir, "meta.json")) as fh:
        return json.load(fh)


def slice_shape(meta, sl):
    Nx, Ny, Nz = meta["Nx"], meta["Ny"], meta["Nz"]
    return {0: (Nx, Nz), 1: (Nx, Ny + 1), 2: (Ny + 1, Nz)}[sl["type"]]


def load_dft(outdir, name, meta=None):
    meta = meta or load_meta(outdir)
    sl = next(s for s in meta["slices"] if s["name"] == name)
    shape = slice_shape(meta, sl)
    raw = np.fromfile(os.path.join(outdir, f"dft_{name}.bin"), dtype=np.float64)
    z = (raw[0::2] + 1j * raw[1::2]).reshape((6,) + shape)
    return {c: z[q] for q, c in enumerate(COMPS)}, sl


def load_aux_dft(outdir, meta=None):
    meta = meta or load_meta(outdir)
    raw = np.fromfile(os.path.join(outdir, "dft_aux.bin"), dtype=np.float64)
    z = raw[0::2] + 1j * raw[1::2]
    lo, hi = meta["aux"]["dft_jlo"], meta["aux"]["dft_jhi"]
    z = z.reshape(6, hi - lo + 1)
    return {c: z[q] for q, c in enumerate(COMPS)}, np.arange(lo, hi + 1)


def load_fields(outdir, meta=None):
    meta = meta or load_meta(outdir)
    shape = (meta["Nx"], meta["Ny"] + 1, meta["Nz"])
    raw = np.fromfile(os.path.join(outdir, "fields_final.bin"), dtype=np.float64).reshape((6,) + shape)
    return {c: raw[q] for q, c in enumerate(COMPS)}


def load_log(outdir):
    return np.genfromtxt(os.path.join(outdir, "log.csv"), delimiter=",", names=True)


def load_snapshots(outdir, meta=None):
    meta = meta or load_meta(outdir)
    raw = np.fromfile(os.path.join(outdir, "snapshots.bin"), dtype=np.float64)
    return raw.reshape(meta["nsnap"], meta["Nx"], meta["Ny"] + 1)


def slice_coords(meta, sl, comp):
    """Physical coordinates (x, y, z) of comp's samples in slice sl, broadcast to the slice shape."""
    D = meta["Delta"]
    ox, oy, oz, _ = meta["offsets"][comp]
    Nx, Ny, Nz = meta["Nx"], meta["Ny"], meta["Nz"]
    if sl["type"] == 0:
        i = np.arange(Nx)[:, None]
        k = np.arange(Nz)[None, :]
        x, z = (i + ox) * D, (k + oz) * D
        y = np.full(np.broadcast(x, z).shape, (sl["idx"] + oy) * D)
        return np.broadcast_to(x, y.shape), y, np.broadcast_to(z, y.shape)
    if sl["type"] == 1:
        i = np.arange(Nx)[:, None]
        j = np.arange(Ny + 1)[None, :]
        x, y = (i + ox) * D, (j + oy) * D
        z = np.full(np.broadcast(x, y).shape, (sl["idx"] + oz) * D)
        return np.broadcast_to(x, z.shape), np.broadcast_to(y, z.shape), z
    j = np.arange(Ny + 1)[:, None]
    k = np.arange(Nz)[None, :]
    y, z = (j + oy) * D, (k + oz) * D
    x = np.full(np.broadcast(y, z).shape, (sl["idx"] + ox) * D)
    return x, np.broadcast_to(y, x.shape), np.broadcast_to(z, x.shape)


def valid_j_mask(meta, sl, comp):
    """Mask of meaningful samples: half-integer-y components have no sample at j = Ny."""
    shape = slice_shape(meta, sl)
    m = np.ones(shape, dtype=bool)
    if meta["offsets"][comp][1] == 0.5:
        if sl["type"] == 1:
            m[:, -1] = False
        elif sl["type"] == 2:
            m[-1, :] = False
    return m
