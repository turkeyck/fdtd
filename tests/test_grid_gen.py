"""Unit tests for grid_gen (SPEC_nonuniform §24: T-G1, T-G2, T-G3, T-L1). numpy only; exit 1 on failure."""
import copy
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import grid_gen as gg  # noqa: E402

FAIL = []


def check(name, cond, msg):
    print(f"[{name}] {'PASS' if cond else 'FAIL'}: {msg}")
    if not cond:
        FAIL.append(name)


def t_g1():
    for struct in ("F1", "F5"):
        g = gg.level2(struct, 40)
        geo = nc.Geo(g)
        itf = g["interfaces_j"]
        ok = True
        for q in range(0, len(itf) - 1):
            L = geo.y[itf[q + 1]] - geo.y[itf[q]]
            want = gg.STACK5[q][1] if struct == "F5" else gg.D_FILM
            ok &= abs(L - want) < 1e-15
        film_cells = itf[1] - itf[0]
        check("T-G1", ok and (struct != "F1" or film_cells == 7),
              f"{struct}: interfaces on primal nodes, layer thicknesses exact (<1e-15), film cells = {film_cells}")


def t_g2():
    worst = 1.0
    for g in gg.all_grids():
        h = np.array(g["y"]["h"])
        rr = h[1:] / h[:-1]
        r = g["generator"]["r_max"]
        worst = max(worst, float(np.max(np.maximum(rr, 1 / rr)) / r))
    check("T-G2", worst <= 1 + 1e-12, f"all grids: max(adjacent ratio / r_max) = {worst:.12f}")


def t_g3():
    g = gg.level1(20, 1.1)
    bad = copy.deepcopy(g)
    c0, n = g["marks"]["down"]
    bad["zones"]["j0"] = c0 + n // 2               # TF/SF plane inside the taper
    bad["hash"] = gg.canonical_hash(bad)
    v = gg.check(bad)
    check("T-G3", any("j0" in s for s in v), f"j0 inside the taper rejected: {v}")
    bad2 = copy.deepcopy(g)
    bad2["y"]["h"][100] *= 1.0000001
    check("T-G3b", any("hash" in s for s in gg.check(bad2)), "hand-edited spacing rejected by hash check")
    bad3 = copy.deepcopy(g)
    bad3["generator"]["ppw"] = 40
    bad3["hash"] = gg.canonical_hash(bad3)
    check("T-G3c", any("lambda0/(n ppw)" in s for s in gg.check(bad3)), "spacing above lambda0/(n ppw) rejected")


def t_l1():
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write("import numpy as np\nD = 0.05\nx = np.arange(10) * D\n")
        path = fh.name
    rc = subprocess.run([sys.executable, os.path.join(nc.ROOT, "tests", "nu_lint_coords.py"), path],
                        capture_output=True, text=True).returncode
    os.unlink(path)
    check("T-L1", rc == 1, "lint flags 'np.arange(n) * D'")
    rc = subprocess.run([sys.executable, os.path.join(nc.ROOT, "tests", "nu_lint_coords.py")],
                        capture_output=True, text=True)
    check("T-L1b", rc.returncode == 0, "lint passes on the nonuniform code base:\n" + rc.stdout.strip())


if __name__ == "__main__":
    t_g1()
    t_g2()
    t_g3()
    t_l1()
    print("test_grid_gen:", "PASS" if not FAIL else f"FAIL {FAIL}")
    sys.exit(1 if FAIL else 0)
