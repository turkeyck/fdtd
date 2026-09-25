"""Gate 0-2: adjointness <E, C_H H>_{W_E} = <C_E E, H>_{W_H} with random fields (SPEC_nonuniform §20.1),
plus a check that the C update kernel applies the same operator (two steps from init=r, dumped)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import analyze_nu as an  # noqa: E402

TH = nc.thresholds()["gate0"]["0-2"]


def used_of(grid):
    """grid_used arrays straight from the solver (short run) so the test sees exactly the solver's geometry."""
    out = os.path.join(nc.RUNS, "gate0", f"geom_{grid}")
    nc.run(out, grid=grid, inc="0", init="r", nsteps=2, dump_at=[1, 2], energy_every=0)
    return out, nc.load_used(out)


def main():
    nc.build()
    T = nc.Table("gate0_adjoint")
    for g in TH["grids"]:
        pec = g + "_pec_thin"
        out, used = used_of(pec)
        geo = nc.Geo(nc.load_grid(pec))
        W = an.weights(used)
        for seed in (1, 2, 3):
            F = an.random_state(used, seed)
            CH = an.curl_H(F, used)
            CE = an.curl_E(F, used)
            lhs = an.inner(F, CH, W, ("Ex", "Ey", "Ez"))
            rhs = an.inner(CE, F, W, ("Hx", "Hy", "Hz"))
            nE = np.sqrt(an.inner(F, F, W, ("Ex", "Ey", "Ez")))
            nCH = np.sqrt(an.inner(CH, CH, W, ("Ex", "Ey", "Ez")))
            rel = abs(lhs - rhs) / (nE * nCH)
            T.row("0-2", f"|<E,C_H H>_WE - <C_E E,H>_WH| / (|E| |C_H H|), seed {seed}", "0", f"{rel:.2e}",
                  f"< {TH['rel']:g}", rel < TH["rel"], grid=nc.grid_info(geo))
            one = {c: np.ones_like(W[c]) for c in W}
            l1 = an.inner(F, CH, one, ("Ex", "Ey", "Ez"))
            r1 = an.inner(CE, F, one, ("Hx", "Hy", "Hz"))
            T.row("0-2", f"control: unit weights, seed {seed}", "≠ 0", f"{abs(l1 - r1) / abs(l1):.2e}", "INFO",
                  None, grid=pec, note="shows the weights are what makes C_E, C_H adjoint")
        # C kernel == Python operator (dumps after steps 1 and 2 of a source-free run)
        meta = nc.load_meta(out)
        dt = meta["dt"]
        F1 = nc.load_fields(out, "fields_n1.bin", meta)
        F2 = nc.load_fields(out, "fields_n2.bin", meta)
        eps = an.eps_fields(used)
        CE1 = an.curl_E(F1, used)
        eH = max(np.abs(F2[c] - (F1[c] - dt * CE1[c])).max() / np.abs(F1[c]).max() for c in ("Hx", "Hy", "Hz"))
        CH2 = an.curl_H(F2, used)
        eE = max(np.abs(F2[c] - (F1[c] + dt * CH2[c] / eps[c])).max() / np.abs(F1[c]).max() for c in ("Ex", "Ey", "Ez"))
        T.row("0-2", "C kernel vs Python operator (one H and one E update)", "0", f"H {eH:.1e}, E {eE:.1e}",
              "< 1e-13", max(eH, eE) < 1e-13, grid=pec, note="ties the tested operator to the solver")
    T.save()
    return 0 if T.ok() else 1


if __name__ == "__main__":
    sys.exit(main())
