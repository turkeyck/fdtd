"""Re-run only the gate-3 (B) Meep jobs (after the geometry-block fix)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_gate3_meep as g3  # noqa: E402

for case in ("vac", "film"):
    grid = f"L3B_mapped_{case}"
    dt = nc.Geo(nc.load_grid(grid)).dt_run()
    for pol in ("s", "p"):
        for variant in ("cont", "disc"):
            out = os.path.join(g3.MROOT, f"B_{case}_{pol}_{variant}")
            g3.meep("meep_ref_transform.py", out, ["--grid", nc.grid_path(grid), "--pol", pol, "--dt", repr(dt),
                                                   "--variant", variant], "meep_B.json")
            print("B", case, pol, variant, flush=True)
