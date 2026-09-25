"""Pre-run every Meep job of gate 3 (single-threaded Meep; results cached under runs/nu/meep, reused by
tests/nu_gate3_meep.py). Useful to overlap the Meep work with the C runs of other gates."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nu_common as nc  # noqa: E402
import nu_gate3_meep as g3  # noqa: E402


def main():
    for pol in ("s", "p"):
        for st in ("F1", "F5"):
            for r in nc.thresholds()["gate3"]["3A-2"]["resolutions_D18"]:
                g3.meep_A(st, pol, r)
                print("A", st, pol, r, flush=True)
    for case in ("vac", "film"):
        grid = f"L3B_mapped_{case}"
        dt = nc.Geo(nc.load_grid(grid)).dt_run()
        for pol in ("s", "p"):
            for variant in ("cont", "disc"):
                out = os.path.join(g3.MROOT, f"B_{case}_{pol}_{variant}")
                g3.meep("meep_ref_transform.py", out, ["--grid", nc.grid_path(grid), "--pol", pol, "--dt", repr(dt),
                                                       "--variant", variant], "meep_B.json")
                print("B", case, pol, variant, flush=True)
    g3.meep_A("F1", "s", 40, "thin")
    g3.meep_A("F1", "s", 40, "full")
    print("full cell done", flush=True)


if __name__ == "__main__":
    main()
