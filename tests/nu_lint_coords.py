"""Lint: nonuniform analysis code must take coordinates from grid files / grid_used.json, never idx*Delta.

    python3 tests/nu_lint_coords.py [files...]      # default: every nonuniform-path Python file

Scans new modules completely and analyze.py / fdtd_io.py / compare.py after their '# ==== nonuniform' marker.
A line may opt out with the comment '# lint-ok: <reason>' (used only by the uniform generators). Exit 1 on hits.
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATTERNS = [
    (re.compile(r"\*\s*(D|Delta|delta_ref)\b"), "multiplication by a uniform spacing"),
    (re.compile(r"\b(D|Delta)\s*\*"), "multiplication by a uniform spacing"),
    (re.compile(r"/\s*(nl|n_lambda)\b"), "spacing derived from cells per wavelength"),
    (re.compile(r"\[\s*[\"']Delta[\"']\s*\]"), "meta['Delta'] used as a coordinate scale"),
]
MARKER = "# ==== nonuniform"


def default_files():
    fs = glob.glob(os.path.join(ROOT, "tests", "nu*.py")) + glob.glob(os.path.join(ROOT, "tests", "test_*.py"))
    for f in ("grid_gen.py", "tmm.py", "rcwa.py", "meep_ref_uniform.py", "meep_ref_transform.py"):
        p = os.path.join(ROOT, f)
        if os.path.exists(p):
            fs.append(p)
    partial = [os.path.join(ROOT, f) for f in ("analyze.py", "fdtd_io.py", "compare.py")]
    return [(f, False) for f in fs if not f.endswith("nu_lint_coords.py")] + [(f, True) for f in partial]


def scan(path, after_marker):
    hits = []
    active = not after_marker
    in_doc = False
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if MARKER in line:
                active = True
                continue
            quotes = line.count('"""') + line.count("'''")
            if in_doc or quotes:                               # docstrings are prose, not code
                if quotes % 2 == 1:
                    in_doc = not in_doc
                continue
            if not active or "lint-ok" in line:
                continue
            code = line.split("#", 1)[0]
            code = re.sub(r"(\"[^\"]*\"|'[^']*')", "''", code)      # ignore string literals
            for rx, why in PATTERNS:
                if rx.search(code):
                    hits.append((path, n, line.rstrip(), why))
    return hits


def main(argv):
    files = [(f, False) for f in argv] if argv else default_files()
    hits = []
    for f, part in files:
        if os.path.exists(f):
            hits += scan(f, part)
    for p, n, line, why in hits:
        print(f"{os.path.relpath(p, ROOT)}:{n}: {why}: {line.strip()}")
    print(f"nu_lint_coords: {len(files)} files, {len(hits)} hits ->", "PASS" if not hits else "FAIL")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
