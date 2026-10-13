"""Point a family's dispatch header at regenerated string tables.

    PYTHONPATH=tools python3 -m angular_grids.switch_to_load --family delley --root .

Rewrites

    detail::copy_grid<delley_50<RealType>>(points, weights);
        ->  delley_50<RealType>::load(points, weights);

and drops the `weights[i] *= 4.0*M_PI` normalisation loop, because a
regenerated table already carries the 4*pi convention while the shipped
literal tables are normalised to one.

This only makes sense once the data headers have actually been regenerated;
run angular_grids.batch first. The change is textual and reversible with git.
"""
import argparse
import re
import sys
from pathlib import Path

FAMILIES = ("lebedev_laikov", "delley", "ahrens_beylkin", "womersley")

# families whose shipped tables are normalised to 1 and scaled at load
SCALED = ("lebedev_laikov", "delley")

# The loop index is declared `auto` in delley.hpp and `size_t` in
# lebedev_laikov.hpp, so do not pin the type.
_SCALE_LOOP = re.compile(
    r"\n[ \t]*//[^\n]*4 ?pi[^\n]*\n"
    r"[ \t]*for\s*\(\s*\w+\s+i\s*=\s*0\s*;\s*i\s*<\s*npts\s*;\s*i\+\+\s*\)\s*\n"
    r"[ \t]*weights\[i\]\s*\*=\s*4\.0\s*\*\s*M_PI\s*;\n",
    re.IGNORECASE)


def switch(family, root, dry_run=False):
    path = Path(root) / "include/integratorxx/quadratures/s2" / f"{family}.hpp"
    text = path.read_text()

    calls = re.findall(r"detail::copy_grid<\s*(\w+)<RealType>\s*>\s*\(\s*points\s*,\s*weights\s*\)",
                       text)
    new = re.sub(r"detail::copy_grid<\s*(\w+)<RealType>\s*>\s*\(\s*points\s*,\s*weights\s*\)",
                 r"\1<RealType>::load(points, weights)", text)

    dropped = 0
    if family in SCALED:
        new, dropped = _SCALE_LOOP.subn("\n", new)

    if not dry_run:
        path.write_text(new)
    return len(calls), dropped, path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--family", required=True, choices=FAMILIES)
    ap.add_argument("--root", default=".")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    n, dropped, path = switch(a.family, a.root, a.dry_run)
    what = "would rewrite" if a.dry_run else "rewrote"
    print(f"{what} {n} dispatch calls in {path}"
          + (f", dropped {dropped} 4*pi scaling loop(s)" if dropped else ""))
    return 0 if n else 1


if __name__ == "__main__":
    sys.exit(main())
