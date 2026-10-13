"""Regression tests: the tool must reproduce the tabulated rules it starts from.

Run with:  PYTHONPATH=tools python3 -m angular_grids.test_refine
"""
import sys
from pathlib import Path

from mpmath import mp, mpf

from . import designs
from .moments import conditions
from .refine import decompose, materialise, refine, verify
from .tables import grid_path, read_grid

CASES = [("lebedev_laikov", 6, 3), ("lebedev_laikov", 50, 11), ("lebedev_laikov", 110, 17),
         ("lebedev_laikov", 194, 23), ("lebedev_laikov", 302, 29),
         ("delley", 50, 11), ("delley", 194, 23), ("delley", 302, 29)]

# Equal-weight spherical designs. Kept small: the inner CGLS is the cost, and
# it grows steeply with N.
DESIGN_CASES = [(14, 4), (18, 5), (26, 6)]


def main(root="."):
    mp.dps = 60
    failed = 0
    for family, npts, order in CASES:
        pts, w = read_grid(grid_path(root, family, npts))
        orbits = decompose(pts, [x * 4 * mp.pi / sum(w) for x in w])
        n_param = sum(1 + o.type.n_param for o in orbits)
        refine(orbits, order, 45)
        P, W = materialise(orbits)

        ok = len(P) == npts
        # verified against the FULL redundant condition set, not the reduced one
        err = verify(P, W, order)
        ok = ok and err < mpf("1e-45")
        # the refined rule must stay close to the tabulated one, not wander to
        # a different rule of the same order
        drift = max(min(max(abs(a - b) for a, b in zip(p, q)) for q in P) for p in pts)
        ok = ok and drift < mpf("1e-13")

        print(f"  {'ok  ' if ok else 'FAIL'} {family}_{npts:<5d} order {order:3d}  "
              f"{n_param:3d} params  err {mp.nstr(err, 3):>10s}  drift {mp.nstr(drift, 3):>10s}")
        failed += not ok
    total = len(CASES)
    for npts, order in DESIGN_CASES:
        pts, _ = read_grid(grid_path(root, "womersley", npts))
        d = designs.Design(pts, order)
        hist = designs.refine(d, 22)
        err = designs.verify(d)

        ok = len(d.points) == npts and err < mpf("1e-20")
        # a design keeps its equal weights and stays on the sphere
        ok = ok and all(abs(sum(v * v for v in p) - 1) < mpf("1e-25") for p in d.points)
        print(f"  {'ok  ' if ok else 'FAIL'} womersley_{npts:<5d} order {order:3d}  "
              f"{2 * npts:3d} unknowns {len(d.conds):3d} conds  "
              f"{len(hist)} iters  err {mp.nstr(err, 3):>10s}")
        failed += not ok
        total += 1

    print(f"\n  {total - failed}/{total} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
