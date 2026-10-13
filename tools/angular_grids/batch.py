"""Regenerate every tabulated size of a family.

    PYTHONPATH=tools python3 -m angular_grids.batch --family delley --digits 40 \
        --out gen/ --jobs 8

Each size is independent, so they run in separate processes. A size that
fails to reach the requested precision is reported and skipped rather than
written, so a partial run never emits a table that was not verified.
"""
import argparse
import multiprocessing as mp_proc
import re
import sys
import time
from pathlib import Path


def tabulated_sizes(family, root):
    """Sizes the family's dispatch header advertises, via its order table."""
    text = (Path(root) / "include/integratorxx/quadratures/s2" / f"{family}.hpp").read_text()
    body = text[text.index("algebraic_order_by_npts"):]
    body = body[:body.index("next_algebraic_order")] if "next_algebraic_order" in body else body
    return sorted({int(m) for m in re.findall(r"case\s+(\d+)\s*:", body)})


def _one(args):
    family, npts, digits, root, out = args
    from .generate import run
    t = time.time()
    try:
        _, err = run(family, npts, digits, root, out, verbose=False)
        from mpmath import mp
        return (npts, True, mp.nstr(err, 4), time.time() - t)
    except SystemExit as e:
        return (npts, False, str(e), time.time() - t)
    except Exception as e:                      # noqa: BLE001 - report, do not abort the batch
        return (npts, False, f"{type(e).__name__}: {e}", time.time() - t)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--family", required=True,
                    choices=["lebedev_laikov", "delley", "ahrens_beylkin", "womersley"])
    ap.add_argument("--digits", type=int, default=40)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=max(1, (mp_proc.cpu_count() or 2) - 1))
    ap.add_argument("--max-npts", type=int, default=None,
                    help="skip sizes above this, for a quick pass")
    a = ap.parse_args(argv)

    sizes = tabulated_sizes(a.family, a.root)
    if a.max_npts:
        sizes = [n for n in sizes if n <= a.max_npts]
    print(f"{a.family}: {len(sizes)} sizes, {a.digits} digits, {a.jobs} jobs", flush=True)

    work = [(a.family, n, a.digits, a.root, a.out) for n in sizes]
    ok = bad = 0
    with mp_proc.Pool(a.jobs) as pool:
        for npts, good, msg, dt in pool.imap_unordered(_one, work):
            if good:
                ok += 1
                print(f"  ok    {a.family}_{npts:<6d} err {msg:<12s} [{dt:5.0f}s]", flush=True)
            else:
                bad += 1
                print(f"  FAIL  {a.family}_{npts:<6d} {msg} [{dt:5.0f}s]", flush=True)
    print(f"\n{ok} regenerated, {bad} failed", flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
