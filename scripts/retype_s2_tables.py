#!/usr/bin/env python3
"""Rewrite the tabulated solid-angle grids to the ixx_real literal policy.

Two transformations, both mechanical and idempotent:

1. Element type `T` -> `ixx_real`. The tables are `static constexpr`, which
   requires a literal type, so they cannot be templated on types such as
   `boost::numeric::interval`. Storing them as `ixx_real` and converting on read
   (see util/copy_grid.hpp) keeps them usable from constant expressions while
   leaving the quadratures type-generic. Under ENABLE_STRING_REALS `ixx_real`
   becomes `std::string_view`, so the tables then carry their exact decimal
   source text rather than a pre-rounded `double`.

2. Every table entry wrapped in `IXX_REAL(...)`, which is what makes (1) work in
   both modes.

Note that entries are wrapped in IXX_REAL *including* the handful whose value is
integral (0 and +/-1 at the axis points). Those are integral literals and would
otherwise be spelled IXX_INT, but a std::array is homogeneous, so they must
share the element type of their neighbors. Users worried about tracking ULPs
would likely want these values as integers to avoid any potential rounding 
errors caused by converting to floating-point values; however, since 0, +/-1 are
exactly representable as floating-point numbers, no precision is lost by storing
them as `ixx_real` instead of `ixx_int`.

The surrounding `template <typename T>` is left in place, though T is no longer
used by the table itself: keeping it means none of the ~30 dispatch branches in
each of the four family headers has to change.

Re-run after syncing tables from upstream.
"""
import argparse
import pathlib
import re
import sys

FAMILIES = ("lebedev_laikov", "delley", "ahrens_beylkin", "womersley")

# Declarations: T (upstream).
POINTS = re.compile(
    r"static constexpr std::array<cartesian_pt_t<T>,(\s*)(\d+)> points")
WEIGHTS = re.compile(
    r"static constexpr std::array<T,(\s*)(\d+)> weights")

# Womersley grids are equal-weight and so compute rather than tabulate their
# weights: `create_array<N, T>(4.0 * M_PI / N.0)`. That form is templated on T
# (so it breaks for non-literal types), pre-divides in `double`, and 
# materializes N identical values. Replace it with the exact integer N, and let 
# copy_grid form 4*pi/N via divide_integer -- one rounding instead of two, 
# correct in string mode.
UNIFORM = re.compile(
    r"[ \t]*static constexpr auto weights = *\n?"
    r"[ \t]*detail::create_array<(\d+), ?T>\([^;]*\);")

# An initializer body: everything between "= {" and the closing "};".
BODY = re.compile(r"(=\s*\{)(.*?)(\}\s*;)", re.S)

# A numeric literal, with or without an exponent.
NUM = re.compile(r"(?<![\w.(])([-+]?\d+\.\d+(?:[eE][-+]?\d+)?)")


def wrap_bodies(src: str) -> str:
    def repl(m: "re.Match[str]") -> str:
        head, body, tail = m.group(1), m.group(2), m.group(3)
        if "IXX_REAL(" in body:          # already wrapped; keep idempotent
            return m.group(0)
        return head + NUM.sub(r"IXX_REAL(\1)", body) + tail
    return BODY.sub(repl, src)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="include/integratorxx/quadratures/s2")
    ap.add_argument("--check", action="store_true",
                    help="report what would change; do not write")
    args = ap.parse_args()

    root = pathlib.Path(args.root)
    if not root.is_dir():
        print(f"error: no such directory: {root}", file=sys.stderr)
        return 2

    files = sorted(p for fam in FAMILIES for p in (root / fam).glob("*.hpp")
                   if not p.name.endswith("_grids.hpp"))
    if not files:
        print(f"error: no table files found under {root}", file=sys.stderr)
        return 2

    changed = []
    for path in files:
        src = path.read_text()
        out = POINTS.sub(
            r"static constexpr std::array<cartesian_pt_t<ixx_real>,\1\2> points", src)
        out = WEIGHTS.sub(
            r"static constexpr std::array<ixx_real,\1\2> weights", out)
        out = wrap_bodies(out)
        out = UNIFORM.sub(
            r"  /// Equal-weight grid: every weight is 4*pi/\1 (sphere area / npts).\n"
            r"  /// copy_grid forms it exactly; see util/copy_grid.hpp.\n"
            r"  static constexpr ixx_int uniform_weight_npts = \1;", out)
        if out != src:
            changed.append(path)
            if not args.check:
                path.write_text(out)

    verb = "would rewrite" if args.check else "rewrote"
    print(f"{verb} {len(changed)} of {len(files)} table files")
    return 1 if (args.check and changed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
