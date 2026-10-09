"""Read the tabulated angular grids out of the generated C++ headers."""
import re
from pathlib import Path

from mpmath import mp, mpf

# Data files mix plain decimals with exponent form, sometimes within one array
# (ahrens_beylkin_552.hpp has both). Requiring an exponent silently drops most
# of the values, so make it optional.
_NUM = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?")


def _numbers(text):
    return [mpf(t.replace("D", "e").replace("E", "e")) for t in _NUM.findall(text)]


def _section(text, name):
    i = text.index(name + " = {")
    j = text.index("};", i)
    return text[i:j]


_EQUAL = re.compile(r"create_array<\s*(\d+)\s*,\s*T\s*>\s*\(\s*4\.0\s*\*\s*M_PI\s*/\s*(\d+)")


def read_grid(path):
    """Return (points, weights) from a *_<npts>.hpp data header.

    Two weight formats occur. Lebedev-Laikov and Delley list the weights
    explicitly. The Womersley grids are equal-weight spherical designs and
    synthesise theirs with ``create_array<N, T>(4.0 * M_PI / N)``, so there is
    no array to read.
    """
    text = Path(path).read_text()
    flat = _numbers(_section(text, "points"))
    points = [tuple(flat[3 * i:3 * i + 3]) for i in range(len(flat) // 3)]

    m = _EQUAL.search(text)
    if m:
        n, denom = int(m.group(1)), int(m.group(2))
        if n != denom or n != len(points):
            raise ValueError(f"{path}: equal-weight size mismatch "
                             f"({n}, {denom}, {len(points)} points)")
        w = 4 * mp.pi / n
        return points, [w] * n

    weights = _numbers(_section(text, "weights"))
    if len(points) != len(weights):
        raise ValueError(f"{path}: {len(points)} points but {len(weights)} weights")
    return points, weights


def grid_path(root, family, npts):
    return Path(root) / "include/integratorxx/quadratures/s2" / family / f"{family}_{npts}.hpp"
