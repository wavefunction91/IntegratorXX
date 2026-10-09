"""Octahedral (O_h) orbit algebra for Lebedev-Laikov and Delley grids.

Both families are invariant under the full octahedral group with inversion, so
their points fall into a handful of orbit types. Each orbit carries one weight
and zero, one or two free angular parameters -- which is why a 5810-point
Lebedev rule has 385 free parameters rather than 17430.

Orbit types follow the standard Lebedev naming:

    a1   6 pts   (1,0,0)      0 parameters
    a2  12 pts   (1,1,0)/r2   0 parameters
    a3   8 pts   (1,1,1)/r3   0 parameters
    bk  24 pts   (l,l,m)      1 parameter,  m = sqrt(1-2l^2)
    ck  24 pts   (p,q,0)      1 parameter,  q = sqrt(1-p^2)
    dk  48 pts   (r,s,t)      2 parameters, t = sqrt(1-r^2-s^2)
"""
from itertools import permutations

from mpmath import mp, mpf, sqrt


class OrbitType:
    """One octahedral orbit type: how to build its points from its parameters."""

    def __init__(self, name, n_param, base, dbase, size):
        self.name = name
        self.n_param = n_param
        self._base = base    # params -> (v0, v1, v2), the |coordinates|
        self._dbase = dbase  # params -> [[dv_j/dparam_k]], shape 3 x n_param
        self.size = size

    def base(self, params):
        return self._base(params)

    def dbase(self, params):
        return self._dbase(params)


def _no_deriv(_):
    return [[], [], []]


A1 = OrbitType("a1", 0, lambda _: (mpf(1), mpf(0), mpf(0)), _no_deriv, 6)
A2 = OrbitType("a2", 0, lambda _: (1 / sqrt(2), 1 / sqrt(2), mpf(0)), _no_deriv, 12)
A3 = OrbitType("a3", 0, lambda _: (1 / sqrt(3),) * 3, _no_deriv, 8)


def _b_base(p):
    l = p[0]
    return (l, l, sqrt(1 - 2 * l * l))


def _b_dbase(p):
    l = p[0]
    m = sqrt(1 - 2 * l * l)
    return [[mpf(1)], [mpf(1)], [-2 * l / m]]


def _c_base(p):
    x = p[0]
    return (x, sqrt(1 - x * x), mpf(0))


def _c_dbase(p):
    x = p[0]
    q = sqrt(1 - x * x)
    return [[mpf(1)], [-x / q], [mpf(0)]]


def _d_base(p):
    r, s = p
    return (r, s, sqrt(1 - r * r - s * s))


def _d_dbase(p):
    r, s = p
    t = sqrt(1 - r * r - s * s)
    return [[mpf(1), mpf(0)], [mpf(0), mpf(1)], [-r / t, -s / t]]


BK = OrbitType("bk", 1, _b_base, _b_dbase, 24)
CK = OrbitType("ck", 1, _c_base, _c_dbase, 24)
DK = OrbitType("dk", 2, _d_base, _d_dbase, 48)

ALL_TYPES = (A1, A2, A3, BK, CK, DK)
BY_NAME = {t.name: t for t in ALL_TYPES}


def signed_permutations(base):
    """The full octahedral orbit of a point given its |coordinates|."""
    out = []
    seen = set()
    for perm in set(permutations(base)):
        for s0 in (1, -1):
            for s1 in (1, -1):
                for s2 in (1, -1):
                    q = (s0 * perm[0], s1 * perm[1], s2 * perm[2])
                    key = tuple(mp.nstr(v, 30) for v in q)
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append(q)
    return out


def classify(coords, tol=mpf("1e-9")):
    """Identify the orbit type and parameters from one point's |coordinates|.

    Returns (OrbitType, [parameters]) with the coordinates put in canonical
    order for that type.
    """
    v = sorted((abs(c) for c in coords), reverse=True)
    zero = [abs(c) < tol for c in v]
    n_zero = sum(zero)

    if n_zero == 2:
        return A1, []
    if n_zero == 1:
        if abs(v[0] - v[1]) < tol:
            return A2, []
        # (p, q, 0): the free parameter is the smaller nonzero coordinate, so
        # that q = sqrt(1-p^2) reproduces the larger one.
        return CK, [v[1]]
    if abs(v[0] - v[1]) < tol and abs(v[1] - v[2]) < tol:
        return A3, []
    if abs(v[1] - v[2]) < tol:
        return BK, [v[1]]          # (l,l,m) with l the repeated value
    if abs(v[0] - v[1]) < tol:
        return BK, [v[0]]
    return DK, [v[0], v[1]]
