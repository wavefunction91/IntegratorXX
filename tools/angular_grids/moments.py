"""Exactness conditions for an O_h-symmetric spherical quadrature.

A rule of algebraic order L integrates every spherical harmonic up to degree L
exactly. For an octahedrally symmetric rule it is equivalent, and much cheaper,
to impose exactness on the even monomials

    x^(2a) y^(2b) z^(2c),   a >= b >= c >= 0,   2(a+b+c) <= L

since the odd moments vanish by symmetry and the monomial conditions span the
same space. The conditions are redundant (x^2+y^2+z^2 = 1 relates them), which
is harmless: the resulting least-squares system is consistent.
"""
from itertools import permutations

from mpmath import mp, mpf


def _dfact(n):
    r = mpf(1)
    while n > 1:
        r *= n
        n -= 2
    return r


def exact_moment(a, b, c):
    """\\int_{S^2} x^(2a) y^(2b) z^(2c) dOmega."""
    return (4 * mp.pi * _dfact(2 * a - 1) * _dfact(2 * b - 1) * _dfact(2 * c - 1)
            / _dfact(2 * (a + b + c) + 1))


def conditions(order, reduced=True):
    """Even monomials whose exactness pins down a rule of the given order.

    With ``reduced`` (the default) only monomials in x and y are used. On the
    unit sphere z^2 = 1 - x^2 - y^2, so every even monomial reduces to a fixed
    linear combination of x^(2a) y^(2b) with constant coefficients; the exact
    integrals obey the same relation, so the residuals do too. Imposing the
    reduced set therefore implies the rest, at roughly a seventh of the cost
    for the larger grids. O_h symmetry lets us also require a >= b.

    Passing ``reduced=False`` returns the full redundant set, which is useful
    as an independent check.
    """
    out = []
    if reduced:
        half = order // 2
        for a in range(half + 1):
            for b in range(min(a, half - a) + 1):
                out.append((a, b, 0))
        return out
    for s in range(order // 2 + 1):
        for a in range(s, -1, -1):
            for b in range(min(a, s - a), -1, -1):
                c = s - a - b
                if 0 <= c <= b:
                    out.append((a, b, c))
    return out


def _distinct_index_perms(base):
    """Index permutations giving distinct coordinate tuples.

    Equal entries in `base` always share the same parameter dependence, so
    keeping one representative per distinct tuple is also correct for the
    derivatives.
    """
    seen = {}
    for perm in permutations(range(3)):
        key = tuple(mp.nstr(base[i], 30) for i in perm)
        seen.setdefault(key, perm)
    return list(seen.values())


def _n_sign(base):
    return 2 ** sum(1 for v in base if v != 0)


def orbit_moment(otype, params, powers):
    """Sum of the monomial over the whole orbit, and its parameter gradient.

    Returns (value, [d/dparam_0, ...]).
    """
    base = otype.base(params)
    dbase = otype.dbase(params)
    ns = _n_sign(base)
    perms = _distinct_index_perms(base)

    total = mpf(0)
    grad = [mpf(0)] * otype.n_param

    for perm in perms:
        # term = prod_j base[perm[j]] ** (2 * powers[j])
        term = mpf(1)
        for j, idx in enumerate(perm):
            e = 2 * powers[j]
            if e:
                term *= base[idx] ** e
        total += term

        for k in range(otype.n_param):
            acc = mpf(0)
            for j, idx in enumerate(perm):
                e = 2 * powers[j]
                if e == 0:
                    continue
                dv = dbase[idx][k]
                if dv == 0:
                    continue
                # d/dparam of base[idx]^e, times the other factors
                rest = mpf(1)
                for j2, idx2 in enumerate(perm):
                    if j2 == j:
                        continue
                    e2 = 2 * powers[j2]
                    if e2:
                        rest *= base[idx2] ** e2
                acc += e * base[idx] ** (e - 1) * dv * rest
            grad[k] += acc

    return ns * total, [ns * g for g in grad]
