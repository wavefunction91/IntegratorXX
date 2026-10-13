"""Gauss-Newton refinement of an O_h-symmetric angular grid.

The tabulated grids are accurate to about 16 digits, which is the ceiling on
what LebedevLaikov<T> can deliver however wide T is. Newton's method converges
quadratically from that starting point, so three or four iterations reach any
precision you ask for. This is a fundamentally easier problem than finding the
rule in the first place, and needs no globalisation.
"""
import collections

from mpmath import fabs, im, matrix, mp, mpf

from .linalg import lstsq_step
from .moments import conditions, exact_moment, orbit_moment
from .orbits import classify, signed_permutations


class Orbit:
    def __init__(self, otype, params, weight):
        self.type = otype
        self.params = list(params)
        self.weight = weight     # per-point weight

    def __repr__(self):
        ps = ", ".join(mp.nstr(p, 12) for p in self.params)
        return f"<{self.type.name} n={self.type.size} w={mp.nstr(self.weight, 12)} [{ps}]>"


def decompose(points, weights, tol=mpf("1e-9")):
    """Group a tabulated grid into octahedral orbits."""
    groups = collections.OrderedDict()
    for p, w in zip(points, weights):
        key = (mp.nstr(w, 12), tuple(mp.nstr(abs(c), 9) for c in sorted((abs(v) for v in p), reverse=True)))
        groups.setdefault(key, []).append((p, w))

    orbits = []
    for members in groups.values():
        p0, w0 = members[0]
        otype, params = classify(p0, tol)
        if len(members) != otype.size:
            raise ValueError(f"orbit of type {otype.name} has {len(members)} points, expected {otype.size}")
        orbits.append(Orbit(otype, params, w0))
    return orbits


def _in_domain(orbits):
    """Are all orbit bases still real?

    The orbit types carry implicit constraints -- bk is (l, l, sqrt(1-2 l^2)),
    valid only for l <= 1/sqrt(2); ck and dk have their own. A Newton step can
    overshoot one, and mpmath then returns a complex square root, which
    propagates until something tries to order two complex numbers. Backtracking
    on the step is cheaper than reparameterising the orbits.
    """
    for o in orbits:
        for v in o.type.base(o.params):
            if im(v) != 0:
                return False
    return True


def _pack(orbits):
    x = []
    for o in orbits:
        x.append(o.weight)
        x.extend(o.params)
    return x


def _unpack(orbits, x):
    k = 0
    for o in orbits:
        o.weight = x[k]; k += 1
        for j in range(o.type.n_param):
            o.params[j] = x[k]; k += 1


def _residual_and_jacobian(orbits, conds, want_jac):
    n = sum(1 + o.type.n_param for o in orbits)
    F = []
    J = matrix(len(conds), n) if want_jac else None

    # orbit moments and gradients, cached per condition
    for ci, powers in enumerate(conds):
        total = mpf(0)
        col = 0
        for o in orbits:
            m, dm = orbit_moment(o.type, o.params, powers)
            total += o.weight * m
            if want_jac:
                J[ci, col] = m; col += 1
                for k in range(o.type.n_param):
                    J[ci, col] = o.weight * dm[k]; col += 1
            else:
                col += 1 + o.type.n_param
        F.append(total - exact_moment(*powers))
    return F, J


def refine(orbits, order, target_dps, max_iter=12, verbose=False):
    """Refine in place. Returns the residual history."""
    conds = conditions(order)
    hist = []
    tol = mpf(10) ** (-target_dps)
    for it in range(max_iter):
        F, J = _residual_and_jacobian(orbits, conds, want_jac=True)
        r = max(fabs(f) for f in F)
        hist.append(r)
        if verbose:
            print(f"      iter {it}: max|residual| = {mp.nstr(r, 4)}")
        if r < tol:
            break
        d = lstsq_step(J, F)
        base = _pack(orbits)
        # Damped Newton. Staying in the orbits' domain is necessary but not
        # sufficient: on the larger Delley grids the undamped step is simply
        # bad, and the residual grows from 1e-16 to O(1) over a few iterations.
        # Halve until the step both keeps every orbit real and actually reduces
        # the residual.
        scale = mpf(1)
        improved = False
        for _ in range(60):
            _unpack(orbits, [x + scale * d[j] for j, x in enumerate(base)])
            if _in_domain(orbits):
                Fn, _ = _residual_and_jacobian(orbits, conds, want_jac=False)
                if max(fabs(f) for f in Fn) < r:
                    improved = True
                    break
            scale /= 2
        if not improved:
            _unpack(orbits, base)
            break
    return hist


def materialise(orbits):
    """Expand refined orbits into the full point and weight lists."""
    points, weights = [], []
    for o in orbits:
        for q in signed_permutations(o.type.base(o.params)):
            points.append(q)
            weights.append(o.weight)
    return points, weights


def verify(points, weights, order):
    """Worst relative error over all even monomials up to `order`."""
    worst = mpf(0)
    for powers in conditions(order):
        a, b, c = powers
        q = mp.fsum(w * p[0] ** (2 * a) * p[1] ** (2 * b) * p[2] ** (2 * c)
                    for p, w in zip(points, weights))
        e = exact_moment(a, b, c)
        worst = max(worst, fabs(q - e) / e)
    return worst
