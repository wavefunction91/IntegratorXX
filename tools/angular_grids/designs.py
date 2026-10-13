"""Matrix-free refinement of equal-weight spherical designs (Womersley grids).

These grids have no symmetry to exploit: every point is its own orbit, so the
unknowns are the 2N tangential degrees of freedom of N points and the Jacobian
has 2N columns. Two consequences shape the implementation.

*Points are carried as unit 3-vectors, not (theta, phi).* The spherical chart is
singular at the poles, and real designs put points there -- the 32-point design
has one at exactly z = 1. In that chart d(p)/d(phi) vanishes, so the Jacobian
acquires a structurally zero column that is a coordinate artifact rather than a
symmetry, and rotations about x and y are not expressible as finite tangent
vectors at all. Stepping in the tangent plane of each point and renormalising
avoids all of this: the parameterisation is regular everywhere on the sphere.

*The Jacobian is rank-deficient by exactly three*, because any rotation of a
spherical design is a spherical design. This needs no special handling: CGLS
started from a zero step keeps every iterate in range(J^T), which is orthogonal
to the null space, so it returns the minimum-norm step -- already the
rotation-projected one. An explicit projection would be redundant here. It is
worth doing when *constructing* a design with a trust-region method, where zero
curvature directions interfere with the model.
"""
from mpmath import mp, mpf, fabs, sqrt

from .moments import _dfact


def exact_monomial(a, b, c):
    """\\int_{S^2} x^a y^b z^c dOmega; zero unless all exponents are even."""
    if a % 2 or b % 2 or c % 2:
        return mpf(0)
    return (4 * mp.pi * _dfact(a - 1) * _dfact(b - 1) * _dfact(c - 1)
            / _dfact(a + b + c + 1))


def conditions(order):
    """A minimal spanning set of monomial conditions for a design of `order`.

    On the unit sphere z^2 = 1 - x^2 - y^2, so any monomial with c >= 2 reduces
    to ones with c in {0,1}. That leaves exactly (order+1)^2 conditions -- the
    dimension of the polynomials of degree <= order restricted to S^2 -- rather
    than the C(order+3,3) of the full monomial set.
    """
    out = []
    for c in (0, 1):
        for d in range(order - c + 1):
            for a in range(d + 1):
                out.append((a, d - a, c))
    return out


def tangent_frame(p):
    """An orthonormal basis of the tangent plane at the unit vector p.

    The seed axis is the one p is least aligned with, so the cross product is
    always well conditioned -- no special case anywhere on the sphere.
    """
    k = min(range(3), key=lambda i: fabs(p[i]))
    seed = [mpf(0)] * 3
    seed[k] = mpf(1)
    e1 = [seed[1] * p[2] - seed[2] * p[1],
          seed[2] * p[0] - seed[0] * p[2],
          seed[0] * p[1] - seed[1] * p[0]]
    n = sqrt(sum(v * v for v in e1))
    e1 = [v / n for v in e1]
    e2 = [p[1] * e1[2] - p[2] * e1[1],
          p[2] * e1[0] - p[0] * e1[2],
          p[0] * e1[1] - p[1] * e1[0]]
    n2 = sqrt(sum(v * v for v in e2))
    return e1, [v / n2 for v in e2]


def retract(p, a, b, e1, e2):
    """Move p by a*e1 + b*e2 and project back onto the sphere."""
    q = [p[i] + a * e1[i] + b * e2[i] for i in range(3)]
    n = sqrt(sum(v * v for v in q))
    return [v / n for v in q]


def _grad(a, b, c, p):
    x, y, z = p
    return (a * x ** (a - 1) * y ** b * z ** c if a else mpf(0),
            b * x ** a * y ** (b - 1) * z ** c if b else mpf(0),
            c * x ** a * y ** b * z ** (c - 1) if c else mpf(0))


class Design:
    """An equal-weight spherical design under refinement."""

    def __init__(self, points, order):
        # Tabulated points are unit vectors only to double precision. Left
        # alone, that ~1e-16 inconsistency between the coordinates and the
        # constraint |p| = 1 corrupts the Jacobian at the same level and stalls
        # the first Newton step. (The (theta, phi) chart normalised implicitly,
        # since every (theta, phi) maps to an exactly unit vector; carrying
        # 3-vectors gives that up and has to do it explicitly.)
        self.points = []
        for p in points:
            n = sqrt(sum(v * v for v in p))
            self.points.append([v / n for v in p])
        self.order = order
        self.n = len(points)
        self.weight = 4 * mp.pi / self.n
        self.conds = conditions(order)
        self._reframe()

    def _reframe(self):
        self.frames = [tangent_frame(p) for p in self.points]
        self._gcache = None

    def _gradients(self):
        """Tangential gradients of every condition at every point.

        The points do not move during a CGLS solve, so these are computed once
        per Newton step rather than once per matrix-vector product. This takes
        every power evaluation out of the inner loop.
        """
        if self._gcache is None:
            tab = []
            for (a, b, c) in self.conds:
                row = []
                for i, p in enumerate(self.points):
                    g = _grad(a, b, c, p)
                    e1, e2 = self.frames[i]
                    row.append((self.weight * sum(g[j] * e1[j] for j in range(3)),
                                self.weight * sum(g[j] * e2[j] for j in range(3))))
                tab.append(row)
            self._gcache = tab
        return self._gcache

    def residual(self):
        F = []
        for (a, b, c) in self.conds:
            s = mp.fsum(p[0] ** a * p[1] ** b * p[2] ** c for p in self.points)
            F.append(self.weight * s - exact_monomial(a, b, c))
        return F

    def Jv(self, v):
        """J @ v, without forming J."""
        out = []
        for row in self._gradients():
            s = mpf(0)
            for i, (d1, d2) in enumerate(row):
                s += v[2 * i] * d1 + v[2 * i + 1] * d2
            out.append(s)
        return out

    def JTu(self, u):
        """J^T @ u, without forming J."""
        out = [mpf(0)] * (2 * self.n)
        for k, row in enumerate(self._gradients()):
            uk = u[k]
            if uk == 0:
                continue
            for i, (d1, d2) in enumerate(row):
                out[2 * i] += uk * d1
                out[2 * i + 1] += uk * d2
        return out

    def rotation_generators(self):
        """The three rotation directions, in tangent-frame coordinates."""
        gens = []
        for axis in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            ax = [mpf(t) for t in axis]
            g = [mpf(0)] * (2 * self.n)
            for i, p in enumerate(self.points):
                v = (ax[1] * p[2] - ax[2] * p[1],
                     ax[2] * p[0] - ax[0] * p[2],
                     ax[0] * p[1] - ax[1] * p[0])
                e1, e2 = self.frames[i]
                g[2 * i] = sum(v[j] * e1[j] for j in range(3))
                g[2 * i + 1] = sum(v[j] * e2[j] for j in range(3))
            gens.append(g)
        return gens

    def step(self, F, iters=None, rtol=None):
        """CGLS for min ||J d + F||; returns the minimum-norm step.

        Every Krylov vector lies in range(J^T), which is orthogonal to null(J),
        so the three rotational directions are never entered and the result is
        the minimum-norm -- already rotation-projected -- step. No explicit
        projection is needed here.

        CG terminates in rank(J) steps only in exact arithmetic, and its rate
        depends on the square root of the condition number rather than on the
        rank, so a cap of rank-many iterations is not enough. Measured on the
        50-point design: capped at 2*2N = 200 the outer iteration stalls at
        2.5e-18, while at 6*2N = 600 the same first step reaches 2.2e-28. The
        cap is therefore generous and the real stopping is by residual, with a
        stagnation guard so a converged solve exits immediately.
        """
        n = 2 * self.n
        iters = iters or 12 * n
        rtol = rtol or mpf(10) ** (-mp.dps + 5)
        d = [mpf(0)] * n
        r = [-f for f in F]
        s = self.JTu(r)
        p = list(s)
        g = mp.fsum(q * q for q in s)
        g0 = g
        stall = 0
        for _ in range(iters):
            q = self.Jv(p)
            qq = mp.fsum(t * t for t in q)
            if qq == 0:
                break
            al = g / qq
            d = [d[j] + al * p[j] for j in range(n)]
            r = [r[k] - al * q[k] for k in range(len(r))]
            s = self.JTu(r)
            gn = mp.fsum(t * t for t in s)
            if gn == 0:
                break
            be = gn / g
            if gn <= rtol * rtol * g0:
                g = gn
                break
            # CG on a rank-deficient system eventually stops making progress;
            # bail rather than grinding out the remaining cap.
            if gn > mpf("0.999999") * g:
                stall += 1
                if stall >= 20:
                    break
            else:
                stall = 0
            g = gn
            p = [s[j] + be * p[j] for j in range(n)]
        return d

    def apply(self, d):
        self.points = [retract(p, d[2 * i], d[2 * i + 1], *self.frames[i])
                       for i, p in enumerate(self.points)]
        self._reframe()


def refine(design, target_dps, max_iter=8, verbose=False):
    """Gauss-Newton with a fully converged inner solve.

    The inner CGLS is run to convergence rather than truncated by a forcing
    term. An Eisenstat-Walker style truncation was tried and its effect could
    not be separated from a Jacobian bug present at the time, so it is simply
    not used: these systems are small enough that solving them properly costs
    little, and correctness of the outer iteration is worth more than inner
    iterations saved.
    """
    hist = []
    tol = mpf(10) ** (-target_dps)
    for it in range(max_iter):
        F = design.residual()
        r = max(fabs(f) for f in F)
        hist.append(r)
        if verbose:
            print(f"      iter {it}: max|residual| = {mp.nstr(r, 4)}", flush=True)
        if r < tol:
            break
        design.apply(design.step(F))
    return hist


def verify(design, order=None):
    """Worst absolute residual over the full monomial set, not the reduced one."""
    order = order or design.order
    worst = mpf(0)
    for d in range(order + 1):
        for a in range(d + 1):
            for b in range(d - a + 1):
                c = d - a - b
                q = mp.fsum(design.weight * p[0] ** a * p[1] ** b * p[2] ** c
                            for p in design.points)
                worst = max(worst, fabs(q - exact_monomial(a, b, c)))
    return worst
