"""Icosahedral orbit algebra for the Ahrens-Beylkin grids.

Unlike Lebedev-Laikov and Delley, the AB grids are invariant under the
icosahedral *rotation* group I (order 60), without inversion -- which is why
they are not antipodally symmetric and why 15012 points divide as 15012/60 =
250.2 against 251 distinct weights.

Orbits under I:

    12 pts   icosahedron vertices        0 parameters
    20 pts   dodecahedron vertices       0 parameters
    30 pts   edge midpoints              0 parameters
    60 pts   generic position            2 parameters

The group is built numerically rather than from memorised generators: it is
exactly the set of rotations carrying the icosahedron's vertex set to itself,
and each is fixed by naming where one vertex and one of its neighbours go.
"""
from mpmath import mp, mpf, sqrt, matrix


def _icosahedron_vertices():
    """The 12 unit vertices: cyclic permutations of (0, +-1, +-phi)."""
    phi = (1 + sqrt(5)) / 2
    n = sqrt(1 + phi * phi)
    out = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            a, b = mpf(s1) / n, mpf(s2) * phi / n
            out.append((mpf(0), a, b))
            out.append((b, mpf(0), a))
            out.append((a, b, mpf(0)))
    return out


def _frame(v, n):
    """Right-handed orthonormal frame from a vertex and one of its neighbours."""
    e1 = list(v)
    d = sum(n[i] * e1[i] for i in range(3))
    e2 = [n[i] - d * e1[i] for i in range(3)]
    m = sqrt(sum(x * x for x in e2))
    e2 = [x / m for x in e2]
    e3 = [e1[1] * e2[2] - e1[2] * e2[1],
          e1[2] * e2[0] - e1[0] * e2[2],
          e1[0] * e2[1] - e1[1] * e2[0]]
    return e1, e2, e3


def rotation_group(tol=mpf("1e-20")):
    """The 60 rotation matrices of the icosahedral group."""
    V = _icosahedron_vertices()
    # nearest neighbours: the 5 vertices at minimum positive distance
    def neighbours(v):
        ds = []
        for w in V:
            d = sum((v[i] - w[i]) ** 2 for i in range(3))
            if d > tol:
                ds.append((d, w))
        ds.sort(key=lambda t: t[0])
        dmin = ds[0][0]
        return [w for d, w in ds if d < dmin * mpf("1.5")]

    v0 = V[0]
    n0 = neighbours(v0)[0]
    F0 = _frame(v0, n0)
    # R maps the reference frame onto (v, n): R = F(v,n)^T . F0  in row form
    G = []
    for v in V:
        for n in neighbours(v):
            F = _frame(v, n)
            R = matrix(3, 3)
            for i in range(3):
                for j in range(3):
                    R[i, j] = mp.fsum(F[k][i] * F0[k][j] for k in range(3))
            G.append(R)
    return G


_GROUP = None


def group():
    global _GROUP
    if _GROUP is None:
        _GROUP = rotation_group()
    return _GROUP


def orbit(p, tol=mpf("1e-9")):
    """The icosahedral orbit of a point, de-duplicated.

    Duplicates are found by distance, not by formatting the coordinates: a
    coordinate that should be zero comes out as +-1e-31, and two such values
    have different decimal representations while being the same point.

    The default tolerance is loose on purpose. Orbits are usually generated
    from tabulated points that are only accurate to double precision, so
    images that should coincide differ by ~1e-16; a tolerance tighter than the
    input accuracy silently splits a 12-point orbit into 60 near-duplicates.
    Distinct points in these grids are separated by far more than 1e-9.
    """
    out = []
    for R in group():
        q = tuple(mp.fsum(R[i, j] * p[j] for j in range(3)) for i in range(3))
        if not any(sum((q[i] - r[i]) ** 2 for i in range(3)) < tol * tol for r in out):
            out.append(q)
    return out


def special_positions(tol=mpf("1e-9")):
    """Exact representatives of the three special orbits.

    A 12-, 20- or 30-point orbit has no free parameter: its position is fixed
    by the symmetry. Taking the representative from tabulated double-precision
    data instead makes it a *generic* point whose 60 images cluster in fives
    (or threes, or twos) rather than coinciding, and de-duplication then keeps
    an arbitrary member of each cluster -- injecting ~1e-16 of asymmetry into
    a grid that should be exactly symmetric.

    Returns {12: vertex, 20: face centre, 30: edge midpoint}.
    """
    V = _icosahedron_vertices()

    def norm(p):
        m = sqrt(sum(x * x for x in p))
        return tuple(x / m for x in p)

    def dist2(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3))

    # nearest-neighbour distance defines adjacency
    ds = sorted(dist2(V[0], w) for w in V[1:])
    dmin = ds[0]
    adj = lambda a, b: dist2(a, b) < dmin * mpf("1.5")

    edge = None
    face = None
    for i, a in enumerate(V):
        for b in V[i + 1:]:
            if not adj(a, b):
                continue
            if edge is None:
                edge = norm([a[k] + b[k] for k in range(3)])
            for c in V:
                if c is a or c is b or not (adj(a, c) and adj(b, c)):
                    continue
                if face is None:
                    face = norm([a[k] + b[k] + c[k] for k in range(3)])
                break
            if face is not None and edge is not None:
                break
        if face is not None and edge is not None:
            break
    return {12: V[0], 20: face, 30: edge}
