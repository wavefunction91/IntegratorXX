"""Small dense least squares shared by the refinement paths.

The Gauss-Newton systems here are tiny but not always full rank, and mpmath's
lu_solve and qr_solve both refuse a singular matrix outright rather than
returning something usable.

Forming the normal equations and then dropping ill-conditioned columns does not
work either: the two failure modes trade off against each other. With a fixed
1e-16 rank threshold, delley_1454 solves at 20 digits but delley_974 stalls at
3.1e-29 instead of 6.9e-51, because columns that still carry information are
being discarded. Tying the threshold to the working precision fixes 974 and
puts 1454 back to "numerically singular". The tension is not in the threshold
-- it is that J^T J squares the condition number.

So the step is solved by a rank-revealing modified Gram-Schmidt QR of J itself.
Columns that are genuinely unresolvable are skipped and left at zero; the rest
are solved without ever squaring anything.
"""
from mpmath import matrix, mp, mpf, sqrt


def _mgs(J, tol):
    """Rank-revealing modified Gram-Schmidt.

    Returns (kept, Q, R) with Q's columns orthonormal and R upper triangular
    over the kept columns, so that J[:, kept] = Q R.
    """
    kept, Q, R = [], [], []
    for j in range(J.cols):
        v = [J[i, j] for i in range(J.rows)]
        nrm0 = sqrt(mp.fsum(x * x for x in v))
        if nrm0 == 0:
            continue
        col = []
        for q in Q:
            d = mp.fsum(v[i] * q[i] for i in range(len(v)))
            col.append(d)
            v = [v[i] - d * q[i] for i in range(len(v))]
        nrm = sqrt(mp.fsum(x * x for x in v))
        if nrm <= tol * nrm0:
            continue                      # not resolvable against what we have
        Q.append([x / nrm for x in v])
        col.append(nrm)
        R.append(col)                     # length len(Q) after the append
        kept.append(j)
    return kept, Q, R


def _default_tol():
    """Cut only what is unresolvable at the working precision.

    QR does not square the condition number, so it carries columns the normal
    equations could not, and the threshold should discard only genuinely null
    directions. Cutting at a fixed 1e-14 instead costs delley_974 its quadratic
    convergence: 11 iterations to 1.97e-33 rather than 2 to 6.93e-51.

    Cases that need more working precision than they are given now fail the
    generator's final verification rather than crashing or diverging silently.
    delley_1454 is one: it refines cleanly at 40 digits and cannot be solved at
    20, whatever the threshold.
    """
    return mpf(10) ** (-(mp.dps - 5))


def independent_columns(J, tol=None):
    """Indices of a maximal independent set of columns."""
    return _mgs(J, tol if tol is not None else _default_tol())[0]


def lstsq_step(J, F, tol=None):
    """Gauss-Newton step for min ||J d + F||, without forming J^T J.

    Directions the Jacobian cannot resolve are left at zero: the conditions
    cannot see them, so moving along them is guesswork.
    """
    # Column equilibration. A weight column and an angular column differ by
    # orders of magnitude here, and that disparity alone wrecks the
    # conditioning: on delley_2354 the undamped step is unusable and damping
    # then crawls, a factor of two per iteration instead of squaring. Scaling
    # every column to unit norm and unscaling the step afterwards costs
    # nothing and is the difference between converging and not.
    scale = []
    Js = matrix(J.rows, J.cols)
    for j in range(J.cols):
        c = sqrt(mp.fsum(J[i, j] ** 2 for i in range(J.rows)))
        scale.append(c if c != 0 else mpf(1))
        for i in range(J.rows):
            Js[i, j] = J[i, j] / scale[j]

    kept, Q, R = _mgs(Js, tol if tol is not None else _default_tol())
    m = len(kept)
    if m == 0:
        return [mpf(0)] * J.cols

    b = [-F[i] for i in range(J.rows)]
    y = [mp.fsum(q[i] * b[i] for i in range(len(b))) for q in Q]

    # back-substitute R x = y; R[j] holds column j of R, entries 0..j
    x = [mpf(0)] * m
    for j in range(m - 1, -1, -1):
        acc = y[j]
        for k in range(j + 1, m):
            acc -= R[k][j] * x[k]
        x[j] = acc / R[j][j]

    step = [mpf(0)] * J.cols
    for a, j in enumerate(kept):
        step[j] = x[a] / scale[j]          # undo the equilibration
    return step
