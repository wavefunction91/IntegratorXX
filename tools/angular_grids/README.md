# Angular grid generation

Regenerates the Lebedev-Laikov and Delley tables at arbitrary precision.

## Why

The tabulated grids carry about 16 significant digits, so `LebedevLaikov<T>`
cannot be more accurate than `double` however wide `T` is. Worse, the tables
are written as unsuffixed floating literals inside `std::array<T,N>`, and an
unsuffixed literal has type `double` whatever `T` is -- so simply writing more
digits into the existing format would change nothing.

Measured on the 110-point rule, integrating `x^4 y^2 z^2`:

| precision          | current table    | regenerated table |
|--------------------|------------------|-------------------|
| `double`           | 1.91e-15         | 1.74e-15          |
| `long double`      | 1.26e-15         | **2.55e-19**      |
| `cpp_bin_float_50` | not achievable   | **1.94e-40**      |

## How

Finding a spherical quadrature rule and refining one are very different
problems, and only the second is needed here:

* **Construction** is a global optimisation with many local minima. It needs
  second-order globalisation -- this is what a trust-region solver such as
  [OpenTrustRegion](https://github.com/eriksen-lab/opentrustregion) is for --
  but `double` precision suffices, since the goal is only to land in the right
  basin.
* **Refinement** starts from a rule that is already correct to 16 digits and is
  purely local. Newton converges quadratically, so two or three iterations
  reach any precision asked for, with no globalisation at all.

For the existing families construction is already done: the tables *are* the
result. Only refinement is needed, and that is what this tool does.

Tractability comes from symmetry. Both families are invariant under the full
octahedral group, so points fall into orbits carrying one weight and zero, one
or two angular parameters:

| grid              | points | orbits | free parameters |
|-------------------|-------:|-------:|----------------:|
| Lebedev-Laikov 50 |     50 |      4 |               5 |
| Lebedev-Laikov110 |    110 |      6 |              10 |
| Lebedev-Laikov590 |    590 |     20 |              44 |
| Lebedev-Laikov5810|   5810 |    144 |             385 |

So the largest Lebedev rule is a 385-parameter dense Newton solve, not a
17430-parameter one.

The step itself is a rank-revealing modified Gram-Schmidt QR of the Jacobian,
never the normal equations. Two things force that. The Jacobian is not always
full rank -- AB-312 comes out 13 of 16 -- and both `lu_solve` and mpmath's
`qr_solve` refuse a singular matrix outright. And forming `J^T J` squares the
condition number, which turns a solvable system into an unsolvable one: with a
fixed rank threshold on the normal equations, delley_974 stalls at 3.1e-29
instead of reaching 6.9e-51, while tightening the threshold to fix that puts
delley_1454 back to "numerically singular". QR removes the tension rather than
trading the two failures against each other.

A size that needs more working precision than it is given now fails the
generator's verification and is not written. delley_1454 is one: 102 parameters
at that conditioning are not resolvable in 20 digits, and it refines cleanly at
40.

The exactness conditions are imposed on even monomials. On the unit sphere
`z^2 = 1 - x^2 - y^2`, so every even monomial reduces to a fixed linear
combination of `x^(2a) y^(2b)`; imposing only those implies the rest, at about
a seventh of the cost. Results are always *verified* against the full redundant
set.

## Usage

```sh
PYTHONPATH=tools python3 -m angular_grids.generate \
    --family lebedev_laikov --npts 110 --digits 40 --out gen/
```

Requires `mpmath`. Output headers store the digits as decimal strings and
convert them through `detail::grid_scalar<T>::parse` at grid construction,
which works for the built-in types as well as Boost.Multiprecision and MPFR.

## Not covered

* **Ahrens-Beylkin** has icosahedral rather than octahedral symmetry (15012
  points / 60 = 250.2 against 251 distinct weights, the rotation group without
  inversion). Same method, different orbit algebra.
* **Constructing new rules.** That is the trust-region half of the problem;
  see below for what is known about the one table that needs it.

## Equal-weight spherical designs (Womersley)

`designs.py` handles the Womersley family, which has no symmetry to exploit:
every point is its own orbit, so the unknowns are the 2N tangential degrees of
freedom of N points. Two things differ from the octahedral path.

**Points are carried as unit 3-vectors, not (theta, phi).** The spherical chart
is singular at the poles and real designs put points there -- the 32-point
design has one at exactly `z = 1`. In that chart `d(p)/d(phi)` vanishes, so the
Jacobian picks up a structurally zero column that is a coordinate artifact
rather than a symmetry, and rotations about x and y cease to be expressible as
finite tangent vectors. Measured on the 32-point design:

| parameterisation      | zero columns | near-null directions | rotation generators `\|Jg\|/\|J\|\|g\|` |
|-----------------------|-------------:|---------------------:|----------------------------------|
| `(theta, phi)`        |            1 |                    3 | 3.4e-02, 1.3e-16, 2.0e-16        |
| tangent frame         |            0 |                    3 | 2.0e-16, 1.4e-16, 2.9e-16        |

Stepping in each point's tangent plane and renormalising (`retract`) is regular
everywhere, and all three rotations then behave alike -- the residual `1e-16`
being the accuracy of the tabulated input, not of the method.

**The Jacobian is rank-deficient by exactly three**, since any rotation of a
design is a design. This needs no special handling: CGLS started from a zero
step keeps every iterate in `range(J^T)`, orthogonal to the null space, so it
returns the minimum-norm step -- already the rotation-projected one. Explicit
projection would be redundant. It *is* worth doing when constructing a design
with a trust-region method, where zero-curvature directions interfere with the
model and with the Hessian diagonal; `Design.rotation_generators()` returns them
for that purpose.

The conditions use monomials with `c` in `{0,1}`: on the sphere `z^2` reduces,
leaving exactly `(order+1)^2` conditions -- the dimension of the polynomials of
degree <= order on S^2 -- rather than `C(order+3,3)`. For order 7 that is 64
rather than 120.

**Two traps, both measured rather than reasoned about.**

*The input points must be renormalised.* The tabulated coordinates are unit
vectors only to double precision. Carrying them as 3-vectors inherits that
~1e-16 inconsistency with the constraint `|p| = 1`, which corrupts the Jacobian
at exactly the level the refinement is trying to work below -- checked against
finite differences, 8.0e-16 relative error before renormalising and 1.7e-21
after. The `(theta, phi)` chart hid this by normalising implicitly, since every
`(theta, phi)` maps to an exactly unit vector.

*The inner CGLS needs a generous iteration cap.* CG terminates in rank-many
steps only in exact arithmetic; its rate goes as the square root of the
condition number. Capping at `2 * 2N` starves it, the outer Gauss-Newton takes
a bad step, and the result looks convincingly like convergence to a spurious
floor:

| CGLS cap        | womersley_50 outer residuals                          |
|-----------------|-------------------------------------------------------|
| `2 * 2N` = 200  | 5.75e-16, 3.28e-18, 2.78e-18, 2.45e-18 (stalls)       |
| `6 * 2N` = 600  | 5.75e-16, 2.16e-28                                     |

The cap is now `12 * 2N`, with the residual tolerance doing the real stopping
and a stagnation guard so a converged solve exits at once. With both fixed:

```
womersley_32 (order 7):  7.39e-16 -> 3.82e-30 -> 9.18e-41   [107s]
womersley_50 (order 9):  5.75e-16 -> 2.16e-28 -> 4.96e-41   [358s]
```

both reaching the 40-digit working precision floor, verified against the full
monomial set rather than the reduced one.

At 7939 points a single `J*v` product is around 1.3e8 high-precision operations
and the iteration count grows too, so the largest designs want a compiled MPFR
kernel rather than this driver.


## The 552-point Ahrens-Beylkin table

This table is not a valid order-39 rule and cannot be recovered by refinement.
The evidence, with AB-612 as a control throughout:

**Its points admit no valid weighting.** The exactness conditions are *linear*
in the weights, so the best possible weights for a fixed set of points follow
from one least-squares solve:

| grid                | best-possible weights, full order-39/41 conditions |
|---------------------|----------------------------------------------------|
| AB-612 (known good) | 8.3e-14, all positive                              |
| AB-552              | **2.4e-02, five negative weights**                 |

So the point *positions* are wrong, not just the weights, and not merely
imprecise. This is the decisive result.

**Refinement diverges.** Gauss-Newton from the tabulated data goes
0.4619 -> 0.03831 -> 145.9. Levenberg-Marquardt globalises the step and does
make progress -- 0.4619 -> 9.99e-4 -> 7.75e-6 -> 1.19e-7 -- but plateaus there,
five orders short, and on a condition subset that does not determine all 28
parameters.

**What is not evidence.** The Jacobian's rank over low-degree conditions is 6
against 28 parameters, which looks alarming until AB-612 is measured the same
way and gives rank 6 against 31. Few icosahedral invariants exist at low
degree; this says nothing about either grid.

The remaining routes are to obtain the original Ahrens-Beylkin data
(Proc. R. Soc. A 465, 3103) and refine that, or to construct the rule from
scratch with a trust-region method. Of the two the first is much cheaper, and
the other 55 sizes in the family are unaffected -- they all refine cleanly.

## Regenerating everything

```sh
PYTHONPATH=tools python3 -m angular_grids.batch \
    --family delley --digits 40 --out gen/ --jobs 8
PYTHONPATH=tools python3 -m angular_grids.switch_to_load --family delley
```

`batch` reads the sizes from the family's own `algebraic_order_by_npts` table,
so the two cannot drift apart, and skips rather than writes any size that
misses the requested precision. `switch_to_load` then repoints the dispatch at
the regenerated tables and drops the `4*M_PI` scaling loop for the families
whose shipped tables are normalised to one.
