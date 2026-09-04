#pragma once

#include <cstddef>

#ifdef ENABLE_STRING_REALS
#include <string_view>
#endif

/** @file config.hpp
 *
 *  Build-time policy for how numeric literals are spelled and carried.
 *
 *  Quadrature code distinguishes two kinds of literal:
 *
 *  - **Integral** values, which are exactly representable in every arithmetic
 *    type the library supports. These are written `IXX_INT(2)` and converted
 *    with `fp_traits<T>::from_integer`. Note that a value is integral by its
 *    *value*, not its spelling: a `2.0` in the source is an integral literal.
 *
 *  - **Non-integral** values, which are not exactly representable. These are
 *    written `IXX_REAL(0.6931471805599453) `and converted with
 *    `fp_traits<T>::from_real`.
 *
 *  Keeping the two apart lets arithmetic be carried in exact integers for as
 *  long as possible, converting once at the end, and lets a ratio of two
 *  integers go through `fp_traits<T>::divide_integer`, which is a single
 *  correctly-rounded operation rather than a pre-rounded constant.
 *
 *  When `ENABLE_STRING_REALS` is defined, `IXX_REAL` captures the *source text*
 *  of its argument instead of a `double`. A `double` literal has already lost
 *  precision by the time any code can inspect it, so a type more precise than
 *  `double` -- or one that must bound its own error -- cannot recover the
 *  intended value from it. Carrying the decimal text lets such a type parse it
 *  directly. See `fp_traits::from_real`.
 */

namespace IntegratorXX {

/** @brief Integer type for counts, indices, and integral literals.
 *
 *  Signed on purpose. Intermediate expressions routinely go negative
 *  (`2*i - 1`, `npts - 1`, `M - 1`), and an unsigned type would wrap silently
 *  rather than produce a negative value.
 */
using ixx_int = std::ptrdiff_t;

#ifdef ENABLE_STRING_REALS
/** @brief Non-integral literals carried as their exact decimal source text.
 *
 *  `std::string_view` rather than `std::string`: the solid-angle tables are
 *  `static constexpr`, which requires a literal type.
 */
using ixx_real = std::string_view;
#else
/// Non-integral literals carried as `double`.
using ixx_real = double;
#endif

} // namespace IntegratorXX

/** @brief Denote an integral literal or an integral expression.
 *
 *  Parenthesised because expressions are passed, e.g. `IXX_INT(2*i - 1)`.
 *
 *  @warning The cast applies to the *result*, not the operands. If the operands
 *  are unsigned, `IXX_INT(2*i - 1)` evaluates in unsigned arithmetic and only
 *  then converts, so an underflow has already happened. Convert the count once
 *  at the top of a scope (`const ixx_int n = IXX_INT(npts);`) and derive the
 *  rest from it.
 */
#define IXX_INT(x) (static_cast<::IntegratorXX::ixx_int>(x))

#ifdef ENABLE_STRING_REALS
#  define IXX_REAL(x) (::IntegratorXX::ixx_real(#x))
#else
/** @brief Denote a non-integral literal.
 *
 *  Deliberately expands to its argument unchanged rather than to `(x)`. The
 *  argument is by contract a single numeric literal token, so no parentheses
 *  are needed for precedence, and the omission makes the default build's token
 *  stream byte-identical to one written with bare literals -- which matters
 *  when the expansion happens ~1.9 million times across the tables.
 */
#  define IXX_REAL(x) x
#endif

namespace IntegratorXX {

/// pi. Spelled as digits because `IXX_REAL(M_PI)` would stringify to "M_PI".
inline constexpr ixx_real ixx_pi = IXX_REAL(3.14159265358979323846);

} // namespace IntegratorXX
