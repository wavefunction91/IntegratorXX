#pragma once

#include <integratorxx/config.hpp>
#include <type_traits>

#ifdef ENABLE_STRING_REALS
#include <charconv>
#endif

namespace IntegratorXX {

/**
 *  @brief Customization point for the floating-point operations used by the
 *         quadratures.
 *
 *  IntegratorXX is templated on the type used to represent quadrature points
 *  and weights. User-supplied types do not in general share
 *  the semantics of the built-in floating point types. This trait is the single
 *  seam through which every operation whose meaning depends on the type is
 *  routed, so that supporting a new type is a matter of specializing one class
 *  rather than auditing every quadrature.
 *
 *  Where an operation can be reached by argument-dependent lookup the primary
 *  template performs a two-step (`using std::fn; fn(x);`) call, so such a type
 *  works with **no specialization at all**. Specialize only when ADL is not
 *  sufficient, or when the default behaviour is wrong for the type -- see
 *  `from_real` below for the motivating case.
 *
 *  @tparam T      The floating-point-like type.
 *  @tparam Enable Hook for SFINAE-constrained partial specializations.
 */
template <typename T, typename Enable = void>
struct fp_traits {
  /** @brief Convert an integral value to @p T.
   *
   *  Integral values are exactly representable in every supported arithmetic
   *  type, so this conversion introduces no error and types that model a set of
   *  values must map it to a degenerate one.
   *
   *  @param[in] v An integral value, normally declared `IXX_INT(...)`.
   *  @return    @p v as a @p T.
   */
  static T from_integer(ixx_int v) { return T(v); }

  /** @brief Convert a non-integral literal to @p T.
   *
   *  The argument is normally declared `IXX_REAL(...)`. Its type depends on the
   *  build: a `double` by default, or the literal's decimal source text when
   *  `ENABLE_STRING_REALS` is defined.
   *
   *  @warning In string mode this default parses via `double`, which discards
   *  exactly the precision the string form exists to preserve. That is the best
   *  a generic implementation can do, and it keeps `float`/`double` correct,
   *  but it means **string mode only pays off for types that specialize this
   *  function**. A type that bounds its own error should parse the text
   *  directly -- for an interval type, twice under directed rounding -- to
   *  obtain a tight enclosure of the decimal rather than of an already-rounded
   *  `double`.
   *
   *  @param[in] v The literal, as `ixx_real`.
   *  @return    @p v as a @p T.
   */
  static T from_real(ixx_real v) {
#ifdef ENABLE_STRING_REALS
    double d{};
    std::from_chars(v.data(), v.data() + v.size(), d);
    return T(d);
#else
    return T(v);
#endif
  }

  /** @brief Convert the exact rational @p num / @p den to @p T.
   *
   *  Prefer this over a pre-divided literal wherever a constant is the ratio of
   *  two integers. Both operands are exact, so the result of the division is
   *  the only rounding. Converting a constant that was already rounded
   *  to `double` by the compiler introduces additional error.
   *
   *  @param[in] num The numerator.
   *  @param[in] den The denominator.
   *  @return    @p num / @p den as a @p T.
   */
  static T divide_integer(ixx_int num, ixx_int den) {
    return from_integer(num) / from_integer(den);
  }
};

}  // namespace IntegratorXX
