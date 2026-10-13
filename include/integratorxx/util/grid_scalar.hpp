#pragma once

#include <cstdlib>
#include <string>

namespace IntegratorXX {
namespace detail {

/**
 *  @brief Convert the decimal digits of a tabulated grid to the working type.
 *
 *  Grid data is stored as decimal strings rather than floating literals: an
 *  unsuffixed literal has type double whatever the template parameter is, so a
 *  literal table cannot carry more precision than double however many digits
 *  are written into it.
 *
 *  The primary template covers extended-precision types that construct from a
 *  string (Boost.Multiprecision, MPFR C++). The built-in types are specialised
 *  onto the corresponding strtoX.
 */
template <typename T>
struct grid_scalar {
  static T parse(const char* s) { return T(s); }
};

template <>
struct grid_scalar<float> {
  static float parse(const char* s) { return std::strtof(s, nullptr); }
};

template <>
struct grid_scalar<double> {
  static double parse(const char* s) { return std::strtod(s, nullptr); }
};

template <>
struct grid_scalar<long double> {
  static long double parse(const char* s) { return std::strtold(s, nullptr); }
};

}  // namespace detail
}  // namespace IntegratorXX
