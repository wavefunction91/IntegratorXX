#pragma once
#include <cmath>
#include <integratorxx/util/factorial.hpp>
#include <integratorxx/util/constants.hpp>
#include <integratorxx/util/pow_2.hpp>

namespace IntegratorXX {

/// GAMMA(n) = (n-1)!
template <typename FloatType, typename IntegralType>
FloatType integer_tgamma(IntegralType n) {
  if(n <= 0) return std::numeric_limits<FloatType>::quiet_NaN();
  else {
    auto val = factorial(n-1);
    assert(val <= std::numeric_limits<FloatType>::max());
    return val;
  }
}

/// GAMMA(n/2) = 
///   CASE n even -> INTEGER_GAMMA(n/2)
///   CASE n odd  -> SQRT_PI / POW_2((n-1)/2) * (n-2)!!
template <typename FloatType, typename IntegralType>
FloatType half_integer_tgamma(IntegralType n) {
  if(n%2 == 0) return integer_tgamma<FloatType>(n/2);
  else {
    constexpr auto c = constants::sqrt_pi_v<FloatType>;
    const auto df = double_factorial(n-2);
    const auto pt = half_integer_pow_2<FloatType>(n-1);
    return  c * df / pt;
  }
}


}
