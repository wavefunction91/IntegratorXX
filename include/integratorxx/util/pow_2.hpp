#pragma once
#include <integratorxx/util/type_traits.hpp>

namespace IntegratorXX {

template <typename IntegralType>
detail::enable_if_integral_t<IntegralType, IntegralType> 
  integer_pow_2(IntegralType n) {
  assert(n >= 0);
  return IntegralType(1) << n;
}

template <typename FloatType, typename IntegralType>
FloatType half_integer_pow_2(IntegralType n) {
  if(n < 0) return FloatType(1) / half_integer_pow_2<FloatType>(-n);
  else if(n%2 == 0) return integer_pow_2(n/2);
  else return M_SQRT2 * integer_pow_2(n/2);
}

}
