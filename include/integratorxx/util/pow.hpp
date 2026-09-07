#pragma once
#include <cmath>
#include <integratorxx/util/type_traits.hpp>

namespace IntegratorXX {

template <typename FloatType, typename IntegralType>
FloatType integer_pow(FloatType x, IntegralType n) {
  if(n == 0) return 1;
  if(n == 1) return x;
  FloatType v = x;
  for(IntegralType i = 1; i < n; ++i) v *= x;
  return v;
}

template <typename FloatType, typename IntegralType>
FloatType half_integer_pow(FloatType x, IntegralType n) {
  assert(n >= 0);
  if(n%2 == 0) return integer_pow(x, n/2);
  else return std::sqrt(x) * integer_pow(x, n/2);
}

}
