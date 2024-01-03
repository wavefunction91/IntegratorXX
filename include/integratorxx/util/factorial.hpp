#pragma once
#include <array>
#include <limits>

#include <integratorxx/util/type_traits.hpp>

namespace IntegratorXX {

namespace detail {

// Factorials for [0,20]
static constexpr std::array<uint64_t,21> cached_factorials = {
  1ul, 1ul, 2ul, 6ul, 24ul, 120ul, 720ul, 5040ul, 40320ul, 362880ul, 3628800ul, 39916800ul,
  479001600ul, 6227020800ul, 87178291200ul, 1307674368000ul, 20922789888000ul,
  355687428096000ul, 6402373705728000ul, 121645100408832000ul, 2432902008176640000ul
};

// Double factorials for [0,32]
static constexpr std::array<uint64_t,33> cached_double_factorials = {
  1ul, 1ul, 2ul, 3ul, 8ul, 15ul, 48ul, 105ul, 384ul, 945ul, 3840ul, 10395ul, 46080ul, 135135ul,
  645120ul, 2027025ul, 10321920ul, 34459425ul, 185794560ul, 654729075ul, 3715891200ul,
  13749310575ul, 81749606400ul, 316234143225ul, 1961990553600ul, 7905853580625ul,
  51011754393600ul, 213458046676875ul, 1428329123020800ul, 6190283353629375ul,
  42849873690624000ul, 191898783962510625ul, 1371195958099968000ul
};

}

template <typename IntegralType>
detail::enable_if_integral_t<IntegralType, IntegralType> 
  factorial(IntegralType n) {

  assert(n >= 0);
  if( n < detail::cached_factorials.size() ) {
    auto val = detail::cached_factorials[n];
    assert(val <= std::numeric_limits<IntegralType>::max());
    return val;
  } else {
    auto val = n * factorial(n-1);
    assert(val <= std::numeric_limits<IntegralType>::max());
    return val;
  }

}

template <typename IntegralType>
detail::enable_if_integral_t<IntegralType, IntegralType> 
  double_factorial(IntegralType n) {

  // If we're signed, check negative arguments
  if(std::is_signed_v<IntegralType>) {
    // Only provide values for which n!! is integral
    if(n == -1) return 1;
    if(n == -3) return -1;
  }

  // This also handles the case when n < 0 is undefied or non-integral
  assert(n >= 0);

  if( n < detail::cached_double_factorials.size() ) {
    auto val = detail::cached_double_factorials[n];
    assert(val <= std::numeric_limits<IntegralType>::max());
    return val;
  } else {
    auto val = n * double_factorial(n-2);
    assert(val <= std::numeric_limits<IntegralType>::max());
    return val;
  }

}

}
