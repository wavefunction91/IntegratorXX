#pragma once
#include <cmath>
#include <iostream>

namespace IntegratorXX {

// Newton-Raphson implementation of Lambert-W function for real arguments
template <typename RealType>
RealType lambert_w_newton(RealType x, RealType w0) {
  
  const auto eps = std::numeric_limits<RealType>::epsilon();
  const size_t max_iter = 100000;
  RealType w;
  size_t i;
  for(i = 0; i < max_iter; ++i) {
    const auto exp_term = std::exp(w0);
    const auto w_exp_term = w0 * exp_term;
    w = w0 - (w_exp_term - x) / (exp_term + w_exp_term);
    if(std::abs(w - w0) < eps) break;
    w0 = w;
  }

  return w;
}

template< typename RealType>
RealType lambert_w0(RealType x) {
  const RealType one = 1;
  const auto e = std::exp(one);
  if(x < -one/e)       return std::numeric_limits<RealType>::quiet_NaN();
  if(x == RealType(0)) return 0.0;
  if(x == -one/e)      return -1.0;
  
  RealType w0;
  if(x > e) {
    w0 = std::log(x) - std::log(std::log(x));
  } else if(x > 0.0) {
    w0 = x / e;
  } else {
    const auto xe = x * e;
    const auto sqrt_term = std::sqrt(1.0 + xe);
    w0 = xe * std::log(1.0 + sqrt_term) / (1.0 + xe + sqrt_term);
  }

  return lambert_w_newton(x, w0);

}

}
