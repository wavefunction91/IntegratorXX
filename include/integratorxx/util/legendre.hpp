#pragma once
#include <tuple>

namespace IntegratorXX {

/**
 *  Evaluates the Legendre polynomial at z with the recurrence
 *  relation \f$(n+1) P_{n+1}(x) = (2n+1) x P_{n}(x) - n P_{n-1}(x)
 *  \f$, and its derivative with the relation \f$ \frac {x^2-1} {n}
 *  \frac {\rm d} {{\rm d}x} P_{n}(x) = x P_{n}(x) - P_{n-1}(x) \f$.
 *
 *  Inputs:
 *  x - x coordinate
 *  n - order of polynomial
 *
 *  Outputs: pair (p_n, dp_n)
 *  p_n - value of the Legendre polynomial, P_{n}(x)
 *  dp_n - value of the derivative of the Legendre polynomial, dP_{n}/dx
 */
template <typename PointType> inline std::tuple<PointType,PointType,PointType> eval_Pn(PointType x, size_t n) {
  // Evaluate P_{n}(x) by iteration starting from P_0(x) = 1
  PointType p_n = 1.0;
  // Values of P_{n-1}(x) and P_{n-2}(x) used in the recursions
  PointType p_n_m1(0.0), p_n_m2(0.0);

  // Recurrence can be rewritten as
  // P_{n}(x) = [(2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)] / n
  for( size_t m = 1; m <= n; ++m ){
    p_n_m2 = p_n_m1;
    p_n_m1 = p_n;
    p_n = ((2.0 * m - 1.0) * x * p_n_m1 - (m-1.0)*p_n_m2)/m;
  }

  // Compute the derivative
  // dP_{n}(x) = (x P_{n}(x) - P_{n-1}(x)) * n/(x^2-1)
  PointType dp_n = (x * p_n - p_n_m1) * (n / (x * x - 1.0));

  return std::make_tuple( p_n, dp_n, p_n_m1 );
}

}
