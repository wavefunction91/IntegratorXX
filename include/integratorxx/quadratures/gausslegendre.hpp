#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/util/bound_transform.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class GaussLegendre :
  public Quadrature<GaussLegendre<PointType,WeightType>> {

  using base_type = Quadrature<GaussLegendre<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussLegendre(size_t npts):
    base_type( npts ) { }

  GaussLegendre( const GaussLegendre& ) = default;
  GaussLegendre( GaussLegendre&& ) noexcept = default;
};

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
  template <typename PointType> inline std::pair<PointType,PointType> eval_Pn(PointType x, size_t n) {
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

    return std::make_pair( p_n, dp_n );
  }


template <typename PointType, typename WeightType>
struct quadrature_traits<
  GaussLegendre<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
  generate( size_t npts ) {
    point_container  points( npts );
    weight_container weights( npts );

    // Absolute precision for the nodes
    const double eps=3.0e-11;

    // Since the rules are symmetric around the origin, we only need
    // to compute one half of the points
    const size_t mid = (npts + 1) / 2;
    for( size_t idx = 0; idx < mid; ++idx ) {
        // Index
        const point_type i = idx+1;

        // Standard initial guess for location of i:th root in [-1, 1]
        point_type z = std::cos( M_PI * (i - 0.25) / (npts + 0.5));
        // Old value of root
        point_type z_old = -2.0;
        // Values of P_n(x) and dP_n/dx
        point_type p_n, dp_n;

        // Solve the root to eps absolute precision
        while(std::abs(z-z_old) > eps) {
          // Evaluate the Legendre polynomial at z and its derivative
          auto pn_eval = eval_Pn(z,npts);
          p_n = pn_eval.first;
          dp_n = pn_eval.second;

          // Newton update for root
          z_old = z;
          z  = z_old - p_n / dp_n;
        } // end while

        // Quadrature node is z
        point_type  pt = z;
        // Quadrature weight is 2 / [ (1-x^2) dPn(x)^2]
        weight_type wgt = 2. / ((1. - z*z) * dp_n*dp_n);

        // Store the symmetric points
        points[idx] = -pt;
        weights[idx] = wgt;

        points[npts-1-idx]  = pt;
        weights[npts-1-idx] = wgt;
    } // Loop over points

    return std::make_tuple( points, weights );
  }

};

}
