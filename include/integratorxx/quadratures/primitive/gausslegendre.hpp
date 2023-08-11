#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/util/legendre.hpp>

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
    const auto eps = std::numeric_limits<double>::epsilon();

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
        bool converged = false;
        const int maxit = 100;
        for(int it = 0; it < maxit; ++it) {
          // Evaluate the Legendre polynomial at z and its derivative
          std::tie(p_n, dp_n, std::ignore) = eval_Pn(z,npts);

          // Newton update for root
          z_old = z;
          z  -= p_n / dp_n;

          // Convergence check
          if(std::abs(z-z_old) <= eps) {
            converged = true;
            break;
          }
        } // end while

        if(not converged) {
          throw std::runtime_error(
            "Gauss-Legendre Newton Iterations Failed to Converge"
          );
        }

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
