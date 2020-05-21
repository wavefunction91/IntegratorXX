#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class GaussChebyshev1 : 
  public Quadrature<GaussChebyshev1<PointType,WeightType> {

  using base_type = typename 
    Quadrature<GaussChebyshev1<PointType,WeightType>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  GaussChebyshev1(size_t npts, point_type lo, point_type up):
    base_type( npts, lo, up ) { }

  GaussChebyshev1( const GaussChebyshev1& ) = default;
  GaussChebyshev1( GaussChebyshev1&& ) noexcept = default;
};







template <typename PointType, typename WeightType>
struct quadrature_traits<
  GaussChebyshev1<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, point_type lo, point_type up ) {

    point_container  points( npts );
    weight_container weights( npts );

    const weight_type pi_ov_npts    = M_PI / npts;
    const weight_type two_npts_x_pi = 2 * npts * M_PI;

    #pragma omp parallel for
    for( size_t i = 0; i < npts; ++i ) {
      // Generate raw points and wghts on (-1,1)
      point_type pt = std::cos( (2.0*(i+1)-1.) / two_npts_x_pi );

      // As we're integrating f(x) not \frac{f(x)}{\sqrt{1 - x^2}}
      // we must factor the \sqrt{1-x^2} into the wghts
      weight_type wgt = pi_ov_npts * std::sqrt( 1. - pt * pt ); 

      // TODO: Tranform on bound
      points[i]  = pt;
      weights[i] = wgt;
    }

    return std::tuple( points, weights );

  }

}

}
