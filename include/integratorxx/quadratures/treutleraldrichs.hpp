#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class TreutlerAldrichs : 
  public Quadrature<TreutlerAldrichs<PointType,WeightType>> {

  using base_type = Quadrature<TreutlerAldrichs<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  TreutlerAldrichs(size_t npts, weight_type R = 1., weight_type alpha = 0.6): 
    base_type( npts, R, alpha ) { }

  TreutlerAldrichs( const TreutlerAldrichs& ) = default;
  TreutlerAldrichs( TreutlerAldrichs&& ) noexcept = default;
};







template <typename PointType, typename WeightType>
struct quadrature_traits<
  TreutlerAldrichs<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, weight_type R, weight_type alpha ) {


    const point_type ln_2 = std::log(2.);
    const point_type pi_ov_npts_p1 = M_PI / (npts + 1);
    
    point_container  points( npts );
    weight_container weights( npts );

    for( size_t i = 0; i < npts; ++i ) {
      const auto xi = std::cos( (i+1) * pi_ov_npts_p1 );

      const auto pow_term  = std::pow( 1. + xi, alpha ) / ln_2;
      const auto log_term  = std::log( (1. - xi)/2. );
      const auto sqrt_term = std::sqrt( (1.+xi)/(1.-xi) );

      points[i]  = -R * pow_term * log_term;
      weights[i] = R * pi_ov_npts_p1 * pow_term * 
        ( sqrt_term - alpha * log_term / sqrt_term );
    } 
    


    return std::tuple( points, weights );

  }

};

}
