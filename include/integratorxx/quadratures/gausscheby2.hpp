#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class GaussChebyshev1 : 
  public Quadrature<GaussChebyshev1<PointType,WeightType>> {

  using base_type = Quadrature<GaussChebyshev1<PointType,WeightType>>;

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

    const weight_type pi_ov_npts_p_1 = M_PI / (npts+1);

    weight_container weights( npts, pi_ov_npts_p_1 );


    point_container  points( npts );
    for( size_t i = 0; i < npts; ++i ) {
      points[i] = std::cos( (i+1) * pi_ov_npts_p_1 ); 
    }

    return std::make_tuple( points, weights );

  }

};

}
