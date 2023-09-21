#pragma once

#include <integratorxx/quadrature.hpp>
#include <vector>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class UniformTrapezoid : 
  public Quadrature<UniformTrapezoid<PointType,WeightType>> {

  using base_type = Quadrature<UniformTrapezoid<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  UniformTrapezoid(size_t npts, point_type lo, point_type up): 
    base_type( npts, lo, up ) { }

  UniformTrapezoid( const UniformTrapezoid& ) = default;
  UniformTrapezoid( UniformTrapezoid&& ) noexcept = default;
};







template <typename PointType, typename WeightType>
struct quadrature_traits<
  UniformTrapezoid<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static constexpr bool bound_inclusive = true;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, point_type lo, point_type up ) {

    const auto npts_use = npts - 1;
    const auto delta_x  = (up - lo) / npts_use;
    point_container  points( npts );
    weight_container weights( npts, delta_x );

    for( size_t i = 0; i <= npts_use; ++i )
      points[i] = lo + i * delta_x;
    
    weights.front() *= 0.5;
    weights.back()  *= 0.5;

    return std::make_tuple( points, weights );

  }

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts ) { return generate(npts, 0.0, 1.0); }
};

}
