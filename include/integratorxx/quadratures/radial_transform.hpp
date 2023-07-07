#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

template <typename BaseQuad, typename RadialTraits>
struct RadialTransformQuadrature :
  public Quadrature<RadialTransformQuadrature<BaseQuad, RadialTraits>> {

  using base_type = Quadrature<RadialTransformQuadrature<BaseQuad, RadialTraits>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  RadialTransformQuadrature(size_t npts, weight_type R = 1.0) :
    base_type( npts, R ) { }

  RadialTransformQuadrature( const RadialTransformQuadrature& )     = default;
  RadialTransformQuadrature( RadialTransformQuadrature&& ) noexcept = default;

};

template <typename BaseQuad, typename RadialTraits>
struct quadrature_traits<RadialTransformQuadrature<BaseQuad, RadialTraits>> {

  using point_type       = typename BaseQuad::point_type;
  using weight_type      = typename BaseQuad::weight_type;
  using point_container  = typename BaseQuad::point_container;
  using weight_container = typename BaseQuad::weight_container;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, weight_type R ) {

    using base_quad_traits = quadrature_traits<BaseQuad>;

    point_container  points( npts );
    weight_container weights( npts );


    const auto npts_base = base_quad_traits::bound_inclusive ? npts+2 : npts;
    auto [base_x, base_w] = base_quad_traits::generate(npts_base);

    const auto ipts_offset = !!base_quad_traits::bound_inclusive;
    for(size_t i = 0; i < npts; ++i) {
      const auto xi = base_x[i + ipts_offset];
      const auto wi = base_w[i + ipts_offset];
      points[i]  = R * RadialTraits::radial_transform(xi);
      weights[i] = R * wi * RadialTraits::radial_jacobian(xi);
    }

    return std::make_tuple(points, weights);
  }

};
  


}
