#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

// Base type for all radial traits
struct RadialTraits {
  virtual ~RadialTraits() noexcept = default;
};

template <typename RadialTraitsType>
const RadialTraitsType& radial_traits_cast( const RadialTraits& traits ) {
  return dynamic_cast<const RadialTraitsType&>(traits);
}

template <typename BaseQuad, typename RadialTraitsType>
struct RadialTransformQuadrature :
  public Quadrature<RadialTransformQuadrature<BaseQuad, RadialTraitsType>> {

  using base_type = Quadrature<RadialTransformQuadrature<BaseQuad, RadialTraitsType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  RadialTransformQuadrature(const RadialTraitsType& traits = RadialTraitsType()) :
    base_type(traits) { }
  RadialTransformQuadrature(const RadialTraits& traits) :
    RadialTransformQuadrature(radial_traits_cast<RadialTraitsType>(traits)) { }

  template <typename... Args,
    typename = std::enable_if_t<(sizeof...(Args) > 1)>
  >
  RadialTransformQuadrature(Args&&... args) :
    RadialTransformQuadrature(RadialTraitsType(std::forward<Args>(args)...)) { }
    

  RadialTransformQuadrature( const RadialTransformQuadrature& )     = default;
  RadialTransformQuadrature( RadialTransformQuadrature&& ) noexcept = default;

};

template <typename BaseQuad, typename RadialTraitsType>
struct quadrature_traits<RadialTransformQuadrature<BaseQuad, RadialTraitsType>> {

  using point_type       = typename BaseQuad::point_type;
  using weight_type      = typename BaseQuad::weight_type;
  using traits_type      = RadialTraitsType;
  using point_container  = typename BaseQuad::point_container;
  using weight_container = typename BaseQuad::weight_container;

  inline static std::tuple<point_container,weight_container>
    generate( const RadialTraitsType& traits ) {

    using base_quad_traits = quadrature_traits<BaseQuad>;

    const auto npts = traits.npts();
    point_container  points( npts );
    weight_container weights( npts );


    const auto npts_base = base_quad_traits::bound_inclusive ? npts+2 : npts;
    auto [base_x, base_w] = base_quad_traits::generate(npts_base);

    const auto ipts_offset = !!base_quad_traits::bound_inclusive;
    for(size_t i = 0; i < npts; ++i) {
      const auto xi = base_x[i + ipts_offset];
      const auto wi = base_w[i + ipts_offset];
      points[i]  = traits.radial_transform(xi);
      weights[i] = wi * traits.radial_jacobian(xi);
    }

    return std::make_tuple(points, weights);
  }

};
  


}
