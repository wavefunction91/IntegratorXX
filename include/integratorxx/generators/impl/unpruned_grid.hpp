#pragma once

#include <integratorxx/generators/spherical_factory.hpp>

#include <integratorxx/generators/impl/radial_types.hpp>
#include <integratorxx/generators/impl/s2_types.hpp>

namespace IntegratorXX {

namespace detail {

template <typename AngularQuadType>
auto generate_unpruned_grid_impl(RadialQuad rq, const RadialTraits& traits,
  AngularQuadType&& ang_quad) {

  switch( rq ) {

    case RadialQuad::Becke:
      return SphericalGridFactory::generate_unpruned_grid( bk_type(traits), std::forward<AngularQuadType>(ang_quad) );

    case RadialQuad::MuraKnowles:
      return SphericalGridFactory::generate_unpruned_grid( mk_type(traits), std::forward<AngularQuadType>(ang_quad) );

    case RadialQuad::MurrayHandyLaming:
      return SphericalGridFactory::generate_unpruned_grid( mhl_type(traits), std::forward<AngularQuadType>(ang_quad) );

    case RadialQuad::TreutlerAhlrichs:
      return SphericalGridFactory::generate_unpruned_grid( ta_type(traits), std::forward<AngularQuadType>(ang_quad) );

    case RadialQuad::LindhMalmqvistGagliardi:
      return SphericalGridFactory::generate_unpruned_grid( lmg_type(traits), std::forward<AngularQuadType>(ang_quad) );

    default:
      throw std::runtime_error("Unsupported Radial Quadrature");
      abort();

  }

}

} // Implementation details

SphericalGridFactory::spherical_grid_ptr 
  SphericalGridFactory::generate_unpruned_grid( RadialQuad rq, 
  const RadialTraits& traits, AngularQuad aq, AngularSize nang) {

  switch(aq) {
    case AngularQuad::AhrensBeylkin:
      return detail::generate_unpruned_grid_impl(rq, traits, ah_type(nang));
    case AngularQuad::Delley:
      return detail::generate_unpruned_grid_impl(rq, traits, de_type(nang));
    case AngularQuad::LebedevLaikov:
      return detail::generate_unpruned_grid_impl(rq, traits, ll_type(nang));
    case AngularQuad::Womersley:
      return detail::generate_unpruned_grid_impl(rq, traits, wo_type(nang));
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }
}

}
