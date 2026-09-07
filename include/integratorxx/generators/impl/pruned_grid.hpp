#pragma once

#include <integratorxx/generators/spherical_factory.hpp>

#include <integratorxx/generators/impl/radial_types.hpp>
#include <integratorxx/generators/impl/s2_types.hpp>

#include <integratorxx/generators/impl/robust_pruning.hpp>
#include <integratorxx/generators/impl/treutler_pruning.hpp>

namespace IntegratorXX {

namespace detail {

template <typename AngularQuadType, typename RadialQuadType>
auto make_pruned_grid_impl(const RadialQuadType& rq, 
  const std::vector<PruningRegion>& pruning_regions) {

  RadialGridPartition<AngularQuadType> rgp;
  for( auto& region : pruning_regions ) {
    rgp.add_quad( rq, region.idx_st, AngularQuadType(region.angular_size) );
  }
  rgp.finalize(rq);

  return SphericalGridFactory::generate_pruned_grid(rq, std::move(rgp));

}

template <typename RadialQuadType>
auto make_pruned_grid(const RadialQuadType& rq, 
  const std::vector<PruningRegion>& pruning_regions) {

  if(pruning_regions.size() == 0)
    throw std::runtime_error("No Pruning Regions");

  auto angular_quad = pruning_regions[0].angular_quad;
  for(auto r : pruning_regions) {
    if(r.angular_quad != angular_quad)
      throw std::runtime_error("Mixed Angular Pruning Not Supported");
  }

  switch(angular_quad) {
    case AngularQuad::AhrensBeylkin:
      return make_pruned_grid_impl<ah_type>(rq, pruning_regions);
    case AngularQuad::Delley:
      return make_pruned_grid_impl<de_type>(rq, pruning_regions);
    case AngularQuad::LebedevLaikov:
      return make_pruned_grid_impl<ll_type>(rq, pruning_regions);
    case AngularQuad::Womersley:
      return make_pruned_grid_impl<wo_type>(rq, pruning_regions);
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }
  

}

} // Implementation Details

SphericalGridFactory::spherical_grid_ptr 
  SphericalGridFactory::generate_pruned_grid( RadialQuad rq, 
  const RadialTraits& traits, 
  const std::vector<PruningRegion>& pruning_regions) { 

  switch( rq ) {

    case RadialQuad::Becke:
      return detail::make_pruned_grid( bk_type(traits), pruning_regions );
    case RadialQuad::MuraKnowles:
      return detail::make_pruned_grid( mk_type(traits), pruning_regions );
    case RadialQuad::MurrayHandyLaming:
      return detail::make_pruned_grid( mhl_type(traits), pruning_regions );
    case RadialQuad::TreutlerAhlrichs:
      return detail::make_pruned_grid( ta_type(traits), pruning_regions );
    case RadialQuad::LindhMalmqvistGagliardi:
      return detail::make_pruned_grid( lmg_type(traits), pruning_regions );

    default:
      throw std::runtime_error("Unsupported Radial Quadrature");
      abort();

  }

}


PrunedSphericalGridSpecification create_pruned_spec(
  PruningScheme scheme, UnprunedSphericalGridSpecification unp
) {

  if(!unp.radial_traits) throw std::runtime_error("RadialTraits Not Set");
  switch(scheme) {
    case PruningScheme::Robust:
      return robust_psi4_pruning_scheme(unp);
    case PruningScheme::Treutler:
      return treutler_pruning_scheme(unp);
    
    // Default to Unpruned Grid
    case PruningScheme::Unpruned:
    default:
      std::vector<PruningRegion> pruning_regions = {
        {0ul, unp.radial_traits->npts(), unp.angular_quad, unp.angular_size}
      };
      return PrunedSphericalGridSpecification(
        unp.radial_quad, unp.radial_traits->clone(), pruning_regions
      );
  }

}

}
