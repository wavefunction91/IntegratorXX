#include <integratorxx/generators/spherical_factory.hpp>
#include <integratorxx/quadratures/radial.hpp>
#include <integratorxx/quadratures/s2.hpp>

#include <algorithm>

namespace IntegratorXX {


UnprunedSphericalGridSpecification::UnprunedSphericalGridSpecification(
  RadialQuad rq, const RadialTraits& traits, AngularQuad aq, AngularSize as) :
  radial_quad(rq), radial_traits(traits.clone()), angular_quad(aq), angular_size(as) { }


RadialQuad radial_from_string(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  if(name == "BECKE")             return RadialQuad::Becke;
  if(name == "MURAKNOWLES")       return RadialQuad::MuraKnowles;
  if(name == "MK")                return RadialQuad::MuraKnowles;
  if(name == "MURRAYHANDYLAMING") return RadialQuad::MurrayHandyLaming;
  if(name == "MHL")               return RadialQuad::MurrayHandyLaming;
  if(name == "TREUTLERAHLRICHS")  return RadialQuad::TreutlerAhlrichs;
  if(name == "TA")                return RadialQuad::TreutlerAhlrichs;

  throw std::runtime_error("Unrecognized Radial Quadrature");
}

AngularQuad angular_from_string(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  if(name == "AHRENSBEYLKIN") return AngularQuad::AhrensBeylkin;
  if(name == "AB")            return AngularQuad::AhrensBeylkin;
  if(name == "DELLEY")        return AngularQuad::Delley;
  if(name == "LEBEDEVLAIKOV") return AngularQuad::LebedevLaikov;
  if(name == "LEBEDEV")       return AngularQuad::LebedevLaikov;
  if(name == "LL")            return AngularQuad::LebedevLaikov;
  if(name == "WOMERSLEY")     return AngularQuad::Womersley;

  throw std::runtime_error("Unrecognized Angular Quadrature");
}


using spherical_grid_ptr = SphericalGridFactory::spherical_grid_ptr;

using bk_type  = Becke<double,double>;
using mk_type  = MuraKnowles<double,double>;
using mhl_type = MurrayHandyLaming<double,double>;
using ta_type  = TreutlerAhlrichs<double,double>;

using ah_type = AhrensBeylkin<double>;
using de_type = Delley<double>;
using ll_type = LebedevLaikov<double>;
using wo_type = Womersley<double>;

/************************/
/**** Unpruned Grids ****/
/************************/
template <typename AngularQuadType>
spherical_grid_ptr generate_unpruned_grid_impl(RadialQuad rq, const RadialTraits& traits,
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

    default:
      throw std::runtime_error("Unsupported Radial Quadrature");
      abort();

  }


}

spherical_grid_ptr SphericalGridFactory::generate_unpruned_grid( RadialQuad rq, 
  const RadialTraits& traits, AngularQuad aq, AngularSize nang) {

  switch(aq) {
    case AngularQuad::AhrensBeylkin:
      return generate_unpruned_grid_impl(rq, traits, ah_type(nang));
    case AngularQuad::Delley:
      return generate_unpruned_grid_impl(rq, traits, de_type(nang));
    case AngularQuad::LebedevLaikov:
      return generate_unpruned_grid_impl(rq, traits, ll_type(nang));
    case AngularQuad::Womersley:
      return generate_unpruned_grid_impl(rq, traits, wo_type(nang));
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }
}

spherical_grid_ptr SphericalGridFactory::generate_grid( 
  UnprunedSphericalGridSpecification gs ) {
  if(!gs.radial_traits) throw std::runtime_error("RadialTraits Not Set");
  return generate_unpruned_grid(gs.radial_quad, *gs.radial_traits, 
    gs.angular_quad, gs.angular_size);
}




/**********************/
/**** Pruned Grids ****/
/**********************/

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

spherical_grid_ptr SphericalGridFactory::generate_pruned_grid( 
  RadialQuad rq, const RadialTraits& traits, 
  const std::vector<PruningRegion>& pruning_regions) { 

  switch( rq ) {

    case RadialQuad::Becke:
      return make_pruned_grid( bk_type(traits), pruning_regions );
    case RadialQuad::MuraKnowles:
      return make_pruned_grid( mk_type(traits), pruning_regions );
    case RadialQuad::MurrayHandyLaming:
      return make_pruned_grid( mhl_type(traits), pruning_regions );
    case RadialQuad::TreutlerAhlrichs:
      return make_pruned_grid( ta_type(traits), pruning_regions );

    default:
      throw std::runtime_error("Unsupported Radial Quadrature");
      abort();

  }

}

spherical_grid_ptr SphericalGridFactory::generate_grid( 
  PrunedSphericalGridSpecification gs ) {
  if(!gs.radial_traits) throw std::runtime_error("RadialTraits Not Set");
  return generate_pruned_grid( gs.radial_quad, *gs.radial_traits, 
    gs.pruning_regions );
}




/***************************
 * Default Pruning Schemes *
 ***************************/


/*** Psi4 Robust ***/


template <typename AngularQuad>
auto get_robust_low_med_sizes(AngularSize asz) {
  using traits = quadrature_traits<AngularQuad>;
  const auto base_order = traits::algebraic_order_by_npts(asz);
  if( base_order < 0 ) throw std::runtime_error("Invalid Base Grid");

  const auto med_order = 
    traits::next_algebraic_order(base_order > 6 ? base_order-6 : base_order);
  const auto low_order = traits::next_algebraic_order(7);

  AngularSize med_sz(traits::npts_by_algebraic_order(med_order));
  AngularSize low_sz(traits::npts_by_algebraic_order(low_order));

  return std::make_pair(low_sz, med_sz);
}


PrunedSphericalGridSpecification robust_psi4_pruning_scheme_impl(
  size_t low_sz, size_t med_sz, AngularQuad angular_quad, 
  UnprunedSphericalGridSpecification unp ) {

  const size_t rsz = unp.radial_traits->npts();
  const size_t r_div_4 = rsz / 4ul + 1ul;
  const size_t r_div_2 = rsz / 2ul + 1ul;
  std::vector<PruningRegion> pruning_regions = {
    {0ul,     r_div_4, angular_quad, low_sz},
    {r_div_4, r_div_2, angular_quad, med_sz},
    {r_div_2,     rsz, angular_quad, unp.angular_size}
  };

  return PrunedSphericalGridSpecification(
    unp.radial_quad, unp.radial_traits->clone(), pruning_regions
  );
  
}

PrunedSphericalGridSpecification robust_psi4_pruning_scheme(
  UnprunedSphericalGridSpecification unp ) {

  size_t low_sz, med_sz;
  const auto angular_quad = unp.angular_quad;
  const auto asz = unp.angular_size;
  switch(angular_quad) {
    case AngularQuad::AhrensBeylkin:
      std::tie(low_sz, med_sz) = get_robust_low_med_sizes<ah_type>(asz);
      break;
    case AngularQuad::Delley:
      std::tie(low_sz, med_sz) = get_robust_low_med_sizes<de_type>(asz);
      break;
    case AngularQuad::LebedevLaikov:
      std::tie(low_sz, med_sz) = get_robust_low_med_sizes<ll_type>(asz);
      break;
    case AngularQuad::Womersley:
      std::tie(low_sz, med_sz) = get_robust_low_med_sizes<wo_type>(asz);
      break;
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }

  return robust_psi4_pruning_scheme_impl(low_sz, med_sz, angular_quad, unp);
}



/*** Treutler-Ahlrichs Pruning ***/


template <typename AngularQuad>
auto get_treutler_low_med_sizes() {
  using traits = quadrature_traits<AngularQuad>;
  AngularSize med_sz(traits::npts_by_algebraic_order(traits::next_algebraic_order(11)));
  AngularSize low_sz(traits::npts_by_algebraic_order(traits::next_algebraic_order(7 )));

  return std::make_pair(low_sz, med_sz);
}

PrunedSphericalGridSpecification treutler_pruning_scheme_impl(
  size_t low_sz, size_t med_sz, AngularQuad angular_quad, 
  UnprunedSphericalGridSpecification unp ) {

  const size_t rsz = unp.radial_traits->npts();
  const size_t r_div_3 = rsz / 3ul + 1ul;
  const size_t r_div_2 = rsz / 2ul + 1ul;
  std::vector<PruningRegion> pruning_regions = {
    {0ul,     r_div_3, angular_quad, low_sz},
    {r_div_3, r_div_2, angular_quad, med_sz},
    {r_div_2, rsz,     angular_quad, unp.angular_size}
  };

  return PrunedSphericalGridSpecification(
    unp.radial_quad, unp.radial_traits->clone(), pruning_regions
  );
  
}



PrunedSphericalGridSpecification treutler_pruning_scheme(
  UnprunedSphericalGridSpecification unp ) {

  size_t low_sz, med_sz;
  const auto angular_quad = unp.angular_quad;
  switch(angular_quad) {
    case AngularQuad::AhrensBeylkin:
      std::tie(low_sz, med_sz) = get_treutler_low_med_sizes<ah_type>();
      break;
    case AngularQuad::Delley:
      std::tie(low_sz, med_sz) = get_treutler_low_med_sizes<de_type>();
      break;
    case AngularQuad::LebedevLaikov:
      std::tie(low_sz, med_sz) = get_treutler_low_med_sizes<ll_type>();
      break;
    case AngularQuad::Womersley:
      std::tie(low_sz, med_sz) = get_treutler_low_med_sizes<wo_type>();
      break;
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }

  return treutler_pruning_scheme_impl(low_sz, med_sz, angular_quad, unp);
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
