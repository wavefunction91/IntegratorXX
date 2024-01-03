#pragma once
#include <integratorxx/generators/spherical_factory.hpp>

#include <integratorxx/generators/impl/s2_types.hpp>

namespace IntegratorXX {

namespace detail {

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

} // Implementation Details


PrunedSphericalGridSpecification robust_psi4_pruning_scheme(
  UnprunedSphericalGridSpecification unp ) {

  size_t low_sz, med_sz;
  const auto angular_quad = unp.angular_quad;
  const auto asz = unp.angular_size;
  switch(angular_quad) {
    case AngularQuad::AhrensBeylkin:
      std::tie(low_sz, med_sz) = detail::get_robust_low_med_sizes<ah_type>(asz);
      break;
    case AngularQuad::Delley:
      std::tie(low_sz, med_sz) = detail::get_robust_low_med_sizes<de_type>(asz);
      break;
    case AngularQuad::LebedevLaikov:
      std::tie(low_sz, med_sz) = detail::get_robust_low_med_sizes<ll_type>(asz);
      break;
    case AngularQuad::Womersley:
      std::tie(low_sz, med_sz) = detail::get_robust_low_med_sizes<wo_type>(asz);
      break;
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  }

  return detail::robust_psi4_pruning_scheme_impl(low_sz, med_sz, angular_quad, unp);
}

}
