#pragma once
#include <memory>
#include <vector>
#include <array>
#include <integratorxx/generators/radial_factory.hpp>
#include <integratorxx/generators/s2_factory.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>

namespace IntegratorXX {

/// High-level specification of pruning schemes for spherical quadratures
enum class PruningScheme {
  Unpruned, /// Unpruned quadrature
  Robust,   /// The "Robust" scheme of Psi4
  Treutler  /// The Treutler-Ahlrichs scheme
};

// TODO: Make these strong (non-convertible) types
//using RadialScale = double;
//using RadialSize  = size_t;
using AngularSize = size_t;
using radial_traits_ptr = std::unique_ptr<RadialTraits>;

/// Generic specification of an unpruned spherical quadrature
struct UnprunedSphericalGridSpecification {
  RadialQuad   radial_quad;   ///< Radial quadrature specification
  radial_traits_ptr radial_traits; ///< Radial traits (order, scaling factors, etc)

  AngularQuad angular_quad; /// Angular quadrature specification
  AngularSize angular_size; /// Number of angular quadrature points

  UnprunedSphericalGridSpecification( RadialQuad rq, const RadialTraits& traits, 
    AngularQuad aq, AngularSize as) : radial_quad(rq), radial_traits(traits.clone()), 
    angular_quad(aq), angular_size(as) {}

  UnprunedSphericalGridSpecification(const UnprunedSphericalGridSpecification& other) :
    radial_quad(other.radial_quad), radial_traits(other.radial_traits ? other.radial_traits->clone() : nullptr),
    angular_quad(other.angular_quad), angular_size(other.angular_size) {}
};



/// Specification of a pruned region of an spherical quadrature
struct PruningRegion {
  size_t idx_st; ///< Starting radial index for pruned region
  size_t idx_en; ///< Ending radial index (exclusive) for the pruned region

  AngularQuad angular_quad;  ///< Angular quadrature for the pruned region
  AngularSize angular_size;  ///< Number of angular quadrature points in the pruned region

  /// Check equality of `PruningRegion` instances
  inline bool operator==(const PruningRegion& other) const noexcept {
    return other.idx_st == idx_st and 
           other.idx_en == idx_en and
           other.angular_size == angular_size;
  }
};

struct PrunedSphericalGridSpecification {
  RadialQuad  radial_quad;  ///< Radial quadrature specification
  radial_traits_ptr radial_traits; ///< Radial traits (order, scaling factors, etc)

  std::vector<PruningRegion> pruning_regions; ///< List of pruning regions over the radial quadrature
  
  PrunedSphericalGridSpecification() = default;

  template <typename... Arg>
  PrunedSphericalGridSpecification(RadialQuad rq, radial_traits_ptr&& traits, Arg&&... arg) :
    radial_quad(rq), radial_traits(std::move(traits)), pruning_regions(std::forward<Arg>(arg)...) { }
  template <typename... Arg>
  PrunedSphericalGridSpecification(RadialQuad rq, const RadialTraits& traits, Arg&&... arg) :
    PrunedSphericalGridSpecification(rq, traits.clone(), std::forward<Arg>(arg)...) { }
  
  PrunedSphericalGridSpecification(const PrunedSphericalGridSpecification& other) :
    PrunedSphericalGridSpecification(other.radial_quad,
      other.radial_traits ? other.radial_traits->clone() : nullptr,
      other.pruning_regions) { }

  PrunedSphericalGridSpecification& operator=(const PrunedSphericalGridSpecification& other) {
    radial_quad = other.radial_quad;
    radial_traits = other.radial_traits ? other.radial_traits->clone() : nullptr;
    pruning_regions = other.pruning_regions;
    return *this;
  }

  inline bool operator==(const PrunedSphericalGridSpecification& other) const noexcept {
    return radial_quad == other.radial_quad and
           (radial_traits ? (other.radial_traits and radial_traits->compare(*other.radial_traits)) : !other.radial_traits) and 
           pruning_regions == other.pruning_regions;
  }
};

/// Generate a "Robust"-Psi4 Pruning specification from an 
/// unpruned quadrature specification
PrunedSphericalGridSpecification robust_psi4_pruning_scheme(
  UnprunedSphericalGridSpecification
);

/// Generate a Pruning specification according to the Treutler-Ahlrichs 
/// scheme from an unpruned specification
PrunedSphericalGridSpecification treutler_pruning_scheme(
  UnprunedSphericalGridSpecification
);

/// Generate a pruning specification from a specificed pruning scheme and 
/// an unpruned grid specification
PrunedSphericalGridSpecification create_pruned_spec(
  PruningScheme, UnprunedSphericalGridSpecification
);



/// Factory for spherical grids
struct SphericalGridFactory {

  template <typename RadialType, typename AngularType>
  using unpruned_sphere_type = SphericalQuadrature< 
      std::decay_t<RadialType>, std::decay_t<AngularType>
    >;
  template <typename RadialType, typename AngularType>
  using pruned_sphere_type = PrunedSphericalQuadrature< 
      std::decay_t<RadialType>, std::decay_t<AngularType>
    >;

  using spherical_grid_ptr = std::shared_ptr<
    SphericalQuadratureBase<
      std::vector<std::array<double,3>>,
      std::vector<double>
    >
  >;

  /**
   *  @brief Generate an unpruned spherical grid given a supplied radial and
   *         angular quadrature.
   *
   *  All arguments are passed with perfect forwarding
   *
   *  @tparam RadialType Type of the radial quadrature
   *  @tparam AngularType Type of the angular quadrature
   * 
   *  @oaram     rq Radial quadrature from which to construct the spherical quadrature. 
   *  @oaram     aq Angular quadrature from which to construct the spherical quadrature. 
   *  @param[in] bsz Batch size for the grid generation.
   */
  template <typename RadialType, typename AngularType>
  static spherical_grid_ptr generate_unpruned_grid( RadialType&& rq, 
    AngularType&& aq ) {
    using sphere_type = unpruned_sphere_type<RadialType,AngularType>;
    return std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<AngularType>(aq) 
    );
  }

  template <typename RadialType, typename RadialPartitionType>
  static spherical_grid_ptr generate_pruned_grid( RadialType&& rq, 
    RadialPartitionType&& rgp ) {
    using angular_type = typename std::decay_t<RadialPartitionType>::angular_type;
    using sphere_type = pruned_sphere_type<RadialType,angular_type>;
    return std::make_shared<sphere_type>( 
      std::forward<RadialType>(rq), std::forward<RadialPartitionType>(rgp)
    );
  }


  static spherical_grid_ptr generate_unpruned_grid( RadialQuad, const RadialTraits&, 
    AngularQuad, AngularSize );
  static spherical_grid_ptr generate_pruned_grid( RadialQuad, const RadialTraits&,
    const std::vector<PruningRegion>&);


  static inline spherical_grid_ptr generate_grid(
    UnprunedSphericalGridSpecification gs) {
    if(!gs.radial_traits) throw std::runtime_error("RadialTraits Not Set");
    return generate_unpruned_grid(gs.radial_quad, *gs.radial_traits, 
      gs.angular_quad, gs.angular_size);
  }


  static inline spherical_grid_ptr generate_grid(
    PrunedSphericalGridSpecification gs) {
    if(!gs.radial_traits) throw std::runtime_error("RadialTraits Not Set");
    return generate_pruned_grid( gs.radial_quad, *gs.radial_traits, 
      gs.pruning_regions );
  } 

};

}
