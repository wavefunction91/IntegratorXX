#pragma once
#include <memory>
#include <vector>
#include <array>
#include <integratorxx/quadratures/radial.hpp>
#include <integratorxx/quadratures/s2.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>

namespace IntegratorXX {

/// High-level specification of radial quadratures
enum class RadialQuad : uint32_t {
  Becke             = 0x0010, 
  MurrayHandyLaming = 0x0020,
  MuraKnowles       = 0x0030, 
  TreutlerAhlrichs  = 0x0040
};

template <typename RadQuadType>
RadialQuad rad_from_type() {
  if constexpr (detail::is_becke_v<RadQuadType>) return RadialQuad::Becke;
  if constexpr (detail::is_mk_v<RadQuadType>   ) return RadialQuad::MuraKnowles;
  if constexpr (detail::is_mhl_v<RadQuadType>)   return RadialQuad::MurrayHandyLaming;
  if constexpr (detail::is_ta_v<RadQuadType>)    return RadialQuad::TreutlerAhlrichs;

  throw std::runtime_error("Unrecognized Radial Quadrature");
};

/// High-level specification of angular quadratures
enum class AngularQuad : uint32_t {
  AhrensBeylkin     = 0x0100,
  Delley            = 0x0200,
  LebedevLaikov     = 0x0300,
  Womersley         = 0x0400
};

template <typename AngQuadType>
AngularQuad ang_from_type() {
  if constexpr (detail::is_ahrens_beyklin_v<AngQuadType>) return AngularQuad::AhrensBeylkin;
  if constexpr (detail::is_delley_v<AngQuadType>   ) return AngularQuad::Delley;
  if constexpr (detail::is_lebedev_laikov_v<AngQuadType>)   return AngularQuad::LebedevLaikov;
  if constexpr (detail::is_womersley_v<AngQuadType>)    return AngularQuad::Womersley;

  throw std::runtime_error("Unrecognized Angular Quadrature");
};

/// High-level specification of pruning schemes for spherical quadratures
enum class PruningScheme {
  Unpruned, /// Unpruned quadrature
  Robust,   /// The "Robust" scheme of Psi4
  Treutler  /// The Treutler-Aldrichs scheme
};

// TODO: Make these strong (non-convertible) types
using RadialScale = double;
using RadialSize  = size_t;
using AngularSize = size_t;


/// Generic specification of an unpruned spherical quadrature
struct UnprunedSphericalGridSpecification {
  RadialQuad  radial_quad;  ///< Radial quadrature specification
  RadialSize  radial_size;  ///< Number of radial quadrature points
  RadialScale radial_scale; ///< Radial scaling factor

  AngularQuad angular_quad; /// Angular quadrature specification
  AngularSize angular_size; /// Number of angular quadrature points
};


/// Speficiation of a pruned region of an spierical quadrature
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
  RadialSize  radial_size;  ///< Number of radial quadrature points
  RadialScale radial_scale; ///< Radial scaling factor

  std::vector<PruningRegion> pruning_regions; ///< List of pruning regions over the radial quadrature

  inline bool operator==(const PrunedSphericalGridSpecification& other) const noexcept {
    return radial_quad == other.radial_quad and
           radial_size == other.radial_size and
           radial_scale == other.radial_scale and 
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
   *  @brief Generate an unpruned spherical grid given a suppled radial and
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


  static spherical_grid_ptr generate_unpruned_grid( RadialQuad, RadialSize, 
    RadialScale, AngularQuad, AngularSize );
  static spherical_grid_ptr generate_pruned_grid( RadialQuad, RadialSize, 
    RadialScale, const std::vector<PruningRegion>&);


  static spherical_grid_ptr generate_grid(UnprunedSphericalGridSpecification gs);
  static spherical_grid_ptr generate_grid(PrunedSphericalGridSpecification gs); 

};

}
