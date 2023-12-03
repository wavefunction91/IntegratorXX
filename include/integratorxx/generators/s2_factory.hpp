#pragma once
#include <integratorxx/quadratures/s2.hpp>

#include <memory>

namespace IntegratorXX {

/// High-level specification of angular quadratures
enum class AngularQuad : uint32_t {
  AhrensBeylkin     = 0x0100,
  Delley            = 0x0200,
  LebedevLaikov     = 0x0300,
  Womersley         = 0x0400
};

template <typename AngQuadType>
AngularQuad angular_from_type() {
  if constexpr (detail::is_ahrens_beyklin_v<AngQuadType>) return AngularQuad::AhrensBeylkin;
  if constexpr (detail::is_delley_v<AngQuadType>   ) return AngularQuad::Delley;
  if constexpr (detail::is_lebedev_laikov_v<AngQuadType>)   return AngularQuad::LebedevLaikov;
  if constexpr (detail::is_womersley_v<AngQuadType>)    return AngularQuad::Womersley;

  throw std::runtime_error("Unrecognized Angular Quadrature");
};

AngularQuad angular_from_string(std::string name);

struct S2Factory {

  using s2_grid_ptr = std::shared_ptr<
    QuadratureBase<
      std::vector<std::array<double,3>>,
      std::vector<double>
    >
  >;

  static s2_grid_ptr generate(AngularQuad aq, size_t npts);
  
};

}
