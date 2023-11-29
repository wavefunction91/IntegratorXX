#pragma once
#include <integratorxx/quadratures/radial.hpp>

namespace IntegratorXX {

/// High-level specification of radial quadratures
enum class RadialQuad : uint32_t {
  Becke             = 0x0010, 
  MurrayHandyLaming = 0x0020,
  MuraKnowles       = 0x0030, 
  TreutlerAhlrichs  = 0x0040
};

template <typename RadQuadType>
RadialQuad radial_from_type() {
  if constexpr (detail::is_becke_v<RadQuadType>) return RadialQuad::Becke;
  if constexpr (detail::is_mk_v<RadQuadType>   ) return RadialQuad::MuraKnowles;
  if constexpr (detail::is_mhl_v<RadQuadType>)   return RadialQuad::MurrayHandyLaming;
  if constexpr (detail::is_ta_v<RadQuadType>)    return RadialQuad::TreutlerAhlrichs;

  throw std::runtime_error("Unrecognized Radial Quadrature");
};

RadialQuad radial_from_string(std::string name);

#if 0
struct RadialFactory {

  using radial_grid_ptr = std::shared_ptr<
    QuadratureBase<
      std::vector<double>,
      std::vector<double>
    >;

  //radial_grid_ptr generate(RadialTraits traits);
  
};
#endif

}
