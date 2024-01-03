#pragma once
#include <integratorxx/quadratures/radial.hpp>

namespace IntegratorXX {

/// High-level specification of radial quadratures
enum class RadialQuad : uint32_t {
  Becke                   = 0x0010, 
  MurrayHandyLaming       = 0x0020,
  MuraKnowles             = 0x0030, 
  TreutlerAhlrichs        = 0x0040,
  LindhMalmqvistGagliardi = 0x0050
};

template <typename RadQuadType>
RadialQuad radial_from_type() {
  if constexpr (detail::is_becke_v<RadQuadType>) return RadialQuad::Becke;
  if constexpr (detail::is_mk_v<RadQuadType>   ) return RadialQuad::MuraKnowles;
  if constexpr (detail::is_mhl_v<RadQuadType>)   return RadialQuad::MurrayHandyLaming;
  if constexpr (detail::is_ta_v<RadQuadType>)    return RadialQuad::TreutlerAhlrichs;
  if constexpr (detail::is_lmg_v<RadQuadType>)   return RadialQuad::LindhMalmqvistGagliardi;

  throw std::runtime_error("Unrecognized Radial Quadrature");
};

RadialQuad radial_from_string(std::string name);

namespace detail {

template <typename RadialTraitsType, typename... Args>
std::unique_ptr<RadialTraits> make_radial_traits(Args&&... args) {
  using traits_type = RadialTraitsType;
  if constexpr (std::is_constructible_v<traits_type,Args...>)
    return std::make_unique<traits_type>(std::forward<Args>(args)...);
  else return nullptr;
}

}

template <typename... Args>
std::unique_ptr<RadialTraits> make_radial_traits(RadialQuad rq, Args&&... args) {
  std::unique_ptr<RadialTraits> ptr;
  switch(rq) {
    case RadialQuad::Becke:
      ptr = 
        detail::make_radial_traits<BeckeRadialTraits>(std::forward<Args>(args)...);
      break;
    case RadialQuad::MurrayHandyLaming:
      ptr = 
        detail::make_radial_traits<MurrayHandyLamingRadialTraits<2>>(std::forward<Args>(args)...);
      break;
    case RadialQuad::MuraKnowles:
      ptr = 
        detail::make_radial_traits<MuraKnowlesRadialTraits>(std::forward<Args>(args)...);
      break;
    case RadialQuad::TreutlerAhlrichs:
      ptr = 
        detail::make_radial_traits<TreutlerAhlrichsRadialTraits>(std::forward<Args>(args)...);
      break;
    case RadialQuad::LindhMalmqvistGagliardi:
      ptr = 
        detail::make_radial_traits<LindhMalmqvistGagliardiRadialTraits>(std::forward<Args>(args)...);
      break;
  }

  if(!ptr) throw std::runtime_error("RadialTraits Construction Failed");
  return ptr;
}


struct RadialFactory {

  using radial_grid_ptr = std::shared_ptr<
    QuadratureBase<
      std::vector<double>,
      std::vector<double>
    >
  >;

  static radial_grid_ptr generate(RadialQuad rq, const RadialTraits& traits);
  
};

}
