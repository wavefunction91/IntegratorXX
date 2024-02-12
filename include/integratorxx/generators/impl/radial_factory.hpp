#pragma once
#include <integratorxx/generators/radial_factory.hpp>

#include <integratorxx/generators/impl/radial_types.hpp>

#include <algorithm>

namespace IntegratorXX {

RadialQuad radial_from_string(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  if(name == "BECKE")             return RadialQuad::Becke;
  if(name == "MURAKNOWLES")       return RadialQuad::MuraKnowles;
  if(name == "MK")                return RadialQuad::MuraKnowles;
  if(name == "MURRAYHANDYLAMING") return RadialQuad::MurrayHandyLaming;
  if(name == "MHL")               return RadialQuad::MurrayHandyLaming;
  if(name == "TREUTLERAHLRICHS")  return RadialQuad::TreutlerAhlrichs;
  if(name == "TA")                return RadialQuad::TreutlerAhlrichs;
  if(name == "LINDHMALMQVISTGAGLIARDI") return RadialQuad::LindhMalmqvistGagliardi;
  if(name == "LMG")               return RadialQuad::LindhMalmqvistGagliardi;

  throw std::runtime_error("Unrecognized Radial Quadrature");
}

RadialFactory::radial_grid_ptr RadialFactory::generate(RadialQuad rq, const RadialTraits& traits) {

  switch(rq) {
    case RadialQuad::Becke:
      return std::make_unique<bk_type>(traits);
    case RadialQuad::MuraKnowles:
      return std::make_unique<mk_type>(traits);
    case RadialQuad::MurrayHandyLaming:
      return std::make_unique<mhl_type>(traits);
    case RadialQuad::TreutlerAhlrichs:
      return std::make_unique<ta_type>(traits);
    case RadialQuad::LindhMalmqvistGagliardi:
      return std::make_unique<lmg_type>(traits);
    default:
      throw std::runtime_error("Unsupported Radial Quadrature");
      abort();
  } 
}


}
