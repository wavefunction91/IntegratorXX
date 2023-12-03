#include <integratorxx/generators/s2_factory.hpp>

#include "s2_types.hpp"

#include <algorithm>

namespace IntegratorXX {

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

using s2_grid_ptr = S2Factory::s2_grid_ptr;

s2_grid_ptr S2Factory::generate(AngularQuad aq, size_t npts) {

  switch(aq) {
    case AngularQuad::AhrensBeylkin:
      return std::make_unique<ah_type>(npts);
    case AngularQuad::Delley:
      return std::make_unique<de_type>(npts);
    case AngularQuad::LebedevLaikov:
      return std::make_unique<ll_type>(npts);
    case AngularQuad::Womersley:
      return std::make_unique<wo_type>(npts);
    default:
      throw std::runtime_error("Unsupported Angular Quadrature");
      abort();
  } 
}

}
