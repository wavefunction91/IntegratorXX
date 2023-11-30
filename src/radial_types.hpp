#pragma once
#include <integratorxx/quadratures/radial.hpp>

namespace IntegratorXX {

using bk_type  = Becke<double,double>;
using mk_type  = MuraKnowles<double,double>;
using mhl_type = MurrayHandyLaming<double,double>;
using ta_type  = TreutlerAhlrichs<double,double>;

}

