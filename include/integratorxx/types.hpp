#pragma once

#include <array>
#include <type_traits>
       
#include "quadrature.hpp"
                                                

namespace IntegratorXX {

template < typename T >
using cartesian_pt_t = std::array<T, 3>;

//template <
//  typename T, 
//  template<typename...> class Container = std::vector, 
//  typename... Args
//>
//using CartesianQuadrature = 
//  QuadratureBase< 
//    Container< cartesian_pt_t<T>, Args... >, 
//    Container< T, Args...> 
//  >;

}
