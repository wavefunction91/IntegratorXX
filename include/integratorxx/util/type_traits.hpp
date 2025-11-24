#pragma once
#include <type_traits>

namespace IntegratorXX::detail {

template <typename T, typename U = void>
using enable_if_integral = std::enable_if<std::is_integral_v<T>,U>;

template <typename T, typename U = void>
using enable_if_integral_t = typename enable_if_integral<T,U>::type;

template <typename T, typename U, typename V = void>
using enable_if_float_and_integral = 
  std::enable_if<std::is_floating_point_v<T> and std::is_integral_v<U>, V>;

template <typename T, typename U, typename V = void>
using enable_if_float_and_integral_t = 
  typename enable_if_float_and_integral<T,U,V>::type; 

}
