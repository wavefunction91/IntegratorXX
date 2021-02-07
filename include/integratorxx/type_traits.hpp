#pragma once

#include "types.hpp"
#include <complex>

namespace IntegratorXX::detail {

template <typename T>
struct real {
  using type = T;
};

template <typename T>
struct real<std::complex<T>> {
  using type = T;
};

template <typename T>
using real_t = typename real<T>::type;






template <typename QuadratureType>
struct is_cartesian_quadrature {
  static constexpr bool value = 
    std::is_same< 
      typename QuadratureType::point_type, 
      cartesian_pt_t< real_t<typename QuadratureType::weight_type> >
    >::value;
};

template <typename QuadratureType>
inline constexpr bool is_cartesian_quadrature_v =
  is_cartesian_quadrature<QuadratureType>::value;

}
