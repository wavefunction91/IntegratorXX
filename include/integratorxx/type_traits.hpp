#pragma once

#include <complex>

#include "types.hpp"

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
  static constexpr bool value = std::is_same<
      typename QuadratureType::point_type,
      cartesian_pt_t<real_t<typename QuadratureType::weight_type>>>::value;
};

template <typename QuadratureType>
inline constexpr bool is_cartesian_quadrature_v =
    is_cartesian_quadrature<QuadratureType>::value;

template <typename T, typename... Args>
using all_are_not = std::conjunction<std::negation<std::is_same<T, Args>>...>;

template <typename T, typename... Args>
using enable_if_all_are_not = std::enable_if<all_are_not<T, Args...>::value>;

template <typename T, typename... Args>
using enable_if_all_are_not_t =
    typename enable_if_all_are_not<T, Args...>::type;

}  // namespace IntegratorXX::detail
