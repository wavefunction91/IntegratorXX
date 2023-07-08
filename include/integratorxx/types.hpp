#pragma once

#include <array>
#include <type_traits>

#include "quadrature.hpp"

namespace IntegratorXX {

template <typename T>
using cartesian_pt_t = std::array<T, 3>;

}  // namespace IntegratorXX
