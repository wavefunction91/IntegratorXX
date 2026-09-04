#pragma once

#include <array>
#include <integratorxx/config.hpp>
#include <type_traits>

namespace IntegratorXX {

template <typename T>
using cartesian_pt_t = std::array<T, 3>;

}
