#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 1 grid with 3 points
 * 
 */
template <typename T>
struct womersley_3 {

  static constexpr std::array<cartesian_pt_t<T>,3> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     8.6602540378443871e-01,      0.0000000000000000e+00,     -4.9999999999999978e-01,
    -8.6602540378443871e-01,     -2.7853501340422215e-16,     -4.9999999999999978e-01
};


static constexpr auto weights = detail::create_array<3, T>(4.0 * M_PI / 3.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
