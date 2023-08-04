#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 2 grid with 6 points
 * 
 */
template <typename T>
struct womersley_6 {

  static constexpr std::array<cartesian_pt_t<T>,6> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     1.0000000000000000e+00,      0.0000000000000000e+00,      0.0000000000000000e+00,
     0.0000000000000000e+00,     -1.0000000000000000e+00,      0.0000000000000000e+00,
     0.0000000000000000e+00,      1.0000000000000000e+00,      0.0000000000000000e+00,
     0.0000000000000000e+00,      0.0000000000000000e+00,     -1.0000000000000000e+00,
    -1.0000000000000000e+00,      0.0000000000000000e+00,      0.0000000000000000e+00
};


static constexpr auto weights = detail::create_array<6, T>(4.0 * M_PI / 6.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
