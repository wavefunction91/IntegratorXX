#pragma once

namespace IntegratorXX {
namespace LebedevLaikovGrids {

/**
 *  \brief Lebedev-Laikov Quadrature specification for Order = 6
 *
 */
template <typename T>
struct lebedev_laikov_6 {
  static constexpr std::array<cartesian_pt_t<T>, 6> points = {
      1.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
      -1.000000000000000e+00, 0.000000000000000e+00,  0.000000000000000e+00,
      0.000000000000000e+00,  1.000000000000000e+00,  0.000000000000000e+00,
      0.000000000000000e+00,  -1.000000000000000e+00, 0.000000000000000e+00,
      0.000000000000000e+00,  0.000000000000000e+00,  1.000000000000000e+00,
      0.000000000000000e+00,  0.000000000000000e+00,  -1.000000000000000e+00};

  static constexpr std::array<T, 6> weights = {
      1.666666666666667e-01, 1.666666666666667e-01, 1.666666666666667e-01,
      1.666666666666667e-01, 1.666666666666667e-01, 1.666666666666667e-01};
};
}  // namespace LebedevLaikovGrids
}  // namespace IntegratorXX
