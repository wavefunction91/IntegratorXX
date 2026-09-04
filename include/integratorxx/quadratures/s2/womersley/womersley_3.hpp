#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 1 grid with 3 points
 * 
 */
template <typename T>
struct womersley_3 {

  static constexpr std::array<cartesian_pt_t<ixx_real>,3> points = {
     IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(1.0000000000000000e+00),
     IXX_REAL(8.6602540378443871e-01),      IXX_REAL(0.0000000000000000e+00),     IXX_REAL(-4.9999999999999978e-01),
    IXX_REAL(-8.6602540378443871e-01),     IXX_REAL(-2.7853501340422215e-16),     IXX_REAL(-4.9999999999999978e-01)
};


  /// Equal-weight grid: every weight is 4*pi/3 (sphere area / npts).
  /// copy_grid forms it exactly; see util/copy_grid.hpp.
  static constexpr ixx_int uniform_weight_npts = 3;
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
