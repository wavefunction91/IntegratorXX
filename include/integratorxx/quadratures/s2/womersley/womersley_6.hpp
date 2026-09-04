#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 2 grid with 6 points
 * 
 */
template <typename T>
struct womersley_6 {

  static constexpr std::array<cartesian_pt_t<ixx_real>,6> points = {
     IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(1.0000000000000000e+00),
     IXX_REAL(1.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),
     IXX_REAL(0.0000000000000000e+00),     IXX_REAL(-1.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),
     IXX_REAL(0.0000000000000000e+00),      IXX_REAL(1.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),
     IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),     IXX_REAL(-1.0000000000000000e+00),
    IXX_REAL(-1.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00)
};


  /// Equal-weight grid: every weight is 4*pi/6 (sphere area / npts).
  /// copy_grid forms it exactly; see util/copy_grid.hpp.
  static constexpr ixx_int uniform_weight_npts = 6;
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
