#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 3 grid with 8 points
 * 
 */
template <typename T>
struct womersley_8 {

  static constexpr std::array<cartesian_pt_t<ixx_real>,8> points = {
     IXX_REAL(0.0000000000000000e+00),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(1.0000000000000000e+00),
     IXX_REAL(9.8198452026910221e-01),      IXX_REAL(0.0000000000000000e+00),      IXX_REAL(1.8896137687861275e-01),
     IXX_REAL(1.8319459063686479e-02),     IXX_REAL(-9.4263104463438019e-01),      IXX_REAL(3.3333333333333326e-01),
    IXX_REAL(-8.2550216058107961e-01),      IXX_REAL(4.5545040538444831e-01),      IXX_REAL(3.3333333333333326e-01),
     IXX_REAL(4.8331619830270472e-01),     IXX_REAL(-4.5545040538444836e-01),     IXX_REAL(-7.4764321751311447e-01),
    IXX_REAL(-8.0718270151739313e-01),     IXX_REAL(-4.8718063924993160e-01),     IXX_REAL(-3.3333333333333337e-01),
    IXX_REAL(-1.7480181875170903e-01),      IXX_REAL(4.8718063924993160e-01),     IXX_REAL(-8.5562804354527955e-01),
     IXX_REAL(3.2386650321468835e-01),      IXX_REAL(9.4263104463438019e-01),      IXX_REAL(8.0976550846447878e-02)
};


  /// Equal-weight grid: every weight is 4*pi/8 (sphere area / npts).
  /// copy_grid forms it exactly; see util/copy_grid.hpp.
  static constexpr ixx_int uniform_weight_npts = 8;
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
