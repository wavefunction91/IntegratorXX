#pragma once

namespace IntegratorXX {
namespace LebedevLaikovGrids {

/**
 *  \brief Lebedev-Laikov Quadrature specification for Order = 6
 * 
 */
template <typename T>
struct lebedev_laikov_6 {

  static constexpr std::array<cartesian_pt_t<ixx_real>,6> points = {
      IXX_REAL(1.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),
     IXX_REAL(-1.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),
      IXX_REAL(0.000000000000000e+00),      IXX_REAL(1.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),
      IXX_REAL(0.000000000000000e+00),     IXX_REAL(-1.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),
      IXX_REAL(0.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),      IXX_REAL(1.000000000000000e+00),
      IXX_REAL(0.000000000000000e+00),      IXX_REAL(0.000000000000000e+00),     IXX_REAL(-1.000000000000000e+00)
  };


  static constexpr std::array<ixx_real,6> weights = {
        IXX_REAL(1.666666666666667e-01),
        IXX_REAL(1.666666666666667e-01),
        IXX_REAL(1.666666666666667e-01),
        IXX_REAL(1.666666666666667e-01),
        IXX_REAL(1.666666666666667e-01),
        IXX_REAL(1.666666666666667e-01)
  };
};
}
}
