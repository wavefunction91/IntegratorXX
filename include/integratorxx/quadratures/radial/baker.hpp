#pragma once

#include <integratorxx/quadratures/primitive/uniform.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Baker-Andzelm-Scheiner-Delley radial
 *  quadrature transformation rules.
 *
 *  Reference:
 *  J. Chem. Phys. 101, 8894 (1994)
 *  DOI: https://doi.org/10.1063/1.468081
 */
class BakerRadialTraits {

  double Rmax_;

public:


/**
 *  Specify Baker-Andzelm-Scheiner-Delley quadrature parameters.
 *
 *  @param[in] Rmax  Maximum radial extent, constant 12 in Equation (12b) of J. Chem. Phys. 101, 8894 (1994).
 */
  BakerRadialTraits(double Rmax = 12.0) :
    Rmax_(Rmax) { }

  /**
   *  @brief Transformation rule for the Baker-Andzelm-Scheiner-Delley radial quadratures
   *  
   *  @param[in] x Point in (-1,1)
   *  @return    r = Rfac * log [1 - x_i^2 ]
   */
  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    const auto log_term = std::log(1.0 - x*x);
    return Rfac_ * log_term * log_term;
  };


  /**
   *  @brief Jacobian of the Baker-Andzelm-Scheiner-Delley radial transformations
   *
   *  @param[in] x Point in (-1,1)
   *  @returns   dr/dx (see `radial_transform`)
   */
  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    const auto log_term = std::log(1.0 - x*x);
    return - 2.0 * Rfac_ * x / (1.0 - x*x);
  }

};


/**
 *  @brief Implementation of the Baker-Andzelm-Scheiner-Delley radial quadrature.
 *
 *  Taken as the convolution of the uniform trapezoidal quadrature
 *  with the Baker-Andzelm-Scheiner-Delley radial transformation. See
 *  BakerRadialTraits for details.
 *
 *  Reference:
 *  J. Chem. Phys. 101, 8894 (1994)
 *  DOI: https://doi.org/10.1063/1.468081
 *
 *  @tparam PointType  Type describing the quadrature points  
 *  @tparam WeightType Type describing the quadrature weights 
 */
template <typename PointType, typename WeightType>
using Baker = RadialTransformQuadrature<
  UniformTrapezoid<PointType, WeightType>,
  BakerRadialTraits 
>;

}
