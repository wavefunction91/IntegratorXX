#pragma once

#include <integratorxx/quadratures/primitive/gausschebyshev2.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Becke radial quadrature
 *  transformation rule.
 *
 *  Reference:
 *  J. Chem. Phys. 88, 2547 (1988)
 *  DOI: https://doi.org/10.1063/1.454033
 */
class BeckeRadialTraits {
  double R_;

 public:
  /**
   *  Specify Becke quadrature parameters
   *
   *  Equation (25) of J. Chem. Phys. 88, 2547 (1988)
   *
   *  @param[in] R     Radial scaling factor
   */
  BeckeRadialTraits(double R = 1.0) : R_(R) {}

  /**
   *  @brief Transformation rule for the Becke radial quadratures
   *
   *  @param[in] x Point in (-1,1)
   *  @return    r = (1+x)/(1-x)
   */
  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    return R_ * (1.0 + x) / (1.0 - x);
  };

  /**
   *  @brief Jacobian of the Becke radial transformations
   *
   *  @param[in] x Point in (-1,1)
   *  @returns   dr/dx (see `radial_transform`)
   */
  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    return R_ * 2.0 / std::pow(1.0 - x, 2);
  }
};

/**
 *  @brief Implementation of the Becke radial quadrature.
 *
 *  Taken as the convolution of the Gauss-Chebyshev (second kind)
 *  quadrature with the Becke radial transformation. See
 *  BeckeRadialTraits for details.
 *
 * Suitable for integrands which tend to zero as their argument tends
 f(r), with * lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  J. Chem. Phys. 88, 2547 (1988)
 *  DOI: https://doi.org/10.1063/1.454033
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */
template <typename PointType, typename WeightType>
using Becke = RadialTransformQuadrature<GaussChebyshev2<PointType, WeightType>,
                                        BeckeRadialTraits>;


namespace detail {

template <typename QuadType>
static constexpr bool is_becke_v = std::is_same_v<
  QuadType, 
  Becke<typename QuadType::point_type, typename QuadType::weight_type>
>;

}
}  // namespace IntegratorXX
