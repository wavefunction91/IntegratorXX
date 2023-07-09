#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/gausscheby2.hpp>
#include <integratorxx/quadratures/radial_transform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Treutler-Ahlrichs M4 radial 
 *  quadrature transformation rules.
 *
 *  Reference:
 *  J. Chem. Phys. 102, 346 (1995)
 *  DOI: https://doi.org/10.1063/1.469408
 */
class TreutlerAhlrichsM4RadialTraits {

  double R_;
  double alpha_;

public:

  inline static constexpr double a     = 1.0;
  inline static constexpr double ln_2  = 0.693147180559945309417232;

  TreutlerAhlrichsM4RadialTraits(double R = 1.0, double alpha = 0.6) :
    R_(R), alpha_(alpha) { }

  /**
   *  @brief Transformation rule for the TA M4 radial quadrature
   *
   *  Equation (19) of J. Chem. Phys. 102, 346 (1995)
   *  
   *  @param[in] x Point in (-1,1)
   *  @return    r = (a+x)^alpha * log((a+1)/(1-x)) / ln(2) 
   */
  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    const auto pow_term = std::pow(a + x, alpha_);
    const auto log_term = std::log((a + 1.0) / (1.0 - x));
    return R_ * pow_term * log_term / ln_2;
  };


  /**
   *  @brief Jacobian of the TA M4 radial transformation
   *
   *  @param[in] x Point in (-1,1)
   *  @returns   dr/dx (see `radial_transform`)
   */
  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    const auto pow_term = std::pow(a + x, alpha_);
    const auto log_term = std::log((a + 1.0) / (1.0 - x));
    return R_ * pow_term / ln_2 * ( alpha_ * log_term / (a+x) + (1./(1.-x)) );
  }

};


/**
 *  @brief Implementation of the Treutler-Ahlrichs M4 radial quadrature.
 *
 *  Taken as the convolution of the Gauss-Chebyshev (second kind) quadrature 
 *  with the TA M4 radial transformation. See TreutlerAhlrichsM4RadialTraits for
 *  details.
 *
 *  Suitable for integrands which tend to zero as their argument tends 
 *  to 0 and inf. Tailored for radial integrands, i.e. r^2 * f(r), with 
 *  lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  J. Chem. Phys. 102, 346 (1995)
 *  DOI: https://doi.org/10.1063/1.469408
 *
 *  @tparam PointType  Type describing the quadrature points  
 *  @tparam WeightType Type describing the quadrature weights 
 */
template <typename PointType, typename WeightType>
using TreutlerAhlrichs = RadialTransformQuadrature<
  GaussChebyshev2<PointType, WeightType>,
  TreutlerAhlrichsM4RadialTraits 
>;

}
