#pragma once

#include <integratorxx/quadratures/primitive/uniform.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Murray-Handy-Laming radial quadrature
 *  transformation rules.
 *
 *  Reference:
 *  Molecular Physics, 78:4, 997-1014,
 *  DOI: https://doi.org/10.1080/00268979300100651
 *
 *  @tparam M Integer to modulate the MHL transformation. 
 *            Typically taken to be 2.
 */
template <size_t M>
class MurrayHandyLamingRadialTraits : public RadialTraits {

  size_t npts_; ///< Number of grid points
  double R_; ///< Radial scaling factor
  

public:

  MurrayHandyLamingRadialTraits(size_t npts, double R = 1.0) : npts_(npts), R_(R) {}

  size_t npts() const noexcept { return npts_; }

  /**
   *  @brief Transformation rule for the MHL radial quadrature
   *  
   *  @param[in] x Point in (0,1)
   *  @return    r = R * (x / (1-x))^M
   */
  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    return R_ * std::pow( x / (1.0 - x), M );
  }

  /**
   *  @brief Jacobian of the MHL radial transformation
   *
   *  @param[in] x Point in (0,1)
   *  @returns   dr/dx (see `radial_transform`)
   */
  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    return R_ * M * std::pow(x, M-1) / std::pow(1.0 - x, M+1);
  }

};


/**
 *  @brief Implementation of the Murray-Handy-Laming radial quadrature.
 *
 *  Taken as the convolution of the Uniform (Trapezoid) quadrature with
 *  the MHL radial transformation. See MurrayHandyLamingRadialTraits for
 *  details.
 *
 *  Suitable for integrands which tend to zero as their argument tends 
 *  to 0 and inf. Tailored for radial integrands, i.e. r^2 * f(r), with 
 *  lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  Molecular Physics, 78:4, 997-1014,
 *  DOI: https://doi.org/10.1080/00268979300100651
 *
 *  @tparam PointType  Type describing the quadrature points  
 *  @tparam WeightType Type describing the quadrature weights 
 */
template <typename PointType, typename WeightType>
using MurrayHandyLaming = RadialTransformQuadrature<
  UniformTrapezoid<PointType,WeightType>,
  MurrayHandyLamingRadialTraits<2>
>;

namespace detail {

template <typename QuadType>
static constexpr bool is_mhl_v = std::is_same_v<
  QuadType, 
  MurrayHandyLaming<typename QuadType::point_type, typename QuadType::weight_type>
>;

}
}

