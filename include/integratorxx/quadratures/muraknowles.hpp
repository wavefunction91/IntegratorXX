#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/uniform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Mura-Knowles radial quadrature.
 *
 *  Generates a quadrature on the bounds (0, inf). Suitable for integrands
 *  which tend to zero as their argument tends to 0 and inf. Tailored for
 *  radial integrands, i.e. r^2 * f(r), with lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  J. Chem. Phys. 104, 9848 (1996);
 *  DOI: https://doi.org/10.1063/1.471749
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */

template <typename PointType, typename WeightType>
class MuraKnowles : public Quadrature<MuraKnowles<PointType, WeightType>> {
  using base_type = Quadrature<MuraKnowles<PointType, WeightType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  /**
   *  @brief Construct the Mura-Knowles radial quadrature
   *
   *  @param[in] npts Number of quadrature points to generate
   *  @param[in] R    Radial scaling factor. Table for suggested
   *                  values is given in the original reference.
   */
  MuraKnowles(size_t npts, weight_type R = 1.) : base_type(npts, R) {}

  MuraKnowles(const MuraKnowles&) = default;
  MuraKnowles(MuraKnowles&&) noexcept = default;
};

/**
 *  @brief Quadrature traits for the Mura-Knowles quadrature
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */
template <typename PointType, typename WeightType>
struct quadrature_traits<
    MuraKnowles<PointType, WeightType>,
    std::enable_if_t<std::is_floating_point_v<PointType> and
                     std::is_floating_point_v<WeightType>>> {
  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  /**
   *  @brief Generator for the Mura-Knowles quadrature
   *
   *  @param[in] npts Number of quadrature points to generate
   *  @param[in] R    Radial scaling factor. Table for suggested
   *                  values is given in the original reference.
   *
   *  @returns Tuple of quadrature points and weights
   */
  inline static std::tuple<point_container, weight_container> generate(
      size_t npts, weight_type R) {
    point_container points(npts);
    weight_container weights(npts);

    using base_quad_traits =
        quadrature_traits<UniformTrapezoid<PointType, WeightType>>;

    /*
     * Generate uniform trapezoid points on [0,1]
     * ux(j) = j/(m-1)
     * uw(j) = 1/(m-1)
     * m = npts + 2, j in [0, npts+2)
     */
    auto [ux, uw] = base_quad_traits::generate(npts + 2, 0., 1.);

    /*
     * Perform Mura and Knowles transformation
     *
     * Original reference:
     * J. Chem. Phys. 104, 9848 (1996);
     * DOI: https://doi.org/10.1063/1.471749
     *
     * Closed form for points / weights obtained from:
     * Journal of Computational Chemistry, 24: 732–740, 2003
     * DOI: https://doi.org/10.1002/jcc.10211
     *
     * x(i) = ux(i+1)
     * r(i) = - log(1 - x(i)^3)
     * w(i) = uw(i+1) * 3 * x(i)^2 / ( 1 - x(i)^3 )
     * i in [0, npts)
     *
     * XXX: i+1 offset on trapezoid points ignores enpoints of trapezoid
     *      quadrature (f(r) = 0 with r in {0,inf})
     */
    for(size_t i = 0; i < npts; ++i) {
      const auto xi = ux[i + 1];
      const auto one_m_xi3 = 1. - xi * xi * xi;
      points[i] = -R * std::log(one_m_xi3);
      weights[i] = R * uw[i + 1] * 3 * xi * xi / one_m_xi3;
    }

    return std::make_tuple(points, weights);
  }
};

}  // namespace IntegratorXX
