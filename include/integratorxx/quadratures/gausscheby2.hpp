#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of Gauss-Chebyshev quadrature of the second kind.
 *
 *  This quadrature is originally derived for integrals of the type
 *
 *  \f$ \displaystyle \int_{-1}^{1} f(x) \sqrt{1-x^2} {\rm d}x \f$.
 *
 *  The original nodes and weights from the rule are given by
 *
 *  \f{eqnarray*}{ x_{i} = & \cos \left( \displaystyle \frac {i} {n+1} \pi
 * \right) \\ w_{i} = & \displaystyle \frac \pi {n+1} \sin^2 \left( \frac {i}
 * {n+1} \pi \right) \f}
 *
 *  Reference:
 *  Abramowitz and Stegun, Handbook of Mathematical Functions with
 *  Formulas, Graphs, and Mathematical Tables, Tenth Printing,
 *  December 1972, p. 889
 *
 *  To transform the rule to the form \f$ \int_{-1}^{1} f(x) {\rm d}x \f$, the
 * weights have been scaled by \f$ w_i \to w_i / \sqrt{1-x_{i}^2} \f$. The
 * result is then simplified using \f$ \sqrt{1-x_i^2} = \sin \left(
 * \displaystyle \frac {i} {n+1} \pi \right) \f$.
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */

template <typename PointType, typename WeightType>
class GaussChebyshev2
    : public Quadrature<GaussChebyshev2<PointType, WeightType>> {
  using base_type = Quadrature<GaussChebyshev2<PointType, WeightType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussChebyshev2(size_t npts, point_type lo, point_type up)
      : base_type(npts, lo, up) {}

  GaussChebyshev2(const GaussChebyshev2 &) = default;
  GaussChebyshev2(GaussChebyshev2 &&) noexcept = default;
};

template <typename PointType, typename WeightType>
struct quadrature_traits<GaussChebyshev2<PointType, WeightType>> {
  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;
 
  inline static constexpr bool bound_inclusive = false;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts, point_type lo, point_type up) {
    const weight_type pi_ov_npts_p_1 = M_PI / (npts + 1);

    weight_container weights(npts);
    point_container points(npts);
    for(size_t idx = 0; idx < npts; ++idx) {
      // Transform index here for two reasons: the mathematical
      // equations are for 1 <= i <= n, and the nodes are generated in
      // decreasing order. This generates them in the right order
      size_t i = npts - idx;

      // The standard nodes are given by
      const auto ti = i * pi_ov_npts_p_1;
      const auto xi = std::cos(ti);

      // The quadrature weight, transformed to unit weight factor in
      // [-1,1] is given by (see comments above for explanation)
      const auto wi = pi_ov_npts_p_1 * std::sqrt(1.0 - xi * xi);

      // Store into memory
      points[idx] = xi;
      weights[idx] = wi;
    }

    return std::make_tuple(points, weights);
  }

  inline static std::tuple<point_container, weight_container>
  generate(size_t npts) {
    return generate(npts, -1.0, 1.0);
  }

};

}  // namespace IntegratorXX
