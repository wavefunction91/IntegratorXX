#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of Gauss-Chebyshev quadrature of the third kind.
 *
 *  This quadrature is originally derived for integrals of the type
 *
 *  \f$ \displaystyle \int_{0}^{1} f(x) \sqrt{\frac {x} {1-x}} {\rm d}x \f$.
 *
 *  The original nodes and weights from the rule are given by
 *
 *  \f{eqnarray*}{ x_{i} = & \cos^2 \displaystyle \frac {(2i-1) \pi} {2(2n+1)}
 * \\ w_{i} = & \displaystyle \frac {2 \pi} {2n+1} x_i \f}
 *
 *  Reference:
 *  Abramowitz and Stegun, Handbook of Mathematical Functions with
 *  Formulas, Graphs, and Mathematical Tables, Tenth Printing,
 *  December 1972, p. 889
 *
 *  To transform the rule to the form \f$ \int_{0}^{1} f(x) {\rm d}x \f$, the
 * weights have been scaled by \f$ w_i \to w_i \sqrt{\frac {1-x_{i}} {x_i}} \f$.
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */

template <typename PointType, typename WeightType>
class GaussChebyshev3
    : public Quadrature<GaussChebyshev3<PointType, WeightType>> {
  using base_type = Quadrature<GaussChebyshev3<PointType, WeightType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussChebyshev3(size_t npts, point_type lo, point_type up)
      : base_type(npts, lo, up) {}

  GaussChebyshev3(const GaussChebyshev3 &) = default;
  GaussChebyshev3(GaussChebyshev3 &&) noexcept = default;
};

template <typename PointType, typename WeightType>
struct quadrature_traits<GaussChebyshev3<PointType, WeightType>> {
  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts, point_type lo, point_type up) {
    const weight_type pi_ov_2n_p_1 = M_PI / (2 * npts + 1);

    weight_container weights(npts);
    point_container points(npts);
    for(size_t idx = 0; idx < npts; ++idx) {
      // Transform index here for two reasons: the mathematical
      // equations are for 1 <= i <= n, and the nodes are generated in
      // decreasing order. This generates them in the right order
      size_t i = npts - idx;

      // The standard nodes and weights are given by
      const auto ti = 0.5 * (2 * i - 1) * pi_ov_2n_p_1;
      const auto cti = std::cos(ti);
      const auto xi = cti * cti;  // cos^2(t)
      auto wi = 2.0 * pi_ov_2n_p_1 * xi;

      // However, since we want the rule with a unit weight factor, we
      // divide the weights by sqrt(x/(1-x)).
      wi *= std::sqrt((1.0 - xi) / xi);

      // Finally, convert from [0,1] to [-1,1]
      points[idx] = 2 * xi - 1.0;
      weights[idx] = 2.0 * wi;
    }

    return std::make_tuple(points, weights);
  }
};

}  // namespace IntegratorXX
