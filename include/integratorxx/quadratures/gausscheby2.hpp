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
 * weights have been scaled by \f$ w_i \to w_i \sqrt{1-x_{i}^2} \f$.
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

  inline static std::tuple<point_container, weight_container>
  generate(size_t npts, point_type lo, point_type up) {

    const weight_type pi_ov_npts_p_1 = M_PI / (npts + 1);

    weight_container weights(npts);
    point_container points(npts);
    for (size_t i = 0; i < npts; ++i) {
      // The standard nodes and weights are given by
      points[i] = std::cos((i + 1) * pi_ov_npts_p_1);
      weights[i] =
          pi_ov_npts_p_1 * std::pow(std::sin((i + 1) * pi_ov_npts_p_1), 2);

      // However, since we want the rule with a unit weight factor, we
      // divide the weights by sqrt(1-x^2).
      weights[i] /= std::sqrt(1.0 - std::pow(points[i], 2));
    }

    return std::make_tuple(points, weights);
  }
};

} // namespace IntegratorXX
