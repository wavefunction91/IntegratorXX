#pragma once

#include <integratorxx/quadrature.hpp>
#include <iostream>

namespace IntegratorXX {

/**
 *  @brief Implementation of Gauss-Chebyshev quadrature of the first kind.
 *
 *  This quadrature is originally derived for integrals of the type
 *
 *  \f$ \displaystyle \int_{-1}^{1} \frac {f(x)} {\sqrt{1-x^2}} {\rm d}x \f$.
 *
 *  The original nodes and weights from the rule are given by
 *
 *  \f{eqnarray*}{ x_{i} = & \cos \left( \displaystyle \frac {2i-1} {2n} \pi
 * \right) \\ w_{i} = & \displaystyle \frac \pi n \f}
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
class GaussChebyshev1
    : public Quadrature<GaussChebyshev1<PointType, WeightType>> {

  using base_type =
      Quadrature<GaussChebyshev1<PointType, WeightType>>;

public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussChebyshev1(size_t npts, point_type lo, point_type up)
      : base_type(npts, lo, up) {}

  GaussChebyshev1(const GaussChebyshev1 &) = default;
  GaussChebyshev1(GaussChebyshev1 &&) noexcept = default;
};

template <typename PointType, typename WeightType>
struct quadrature_traits<GaussChebyshev1<PointType, WeightType>> {

  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container>
  generate(size_t npts, point_type lo, point_type up) {

    point_container points(npts);
    weight_container weights(npts);

    const weight_type pi_ov_npts = M_PI / npts;
    const weight_type pi_ov_2npts = pi_ov_npts / 2;

    for (size_t i = 0; i < npts; ++i) {
      const auto ti = (2.0 * (i + 1) - 1.) * pi_ov_2npts;
      // The standard nodes and weights are given by
      const auto xi = std::cos(ti); 
      auto wi = pi_ov_npts;

      // However, as we're integrating f(x) not \frac{f(x)}{\sqrt{1 -
      // x^2}}, we must factor the \sqrt{1-x^2} into the weight
      wi *= std::sqrt(1. - xi*xi);

      // Store into memory
      points[i]  = xi;
      weights[i] = wi;
    }

    return std::make_tuple(points, weights);
  }
};
} // namespace IntegratorXX
