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
 *  This method eliminates the sqrt(1-x^2) weight from the Chebyshev
 *  quadrature rule by a change of variables. The nodes and weights
 *  are still analytic.
 *
 *  References:
 *  Comput. Phys. Commun. 70, 271 (1992)
 *  doi:10.1016/0010-4655(92)90192-2
 *
 *  J. Chem. Phys. 100, 6520–6534 (1994)
 *  doi:10.1063/1.467061
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */

template <typename PointType, typename WeightType>
class GaussChebyshev2Modified
    : public Quadrature<GaussChebyshev2Modified<PointType, WeightType>> {
  using base_type = Quadrature<GaussChebyshev2Modified<PointType, WeightType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussChebyshev2Modified(size_t npts, point_type lo, point_type up)
      : base_type(npts, lo, up) {}

  GaussChebyshev2Modified(const GaussChebyshev2Modified &) = default;
  GaussChebyshev2Modified(GaussChebyshev2Modified &&) noexcept = default;
};

template <typename PointType, typename WeightType>
struct quadrature_traits<GaussChebyshev2Modified<PointType, WeightType>> {
  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts, point_type lo, point_type up) {
    const weight_type oonpp = 1.0 / (npts + 1);

    point_container points(npts);
    weight_container weights(npts);
    for(size_t i = 1; i <= npts; ++i) {
      auto sine = sin(i * M_PI * oonpp);
      auto sinesq = sine * sine;
      auto cosine = cos(i * M_PI * oonpp);

      points[i - 1] = 1.0 - 2.0 * i * oonpp +
                      M_2_PI * (1.0 + 2.0 / 3.0 * sinesq) * cosine * sine;
      weights[i - 1] = 16.0 / 3.0 / (npts + 1.0) * sinesq * sinesq;
    }

    return std::make_tuple(points, weights);
  }
};

}  // namespace IntegratorXX
