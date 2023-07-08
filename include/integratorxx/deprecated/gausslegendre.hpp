#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_GAUSSLEGENDRE_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_GAUSSLEGENDRE_HPP__

#include "util/bound_transform.hpp"
namespace IntegratorXX::deprecated {

/**
 *  \brief Gauss-Legendre quadrature factory
 */
template <typename PointType, typename wght_t,
          template <typename...> class ContiguousContainer>
class GenerateQuadrature<
    GaussLegendre<PointType, wght_t, ContiguousContainer> > {
  using point_container =
      typename GaussLegendre<PointType, wght_t,
                             ContiguousContainer>::point_container;
  using weight_container =
      typename GaussLegendre<PointType, wght_t,
                             ContiguousContainer>::weight_container;

  /**
   *  \brief Generate the Gauss-Legendre quadrature rule of a specific order
   * (impl)
   *
   *  \param[in] nPts Number of quadrature points
   *  \param[in] lo   Lower bound of the integration
   *  \param[in] up   Upper bound of the integration
   *
   *  \returns [points,weights] tuple of quadrature points and weights
   *
   */
  static auto generate_impl(const size_t nPts, const PointType lo,
                            const PointType up);

 public:
  /**
   *  \brief Generate the Gauss-Legendre quadrature rule of a specific order
   * (interface)
   */
  template <typename... Args>
  inline static auto generate(Args&&... args) {
    return generate_impl(std::forward<Args>(args)...);
  }
};

// Implementation of Gauss-Legendre quadrature factory
template <typename PointType, typename wght_t,
          template <typename...> class ContiguousContainer>
auto GenerateQuadrature<
    GaussLegendre<PointType, wght_t, ContiguousContainer> >::
    generate_impl(const size_t nPts, const PointType lo, const PointType up) {
  assert(nPts % 2 == 0);
  assert(lo != -std::numeric_limits<double>::infinity());
  assert(up != std::numeric_limits<double>::infinity());

  point_container pts(nPts);
  weight_container wghts(nPts);

  auto mid = (nPts + 1) / 2;

  const PointType eps(3.e-11);  // Convergence tolerance

  // Legendre indexing starts at 1
  for(size_t i = 1; i <= mid; ++i) {
    PointType z =
        std::cos(M_PI * (PointType(i) - 0.25) / (PointType(nPts) + 0.5));
    PointType pp(0), z1(0);

    // Iteratively determine the i-th root
    while(std::abs(z - z1) > eps) {
      PointType p1(1.), p2(0.);

      // Loop over the recurrence relation to evaluate the
      // Legendre polynomial at position z
      for(size_t j = 1; j <= nPts; ++j) {
        PointType p3 = p2;
        p2 = p1;
        p1 = ((2. * PointType(j) - 1.) * z * p2 - (PointType(j) - 1.) * p3) /
             PointType(j);
      }  // end j for

      // p1 is now the desired Legrendre polynomial. We next compute
      // pp, its derivative, by a standard relation involving also p2,
      // the polynomial of one lower order
      pp = PointType(nPts) * (z * p1 - p2) / (z * z - 1.);
      z1 = z;
      z = z1 - p1 / pp;
    }  // end while

    PointType pt = z;
    PointType wgt = 2. / (1. - z * z) / pp / pp;

    // Transform points and populate arrays
    std::tie(pt, wgt) = transform_minus_one_to_one(lo, up, pt, wgt);
    pts[i - 1] = pt;
    wghts[i - 1] = wgt;

    // Reflect the points
    pts[nPts - i] = (lo + up) - pt;
    wghts[nPts - i] = wgt;

  };  // Loop over points

  return std::tuple(std::move(pts), std::move(wghts));
}

};  // namespace IntegratorXX::deprecated

#endif
