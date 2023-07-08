#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_EULERMACLAURIN_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_EULERMACLAURIN_HPP__

namespace IntegratorXX::deprecated {

/**
 *  \brief Euler-Maclaurin quadrature factory
 */
template <typename PointType, typename wght_t,
          template <typename...> class ContiguousContainer>
class GenerateQuadrature<
    EulerMaclaurin<PointType, wght_t, ContiguousContainer> > {
  using point_container =
      typename EulerMaclaurin<PointType, wght_t,
                              ContiguousContainer>::point_container;
  using weight_container =
      typename EulerMaclaurin<PointType, wght_t,
                              ContiguousContainer>::weight_container;

  /**
   *  \brief Generate the Euler-Maclaurin quadrature rule of a specific order
   * (impl)
   *
   *  \param[in] nPts Number of quadrature points
   *  \returns [points,weights] tuple of quadrature points and weights
   *
   */
  static auto generate_impl(const size_t nPts);

 public:
  /**
   *  \brief Generate the Euler-Maclaurin quadrature rule of a specific order
   * (interface)
   */
  template <typename... Args>
  inline static auto generate(Args&&... args) {
    return generate_impl(std::forward<Args>(args)...);
  }
};

// Implementation of Euler-Maclaurin quadrature factory
template <typename PointType, typename wght_t,
          template <typename...> class ContiguousContainer>
auto GenerateQuadrature<EulerMaclaurin<
    PointType, wght_t, ContiguousContainer> >::generate_impl(const size_t
                                                                 nPts) {
  point_container pts(nPts);
  weight_container wghts(nPts);

  for(size_t i = 0; i < nPts; i++) {
    pts[i] = (i + 1.0) * (i + 1.0) / ((nPts - i) * (nPts - i));

    wghts[i] =
        2.0 * (i + 1.0) * (nPts + 1.0) / ((nPts - i) * (nPts - i) * (nPts - i));

  };  // Loop over points

  return std::tuple(std::move(pts), std::move(wghts));
}

}  // namespace IntegratorXX::deprecated

#endif
