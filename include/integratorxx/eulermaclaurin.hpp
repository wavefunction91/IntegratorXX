#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_EULERMACLAURIN_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_EULERMACLAURIN_HPP__

namespace IntegratorXX {

  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer
  >
  class GenerateQuadrature< EulerMaclaurin<PointType,wght_t,ContiguousContainer> > {

    using point_container  = typename EulerMaclaurin<PointType,wght_t,ContiguousContainer>::point_container;
    using weight_container = typename EulerMaclaurin<PointType,wght_t,ContiguousContainer>::weight_container;

    static auto generate_impl( const size_t nPts );

  public:

    template <typename... Args>
    inline static auto generate(Args&&... args){
      return generate_impl( std::forward<Args>(args)... );
    }

  };

  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer
  >
  auto GenerateQuadrature< EulerMaclaurin<PointType,wght_t,ContiguousContainer> >::
    generate_impl(const size_t nPts) {

    point_container   pts(nPts);
    weight_container  wghts(nPts);

    for(size_t i = 0; i < nPts; i++) {

      pts[i] = (i + 1.0) * (i + 1.0) / ((nPts - i) * (nPts - i));

      wghts[i] = 2.0 * (i + 1.0) * (nPts + 1.0) / 
        ((nPts - i) * (nPts - i) * (nPts - i));

    }; // Loop over points

    return std::tuple( std::move(pts), std::move(wghts) );

  }

}

#endif
