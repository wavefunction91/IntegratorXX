#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_GAUSSCHEBY2_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_GAUSSCHEBY2_HPP__

namespace IntegratorXX {

  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer
  >
  class GenerateQuadrature< GaussChebyshev2<PointType,wght_t,ContiguousContainer> > {

    using point_container  = typename GaussLegendre<PointType,wght_t,ContiguousContainer>::point_container;
    using weight_container = typename GaussLegendre<PointType,wght_t,ContiguousContainer>::weight_container;

    static auto generate_impl( const size_t nPts, const PointType lo, const PointType up );

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
  auto GenerateQuadrature< GaussChebyshev2<PointType,wght_t,ContiguousContainer> >::
    generate_impl(const size_t nPts, const PointType lo, const PointType up) {

    point_container   pts(nPts);
    weight_container  wghts(nPts);

    for(size_t i = 0; i < nPts; i++) {

      // Generate raw points and wghts on (-1,1)
      PointType pt = std::cos(M_PI * (i+1) / (nPts +1));
      PointType wgt = M_PI / (nPts + 1) * 
        std::pow(
            std::sin(M_PI * (i+1) / (nPts +1)),
            2.0
        );

      // As we're integrating f(x) not f(x)\sqrt{1 - x^2}
      // we must factor the \sqrt{1-x^2} into the wghts
      wgt /= std::sqrt(1 - pt * pt);

      // Transform points and populate arrays
      //std::tie(pts[i],wghts[i]) = unitBoundTransform(lowBound,upBound,pt,wgt);
      pts[i] = pt;
      wghts[i] = wgt;

    }; // Loop over points

    return std::tuple( std::move(pts), std::move(wghts) );

  }

}

#endif
