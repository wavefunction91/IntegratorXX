#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_ALDRICHS_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_ALDRICHS_HPP__

namespace IntegratorXX::deprecated {

  /**
   *  \brief Euler-Maclaurin quadrature factory
   */ 
  template <
    typename PointType, 
    typename wght_t,
    template<typename...> class ContiguousContainer
  >
  class GenerateQuadrature< Ahlrichs<PointType,wght_t,ContiguousContainer> > {

    using point_container  = typename Ahlrichs<PointType,wght_t,ContiguousContainer>::point_container;
    using weight_container = typename Ahlrichs<PointType,wght_t,ContiguousContainer>::weight_container;

    /**
     *  \brief Generate the Euler-Maclaurin quadrature rule of a specific order (impl)
     *
     *  \param[in] nPts Number of quadrature points
     *  \returns [points,weights] tuple of quadrature points and weights
     *
     */ 
    static auto generate_impl( const size_t nPts, const wght_t alpha = 0.6 );

  public:

    /**
     *  \brief Generate the Euler-Maclaurin quadrature rule of a specific order (interface)
     */ 
    template <typename... Args>
    inline static auto generate(Args&&... args){
      return generate_impl( std::forward<Args>(args)... );
    }

  };

  // Implementation of Euler-Maclaurin quadrature factory
  template <
    typename PointType, 
    typename wght_t,
    template<typename...> class ContiguousContainer
  >
  auto GenerateQuadrature< Ahlrichs<PointType,wght_t,ContiguousContainer> >::
    generate_impl(const size_t nPts, const wght_t alpha) {

    point_container   pts(nPts);
    weight_container  wghts(nPts);

    const PointType pi_by_npts = M_PI / (nPts+1);
    const PointType lnt = std::log(2.);

    for(size_t i = 0; i < nPts; i++) {

      const PointType x = std::cos( (i+1) * pi_by_npts );
      const PointType y = std::pow( 1. + x, alpha ) / lnt;
      const PointType z = std::log(1. - x) - lnt;
      const PointType r = std::sqrt((1. + x) / (1. - x));

      pts[i] = y * z;
      wghts[i] = pi_by_npts * y * ( r - alpha*z/r );

    }; // Loop over points

    return std::tuple( std::move(pts), std::move(wghts) );

  }

}

#endif

