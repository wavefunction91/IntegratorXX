#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_KNOWLES_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_KNOWLES_HPP__

namespace IntegratorXX::deprecated {

  /**
   *  \brief Euler-Maclaurin quadrature factory
   */ 
  template <
    typename PointType, 
    typename wght_t,
    template<typename...> class ContiguousContainer
  >
  class GenerateQuadrature< Knowles<PointType,wght_t,ContiguousContainer> > {

    using point_container  = typename Knowles<PointType,wght_t,ContiguousContainer>::point_container;
    using weight_container = typename Knowles<PointType,wght_t,ContiguousContainer>::weight_container;

    /**
     *  \brief Generate the Euler-Maclaurin quadrature rule of a specific order (impl)
     *
     *  \param[in] nPts Number of quadrature points
     *  \returns [points,weights] tuple of quadrature points and weights
     *
     */ 
    static auto generate_impl( const size_t nPts );

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
  auto GenerateQuadrature< Knowles<PointType,wght_t,ContiguousContainer> >::
    generate_impl(const size_t nPts) {

    point_container   pts(nPts);
    weight_container  wghts(nPts);

    for(size_t i = 0; i < nPts; i++) {

      const PointType x = PointType(i+1) / (nPts + 1);
      const PointType x2 = x*x;
      const PointType omx3 = 1. - x2*x;

      pts[i]   = -std::log( omx3 );
      wghts[i] = 3. * x2 / (nPts + 1) / omx3;

    }; // Loop over points

    return std::tuple( std::move(pts), std::move(wghts) );

  }

}

#endif

