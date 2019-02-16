#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_HPP__

#include <array>
#include <vector>
#include <tuple>
#include <iostream>
#include <cassert>
#include <cmath>

namespace IntegratorXX {


  template <typename T = double>
  using cartesian_pt_t = std::array<T, 3>;


  /**
   *  \brief The Quatrature Interface
   *
   *  Details the interface for the instantiation of
   *  quadrature rules
   */ 
  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer,
    typename Derived 
  >
  struct Quadrature;

  /**
   *  \brief Placeholder class to perform template specialization
   *  for quadrature generation
   */ 
  template <typename QuadratureType> class GenerateQuadrature;







  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer,
    typename Derived
  >
  class Quadrature 
  {

  public:
    using point_type = PointType; 
      ///< Storage type for the point
    using weight_type = wght_t;

    template <typename T>
    using container_type = ContiguousContainer<T>;

    using point_container = 
      container_type< point_type >;

    using weight_container = 
      container_type< weight_type >;

  protected:

    point_container  pts; ///< Integration pts
    weight_container weights; ///< Quadrature weights

    
    // Private constructors

    Quadrature(
      const point_container&  _pts,
      const weight_container& _wgt
    ): pts( _pts ), weights( _wgt ){ }

    Quadrature(
      point_container&&  _pts,
      weight_container&& _wgt
    ): pts( std::move(_pts) ), weights( std::move(_wgt) ){ }


    Quadrature(
      const std::tuple<point_container, weight_container>& tup
    ) : Quadrature( std::get<0>(tup), std::get<1>(tup) ){ };

    Quadrature(
      std::tuple<point_container, weight_container>&& tup
    ) : Quadrature( 
          std::move(std::get<0>(tup)), 
          std::move(std::get<1>(tup)) 
        ){ }


    public:


    template <typename... Args>
    Quadrature( const size_t nPts, Args&&... args ) :
      Quadrature( 
        std::move(GenerateQuadrature<Derived>::generate(nPts,std::forward<Args>(args)...)) 
      ) { }



    /**
     *  \brief Return number of integration points
     */ 
    inline auto nPts() const { return pts.size(); };
    

  };


  #define QuadratureImpl(NAME)                                                            \
  template <                                                                              \
    typename PointType,                                                                   \
    typename wght_t = double,                                                             \
    template<typename> class ContiguousContainer = std::vector                            \
  >                                                                                       \
  class NAME : public Quadrature<PointType,wght_t,ContiguousContainer,NAME<PointType>>    \
  {                                                                                       \
    using Base = Quadrature<PointType,wght_t,ContiguousContainer,NAME<PointType>>;        \
    using Base::Quadrature;                                                               \
  };



  QuadratureImpl( GaussLegendre );
  QuadratureImpl( GaussChebyshev1 );
  QuadratureImpl( GaussChebyshev2 );

};

#include "gausslegendre.hpp"
#include "gausscheby1.hpp"
#include "gausscheby2.hpp"

#endif
