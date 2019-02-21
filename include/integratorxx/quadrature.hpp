#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_HPP__

#include <array>
#include <vector>
#include <tuple>
#include <iostream>
#include <cassert>
#include <cmath>

namespace IntegratorXX {


  /**
   *  \brief 3D Cartesian Point
   */ 
  template <typename T = double>
  using cartesian_pt_t = std::array<T, 3>;


  /**
   *  \brief The Quatrature Interface
   *
   *  Details the interface for the instantiation of
   *  quadrature rules. All definitions for quadrature
   *  rules should derive from this class
   */ 
  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer,
    typename Derived 
  >
  class Quadrature;

  /**
   *  \brief Placeholder class to perform template specialization
   *  for quadrature generation
   */ 
  template <typename QuadratureType> class GenerateQuadrature;







  // Instantiation of Quadrature base class
  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer
  >
  class QuadratureBase
  {

  public:
    using point_type = PointType; 
      ///< Storage type for the points
    using weight_type = wght_t;
      ///< Storage type for the weights

    template <typename T>
    using container_type = ContiguousContainer<T>;
      ///< Contiguous container type

    using point_container = 
      container_type< point_type >;
      ///< Contiguous point storage

    using weight_container = 
      container_type< weight_type >;
      ///< Contiguous weight storage

  protected:

    point_container  pts; ///< Integration points
    weight_container wghts; ///< Quadrature weights

    
    // Internal constructors

    // Construct from const points / weights
    QuadratureBase(
      const point_container&  _pts,
      const weight_container& _wgt
    ): pts( _pts ), wghts( _wgt ){ }

    // Construct from rvalue points / weights
    QuadratureBase(
      point_container&&  _pts,
      weight_container&& _wgt
    ): pts( std::move(_pts) ), wghts( std::move(_wgt) ){ }


    // Construct from const tuple of points / weights
    QuadratureBase(
      const std::tuple<point_container, weight_container>& tup
    ) : QuadratureBase( std::get<0>(tup), std::get<1>(tup) ){ };

    // Construct from rvalue tuple of points / weights
    QuadratureBase(
      std::tuple<point_container, weight_container>&& tup
    ) : QuadratureBase( 
          std::move(std::get<0>(tup)), 
          std::move(std::get<1>(tup)) 
        ){ }


    public:

    /**
     *  \brief Return number of integration points
     */ 
    inline auto nPts() const { return pts.size(); };
    
    /**
     *  \brief Return const reference to the quadrature points.
     */ 
    inline const point_container&  points()  const { return pts; };

    /**
     *  \brief Return const reference to the quadrature weights.
     */ 
    inline const weight_container& weights() const { return wghts; };

  }; // class QuadratureBase







  // Instantiation of Quadrature base class
  template <
    typename PointType, 
    typename wght_t,
    template<typename> class ContiguousContainer,
    typename Derived
  >
  class Quadrature : public QuadratureBase< PointType, wght_t, ContiguousContainer >
  {

    public:

    Quadrature() = delete; // no default ctor

    /**
     *  \brief Construct a Quadrature object
     *
     *  Delagate Quadrature constuction to specific
     *  GenerateQuadrature<Derived> implementation.
     */ 
    template <typename... Args>
    Quadrature( const size_t nPts, Args&&... args ) :
      QuadratureBase< PointType, wght_t, ContiguousContainer >( 
        std::move(GenerateQuadrature<Derived>::generate(nPts,std::forward<Args>(args)...)) 
      ) { }

  }; // class Quadrature


  // Define quadrature rules

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


  /**
   *  \brief Gauss-Legendre quadrature rule
   */ 
  QuadratureImpl( GaussLegendre );

  /**
   *  \brief Gauss-Chebyshev (1st kind) quadrature rule
   */ 
  QuadratureImpl( GaussChebyshev1 );

  /**
   *  \brief Gauss-Chebyshev (2nd kind) quadrature rule
   */ 
  QuadratureImpl( GaussChebyshev2 );

  /**
   *  \brief Euler-Maclaurin (1st kind) quadrature rule
   */ 
  QuadratureImpl( EulerMaclaurin );


  /**
   *  \brief Lebedev-Laikov spherical quadrature rule
   */ 
  template <                  
    typename PointType,       
    typename wght_t = double, 
    template<typename> class ContiguousContainer = std::vector
  >                                                           
  class Lebedev : public Quadrature<cartesian_pt_t<PointType>,wght_t,ContiguousContainer,Lebedev<PointType>>
  {                                                                                   
    using Base = Quadrature<cartesian_pt_t<PointType>,wght_t,ContiguousContainer,Lebedev<PointType>>;    
    using Base::Quadrature;                                                           
  };


};

// Instantiations of quadrature rules
#include "gausslegendre.hpp"
#include "gausscheby1.hpp"
#include "gausscheby2.hpp"
#include "eulermaclaurin.hpp"
#include "lebedev.hpp"

#endif
