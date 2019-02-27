#ifndef __INCLUDED_INTEGRATORXX_INTEGRATOR_HPP__
#define __INCLUDED_INTEGRATORXX_INTEGRATOR_HPP__

#include "batch.hpp"

namespace IntegratorXX {


/*
  template <typename Func, typename RetType, typename, typename... Args>
  struct is_increment_function_sig : std::false_type { };

  template <typename Func, typename RetType, typename... Args>
  struct is_increment_function_sig< 
    Func,
    RetType,
    std::void_t<decltype(std::declval<Func>(std::declval<RetType&>(), std::declval<Args>()...))>,
    Args... > : std::true_type { }; 
*/
    
    



  template <typename QuadType>
  class Integrator1D {

    const QuadratureBatch<QuadType> qb;

  public:

    Integrator1D() = delete;
    Integrator1D( const Integrator1D& ) = default;
    Integrator1D( Integrator1D&& )      = default;

    Integrator1D( const QuadType& q, const size_t bsz = 1 ) : 
      qb(q, bsz) { }


    // TODO: SFINAE on function signature
    template <typename RetType, typename Func, typename... Args>
/*
    std::enable_if_t<
      is_increment_function_sig< Func, RetType, const typename QuadType::point_container&, const typename QuadType::weight_container&,
        Args... >::value,
      RetType
    > 
*/
    RetType
    integrate( const Func& func, RetType&& init, Args&&... args ) {

      for( auto&& [points,weights] : qb ) {
        func( init, points, weights, std::forward<Args>(args)... );
      }

      return std::forward<RetType>(init);

    }
    
    
  };


};
#endif
