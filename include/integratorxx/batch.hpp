#ifndef __INCLUDED_INTEGRATORXX_BATCH_HPP__
#define __INCLUDED_INTEGRATORXX_BATCH_HPP__

#include <iterator>

#include "quadrature.hpp"

namespace IntegratorXX {

  template <typename QuadratureType>
  class QuadratureBatch {

    const QuadratureType& quad; ///< Base quadrature
    const size_t          batch_sz; ///< Batch size

    class iterator {

      using point_type  = typename QuadratureType::point_type;
      using weight_type = typename QuadratureType::weight_type;
      using point_container  = typename QuadratureType::point_container;
      using weight_container = typename QuadratureType::weight_container;

      const QuadratureType& quad; ///< Base quadrature

      size_t cur_indx;
      size_t batch_sz;

    public:

      using value_type  = 
        std::tuple< point_container, weight_container >;
      using pointer           = const value_type*;
      using reference         = value_type&&;
      using difference_type   = size_t;
      using iterator_catagory = std::forward_iterator_tag;
      
      iterator() = delete;

      iterator( 
        const QuadratureType &q, 
        const size_t indx,
        const size_t bsz = 1 
      ) : quad(q), cur_indx(indx), batch_sz(bsz) { }

      iterator& operator++(){ cur_indx += std::min(batch_sz,quad.nPts()-cur_indx); return *this; }
      iterator  operator++(int){ iterator retval = *this; ++(*this); return retval; }

      // TODO, check that quads are the same through hash
      bool operator==(iterator other){ return cur_indx == other.cur_indx; }
      bool operator!=(iterator other){ return not (*this == other); }

      auto operator*(){ 

        size_t ncpy = 
          std::min(batch_sz, quad.nPts() - cur_indx);

        point_container  pts( ncpy );
        weight_container wgt( ncpy );

        std::copy( quad.points().begin() + cur_indx,
                   quad.points().begin() + cur_indx + ncpy,
                   pts.begin() );

        std::copy( quad.weights().begin() + cur_indx,
                   quad.weights().begin() + cur_indx + ncpy,
                   wgt.begin() );

        return std::tuple( std::move(pts),std::move(wgt) ); 

      }

    };

  public:

    using point_type  = typename QuadratureType::point_type;
    using weight_type = typename QuadratureType::weight_type;
    using point_container = 
      typename QuadratureType::point_container;
    using weight_container = 
      typename QuadratureType::weight_container;

    QuadratureBatch( 
      const QuadratureType &q, 
      const size_t bsz = 1  
    ): batch_sz(bsz), quad(q){ };


    iterator begin(){ return iterator(quad,0, batch_sz); }
    iterator end(){   return iterator(quad,quad.nPts(),batch_sz); }

    inline size_t n_batches() const {
      if( batch_sz > quad.nPts() ) return 1;

      size_t n = quad.nPts() / batch_sz;
      if( quad.nPts() % batch_sz ) n++;
      return n;
    }

  };

};

#endif
