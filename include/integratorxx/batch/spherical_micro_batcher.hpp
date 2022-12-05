#pragma once

#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/type_traits.hpp>
#include <integratorxx/batch/octree_partition.hpp>
#include <integratorxx/batch/hilbert_partition.hpp>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>


namespace IntegratorXX {



template <typename PointContainer, typename WeightContainer, 
  typename GridPartitioner >
class SphericalMicroBatcher {


  using quad_type   = SphericalQuadratureBase<PointContainer,WeightContainer>;
  static_assert( detail::is_cartesian_quadrature_v<quad_type> );


  using point_type  = typename quad_type::point_type;
  using weight_type = typename quad_type::weight_type;
  using point_container  = typename quad_type::point_container;
  using weight_container = typename quad_type::weight_container;

  using index_type      = size_t;
  using index_container = std::vector<index_type>;


  using point_iterator  = typename point_container::iterator;
  using weight_iterator = typename weight_container::iterator;
  using index_iterator  = typename index_container::iterator;

  using const_point_iterator  = typename point_container::const_iterator;
  using const_weight_iterator = typename weight_container::const_iterator;
  using const_index_iterator  = typename index_container::const_iterator;


  size_t max_batch_sz_;

  std::shared_ptr<quad_type> quad_;
  index_container            partition_idx_;

  inline index_container generate_batches() {
    return GridPartitioner::generate_batches(*quad_, max_batch_sz_);
  };


  struct iterator {
  
    using value_type        = 
      std::tuple<point_type,point_type,point_container,weight_container>;
    using different_type    = size_t;
    using pointer           = value_type*;
    using reference         = value_type&;
    using iterator_catagory = std::input_iterator_tag;

    index_iterator idx_it;
    point_iterator  point_begin;
    weight_iterator weight_begin;

    iterator() = delete;
    iterator( index_iterator it, point_iterator pb, weight_iterator wb ) : 
      idx_it(it), point_begin(pb), weight_begin(wb) { }

    iterator& operator++(){ idx_it++; return *this; }
    iterator  operator++(int) { 
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    iterator& operator+(int i) {
      idx_it += i;
      return *this;
    } 

    bool operator==( iterator other ){ return idx_it == other.idx_it; }
    bool operator!=( iterator other ){ return not (*this == other);   }



    auto range() {

      const auto idx      = *idx_it;
      const auto idx_next = *(idx_it+1);

      const auto npts = idx_next - idx;
    
      return std::make_tuple(
        npts,
        point_begin + idx,
        point_begin + idx + npts,
        weight_begin + idx,
        weight_begin + idx + npts 
      );

    }

    value_type operator*() {

      auto [npts,pb,pe,wb,we] = range();
      auto [box_lo, box_up]   = detail::get_box_bounds_points(pb, pe);

      return std::make_tuple(
        box_lo,
        box_up,
        std::vector<point_type>( pb, pe ),
        std::vector<weight_type>( wb, we )
      );

    }

    
  };



  struct const_iterator {
  
    using value_type        = 
      std::tuple<point_type,point_type,point_container,weight_container>;
    using different_type    = size_t;
    using pointer           = value_type*;
    using reference         = value_type&;
    using iterator_catagory = std::input_iterator_tag;

    const_index_iterator idx_it;
    const_point_iterator  point_begin;
    const_weight_iterator weight_begin;

    const_iterator() = delete;
    const_iterator( const_index_iterator it, const_point_iterator pb, 
      const_weight_iterator wb ) : 
      idx_it(it), point_begin(pb), weight_begin(wb) { }

    const_iterator& operator++(){ idx_it++; return *this; }
    const_iterator  operator++(int) { 
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    const_iterator& operator+(int i) {
      idx_it += i;
      return *this;
    } 

    bool operator==( iterator other ){ return idx_it == other.idx_it; }
    bool operator!=( iterator other ){ return not (*this == other);   }



    auto range() {

      const auto idx      = *idx_it;
      const auto idx_next = *(idx_it+1);

      const auto npts = idx_next - idx;
    
      return std::make_tuple(
        npts,
        point_begin + idx,
        point_begin + idx + npts,
        weight_begin + idx,
        weight_begin + idx + npts 
      );

    }

    value_type operator*() {

      auto [npts,pb,pe,wb,we] = range();
      auto [box_lo, box_up]   = detail::get_box_bounds_points(pb, pe);

      return std::make_tuple(
        box_lo,
        box_up,
        std::vector<point_type>( pb, pe ),
        std::vector<weight_type>( wb, we )
      );

    }

    
  };

public:

  SphericalMicroBatcher( size_t batch_size, std::shared_ptr<quad_type> quad,
    index_container part_idx ) : 
    max_batch_sz_(batch_size), quad_(quad), partition_idx_(part_idx) { }
    

  SphericalMicroBatcher( size_t batch_size, std::shared_ptr<quad_type> quad ) :
    max_batch_sz_(batch_size), quad_(quad) { 

    partition_idx_ = generate_batches();

  }
    
  template <typename Quad>
  using enable_if_spherical_quad = 
    std::enable_if< std::is_base_of_v<quad_type,Quad> >;

  template <typename Quad>
  using enable_if_spherical_quad_t = 
    typename enable_if_spherical_quad<Quad>::type;
  

  template <typename Quad, typename = enable_if_spherical_quad_t<Quad> >
  SphericalMicroBatcher( size_t batch_size, const Quad& quad ):
    SphericalMicroBatcher( 
      batch_size, 
      std::make_shared<Quad>( quad ) 
    ) { }


  template <typename Quad, typename = enable_if_spherical_quad_t<Quad> >
  SphericalMicroBatcher( size_t batch_size, const Quad& quad, index_container idx ):
    SphericalMicroBatcher( 
      batch_size, 
      std::make_shared<Quad>( quad ),
      idx
    ) { }

  template <typename Quad, typename = enable_if_spherical_quad_t<Quad> >
  SphericalMicroBatcher( size_t batch_size, Quad&& quad ):
    SphericalMicroBatcher( 
      batch_size, 
      std::make_shared<Quad>( std::move(quad) ) 
    ) { }


  template <typename Quad, typename = enable_if_spherical_quad_t<Quad> >
  SphericalMicroBatcher( size_t batch_size, Quad&& quad, index_container idx ):
    SphericalMicroBatcher( 
      batch_size, 
      std::make_shared<Quad>( std::move(quad) ),
      idx
    ) { }



  iterator begin() { 
    return iterator( partition_idx_.begin(), quad_->points().begin(),
      quad_->weights().begin() );
  }
  iterator end() { 
    return iterator( partition_idx_.end()-1, quad_->points().begin(),
      quad_->weights().begin() );
  }

  const_iterator cbegin() const { 
    return const_iterator( partition_idx_.cbegin(), quad_->points().cbegin(),
      quad_->weights().cbegin() );
  }
  const_iterator cend() const { 
    return const_iterator( partition_idx_.cend()-1, quad_->points().cbegin(),
      quad_->weights().cbegin() );
  }

  typename iterator::value_type at( size_t i ) {
    if( i >= nbatches() )
      throw std::runtime_error("Index out of bounds");

    return *(begin() + i);
  }

  typename const_iterator::value_type at( size_t i ) const {
    if( i >= nbatches() )
      throw std::runtime_error("Index out of bounds");

    return *(cbegin() + i);
  }


  SphericalMicroBatcher( const SphericalMicroBatcher& ) = default;
  SphericalMicroBatcher( SphericalMicroBatcher&& ) noexcept = default;


  SphericalMicroBatcher clone() const {
    return SphericalMicroBatcher( max_batch_sz_, quad_->clone(), partition_idx_ );
  }



  const quad_type& quadrature() const { return *quad_; }
  quad_type&       quadrature()       { return *quad_; }

  const point_container&  points()  const { return quad_->points();  }
  point_container&        points()        { return quad_->points();  }
  const weight_container& weights() const { return quad_->weights(); }
  weight_container&       weights()       { return quad_->weights(); }




  size_t npts()     const { return quad_->npts();           }
  size_t nbatches() const { return partition_idx_.size()-1; }
  size_t max_batch_size() const { return max_batch_sz_;     }


  inline void recenter( point_type new_center ) {
    quad_->recenter(new_center);
  }

};
















template <typename GridPartitioner, typename QuadType> 
SphericalMicroBatcher< 
  typename QuadType::point_container, 
  typename QuadType::weight_container,
  GridPartitioner
> make_batcher( size_t batch_size, const QuadType& quad ) {

  using point_container = typename QuadType::point_container;
  using weight_container = typename QuadType::weight_container;

  return SphericalMicroBatcher<point_container,weight_container,GridPartitioner>(
    batch_size, quad
  );

}

template <typename GridPartitioner, typename QuadType> 
SphericalMicroBatcher< 
  typename QuadType::point_container, 
  typename QuadType::weight_container,
  GridPartitioner
> make_batcher( size_t batch_size, QuadType&& quad ) {

  using point_container = typename QuadType::point_container;
  using weight_container = typename QuadType::weight_container;

  return SphericalMicroBatcher<point_container,weight_container,GridPartitioner>(
    batch_size, std::move(quad)
  );

}

template <typename GridPartitioner, typename QuadType> 
SphericalMicroBatcher< 
  typename QuadType::point_container, 
  typename QuadType::weight_container,
  GridPartitioner
> make_batcher( size_t batch_size, std::shared_ptr<QuadType> quad ) {

  using point_container = typename QuadType::point_container;
  using weight_container = typename QuadType::weight_container;

  return SphericalMicroBatcher<point_container,weight_container,GridPartitioner>(
    batch_size, quad
  );

}



}
