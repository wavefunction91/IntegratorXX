#pragma once

#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <chrono>
#include <iostream>
#include <iomanip>


namespace IntegratorXX {

namespace detail {

template <typename T>
bool point_in_box( 
  const cartesian_pt_t<T>& lo, 
  const cartesian_pt_t<T>& up,
  const cartesian_pt_t<T>& pt
) {

  if( pt[0] > up[0] or pt[0] < lo[0] )      return false;
  else if( pt[1] > up[1] or pt[1] < lo[1] ) return false;
  else if( pt[2] > up[2] or pt[2] < lo[2] ) return false;
  else                                      return true;

}


template <size_t I, typename T>
inline constexpr T get_coordinate( const cartesian_pt_t<T>& pt ) {
  return pt[I];
}

template <size_t I, typename T, typename BinaryPred>
inline constexpr bool comp_coordinate( 
  const cartesian_pt_t<T>& a, 
  const cartesian_pt_t<T>& b, 
  const BinaryPred& pred 
) {

  return pred( get_coordinate<I>(a), get_coordinate<I>(b) );

}


template <typename T>
using pw_type = std::pair< cartesian_pt_t<T>, T >;

template <typename T>
using pw_container = typename std::vector< pw_type<T> >;

template <typename T>
using pw_iterator  = typename pw_container<T>::iterator;


template <typename PWIterator>
auto get_box_bounds(
  PWIterator  pw_batch_begin,
  PWIterator  pw_batch_end
) {

  using traits = std::iterator_traits<PWIterator>;
  using cart_type = typename traits::value_type::first_type::value_type;
  using less = std::less<cart_type>;

  auto comp_x = [](const auto& a, const auto& b) {
    return comp_coordinate<0>( a.first, b.first, less());
  };
  auto comp_y = [](const auto& a, const auto& b) {
    return comp_coordinate<1>( a.first, b.first, less());
  };
  auto comp_z = [](const auto& a, const auto& b) {
    return comp_coordinate<2>( a.first, b.first, less());
  };

  auto [min_x_it, max_x_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_x );
  auto [min_y_it, max_y_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_y );
  auto [min_z_it, max_z_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_z );

  const auto min_x = min_x_it->first[0];
  const auto max_x = max_x_it->first[0];
  const auto min_y = min_y_it->first[1];
  const auto max_y = max_y_it->first[1];
  const auto min_z = min_z_it->first[2];
  const auto max_z = max_z_it->first[2];

  return std::tuple(
    cartesian_pt_t<cart_type>{ min_x, min_y, min_z },
    cartesian_pt_t<cart_type>{ max_x, max_y, max_z }
  );
}

template <typename PIterator>
auto get_box_bounds_points(
  PIterator  p_batch_begin,
  PIterator  p_batch_end
) {

  using traits = std::iterator_traits<PIterator>;
  using cart_type = typename traits::value_type::value_type;
  using less = std::less<cart_type>;

  auto comp_x = [](const auto& a, const auto& b) {
    return comp_coordinate<0>( a, b, less());
  };
  auto comp_y = [](const auto& a, const auto& b) {
    return comp_coordinate<1>( a, b, less());
  };
  auto comp_z = [](const auto& a, const auto& b) {
    return comp_coordinate<2>( a, b, less());
  };

  auto [min_x_it, max_x_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_x );
  auto [min_y_it, max_y_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_y );
  auto [min_z_it, max_z_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_z );

  const auto min_x = (*min_x_it)[0];
  const auto max_x = (*max_x_it)[0];
  const auto min_y = (*min_y_it)[1];
  const auto max_y = (*max_y_it)[1];
  const auto min_z = (*min_z_it)[2];
  const auto max_z = (*max_z_it)[2];

  return std::tuple(
    cartesian_pt_t<cart_type>{ min_x, min_y, min_z },
    cartesian_pt_t<cart_type>{ max_x, max_y, max_z }
  );
}





template <typename T>
auto partition_box( 
  size_t                   npart, 
  const cartesian_pt_t<T>& bbox_lo,
  const cartesian_pt_t<T>& bbox_up,
  pw_iterator<T>           pw_batch_begin,
  pw_iterator<T>           pw_batch_end
) {

  const auto extent_x = bbox_up[0] - bbox_lo[0];
  const auto extent_y = bbox_up[1] - bbox_lo[1];
  const auto extent_z = bbox_up[2] - bbox_lo[2];

  std::vector<T> x_part(npart+1), y_part(npart+1), z_part(npart+1);

  const auto delta_x = extent_x / npart;
  const auto delta_y = extent_y / npart;
  const auto delta_z = extent_z / npart;
  for( auto i = 0; i < npart; ++i ) {
    x_part[i] = bbox_lo[0] + i * delta_x;
    y_part[i] = bbox_lo[1] + i * delta_y;
    z_part[i] = bbox_lo[2] + i * delta_z;
  }
  x_part.back() = bbox_up[0];
  y_part.back() = bbox_up[1];
  z_part.back() = bbox_up[2];

  std::vector partition_its = { pw_batch_begin };

  for( int i = 0; i < npart; ++i )
  for( int j = 0; j < npart; ++j )
  for( int k = 0; k < npart; ++k ) {
    std::array box_lo = {
      x_part[i], y_part[j], z_part[k]
    };

    std::array box_up = {
      x_part[i+1], y_part[j+1], z_part[k+1]
    };

    auto point_in_this_box = [&]( const auto& pw ) {
      return detail::point_in_box( box_lo, box_up, pw.first );
    };

    partition_its.emplace_back(
      std::partition( partition_its.back(), pw_batch_end,
        point_in_this_box )
    );
  }

  // Remove batches with no points
  auto uniq_end = std::unique( partition_its.begin(), partition_its.end() );
  partition_its.erase( uniq_end, partition_its.end() );
  partition_its.shrink_to_fit();

  return partition_its;

}



template <typename PWIterator>
auto partition_box( 
  size_t          npart, 
  PWIterator  pw_batch_begin,
  PWIterator  pw_batch_end
) {

  auto [bbox_lo, bbox_up] = get_box_bounds( pw_batch_begin, pw_batch_end );
  return partition_box( npart, bbox_lo, bbox_up, pw_batch_begin, pw_batch_end );

}



}

template <typename RadialQuad, typename AngularQuad>
class SphericalMicroBatcher {

  using quad_type   = SphericalQuadrature<RadialQuad,AngularQuad>;
  using point_type  = typename quad_type::point_type;
  using weight_type = typename quad_type::weight_type;
  using point_container  = typename quad_type::point_container;
  using weight_container = typename quad_type::weight_container;

  using index_type      = size_t;
  using index_container = std::vector<index_type>;


  using point_iterator  = typename point_container::iterator;
  using weight_iterator = typename weight_container::iterator;
  using index_iterator  = typename index_container::iterator;

  size_t max_batch_sz_;

  std::shared_ptr<quad_type> quad_;
  index_container            partition_idx_;

  index_container generate_batches();


  struct iterator {
  
    using value_type        = std::tuple<point_type,point_type,point_container,weight_container>;
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
    
      return std::tuple(
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

      return std::tuple(
        box_lo,
        box_up,
        std::vector( pb, pe ),
        std::vector( wb, we )
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
    
  SphericalMicroBatcher( size_t batch_size, const quad_type& quad ):
    SphericalMicroBatcher( 
      batch_size, std::make_shared<quad_type>( quad ) 
    ) { }
  SphericalMicroBatcher( size_t batch_size, const quad_type& quad, 
    index_container idx ):
    SphericalMicroBatcher( 
      batch_size, std::make_shared<quad_type>( quad ), idx 
    ) { }

  SphericalMicroBatcher( size_t batch_size, quad_type&& quad ):
    SphericalMicroBatcher( 
      batch_size, std::make_shared<quad_type>( std::move(quad) ) 
    ) { }
  SphericalMicroBatcher( size_t batch_size, quad_type&& quad, 
    index_container idx ):
    SphericalMicroBatcher( 
      batch_size, std::make_shared<quad_type>( std::move(quad) ), idx 
    ) { }


  iterator begin() { 
    return iterator( partition_idx_.begin(), quad_->points().begin(),
      quad_->weights().begin() );
  }
  iterator end() { 
    return iterator( partition_idx_.end()-1, quad_->points().begin(),
      quad_->weights().begin() );
  }

  typename iterator::value_type at( size_t i ) {
    if( i >= nbatches() )
      throw std::runtime_error("Index out of bounds");

    return *(begin() + i);
  }


  SphericalMicroBatcher( const SphericalMicroBatcher& ) = default;
  SphericalMicroBatcher( SphericalMicroBatcher&& ) noexcept = default;


  SphericalMicroBatcher clone() const {
    return SphericalMicroBatcher( max_batch_sz_, *quad_, partition_idx_ );
  }



  const quad_type& quadrature() const { return *quad_; }
  quad_type&       quadrature()       { return *quad_; }

  const point_container&  points()  const { return quad_->points();  }
  point_container&        points()        { return quad_->points();  }
  const weight_container& weights() const { return quad_->weights(); }
  weight_container&       weights()       { return quad_->weights(); }




  size_t npts()     const { return quad_->npts();           }
  size_t nbatches() const { return partition_idx_.size()-1; }

};

template <typename RadialQuad, typename AngularQuad>
std::vector<size_t> SphericalMicroBatcher<RadialQuad,AngularQuad>::generate_batches() {



  std::vector<std::pair<point_type, weight_type>> q( quad_->npts() );
  for( size_t i = 0; i < quad_->npts(); ++i ) {
    q[i] = { quad_->points()[i], quad_->weights()[i] };
  }

  auto part_st = std::chrono::high_resolution_clock::now();



  auto partition_its = detail::partition_box( 3, q.begin(), q.end() );
  auto nparts = partition_its.size()-1;
  std::vector<size_t> part_counts( nparts );
  for( auto ip = 0; ip < nparts; ++ip ) {
    part_counts[ip] = std::distance(partition_its[ip], partition_its[ip+1]);
  }

  
  size_t refine_iter = 0;
  auto check_max_batch = [&]( size_t bsz ){ return bsz > max_batch_sz_; };
  while( std::any_of( part_counts.begin(), part_counts.end(), 
         check_max_batch ) ) {

    std::cout << "REFINE ITER = " << std::setw(4) << refine_iter++ 
              << " NBATCH = " << std::setw(6) << nparts 
              << " MAX BATCH  = " << *std::max_element(part_counts.begin(), part_counts.end() ) << std::endl;


    // Partition large boxes
    decltype( partition_its ) new_part_its;

    for( auto ip = 0; ip < nparts; ++ip ) 
    if( check_max_batch( part_counts[ip] ) ) {

      auto batch_st = partition_its[ip];
      auto batch_en = partition_its[ip+1];

      auto [bbox_lo, bbox_up] = detail::get_box_bounds( batch_st, batch_en );
 
      auto npart_batch = 
        detail::point_in_box( bbox_lo, bbox_up, {0., 0., 0.} ) ? 3 : 2;

      auto batch_part_its = 
        detail::partition_box( npart_batch, bbox_lo, bbox_up, batch_st, batch_en );

      new_part_its.insert( new_part_its.end(), batch_part_its.begin(), 
        batch_part_its.end() );
    } 

    // Sort and remove duplicate iterators of newly added partitions
    std::sort( new_part_its.begin(), new_part_its.end() );
    auto uniq_new_part_end = 
      std::unique( new_part_its.begin(), new_part_its.end() );
    new_part_its.erase( uniq_new_part_end, new_part_its.end() );

    // Take union of new and old (XXX Should already by sorted)
    decltype( partition_its ) union_part_its;
    std::set_union( 
      partition_its.begin(), partition_its.end(),
      new_part_its.begin(),  new_part_its.end(),
      std::back_inserter( union_part_its )
    );
    
    // Compute new partition sizes
    partition_its = std::move( union_part_its );
    nparts = partition_its.size()-1;
    part_counts.resize( nparts );
    for( auto ip = 0; ip < nparts; ++ip ) {
      part_counts[ip] = std::distance(partition_its[ip], partition_its[ip+1]);
    }

    

  }

  std::cout << std::endl;
  std::cout << "FINAL NBATCH = " <<  nparts << "     FINAL MAX BATCH = " << *std::max_element(part_counts.begin(), part_counts.end() ) << std::endl;











  auto part_en = std::chrono::high_resolution_clock::now();

  double part_dur = std::chrono::duration<double,std::milli>( part_en-part_st ).count();
  std::cout << "Part Dur = " << part_dur << std::endl;
  
  for( size_t i = 0; i < quad_->npts(); ++i )
    std::tie( quad_->points()[i], quad_->weights()[i] ) = q[i];


  auto npart = partition_its.size();
  std::vector<size_t> partition_idx( partition_its.size() );
  for( size_t i = 0; i < npart; ++i )
    partition_idx[i] = std::distance( q.begin(), partition_its[i] );
  
  return partition_idx;

}

}
