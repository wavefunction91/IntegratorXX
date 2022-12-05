#pragma once
#include <integratorxx/batch/batching_util.hpp>
#include <chrono>
#include <iomanip>
//#define INTEGRATORXX_VERBOSITY_LEVEL 1

namespace IntegratorXX {

namespace detail {

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
  for( auto i = 0ul; i < npart; ++i ) {
    x_part[i] = bbox_lo[0] + i * delta_x;
    y_part[i] = bbox_lo[1] + i * delta_y;
    z_part[i] = bbox_lo[2] + i * delta_z;
  }
  x_part.back() = bbox_up[0];
  y_part.back() = bbox_up[1];
  z_part.back() = bbox_up[2];

  std::vector<pw_iterator<T>> partition_its = { pw_batch_begin };

  for( auto i = 0ul; i < npart; ++i )
  for( auto j = 0ul; j < npart; ++j )
  for( auto k = 0ul; k < npart; ++k ) {
    std::array<T,3> box_lo = {
      x_part[i], y_part[j], z_part[k]
    };

    std::array<T,3> box_up = {
      x_part[i+1], y_part[j+1], z_part[k+1]
    };

    auto point_in_this_box = [&]( const auto& pw ) {
      return detail::point_in_box( box_lo, box_up, pw.first );
    };
#if 0
    auto point_on_boundary = [&]( const auto& pw ) {
      auto min_x = std::min( 
        std::abs(pw.first[0] - box_lo[0]), 
        std::abs(pw.first[0] - box_up[0])
      );
      auto min_y = std::min( 
        std::abs(pw.first[1] - box_lo[2]), 
        std::abs(pw.first[1] - box_up[2])
      );
      auto min_z = std::min( 
        std::abs(pw.first[2] - box_lo[2]), 
        std::abs(pw.first[2] - box_up[2])
      );

      return min_x < 1e-18 or min_y < 1e-18 or min_z < 1e-18;
    };
    if( std::any_of(partition_its.back(), pw_batch_end, point_on_boundary) ) {
      std::cout << "POINT ON BOUNDARY" << std::endl;
    }
#endif

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
} // namespace detail


struct OctreeGridPartitioner {

  using index_type      = size_t;
  using index_container = std::vector<index_type>;

  template <typename QuadType>
  static index_container generate_batches( QuadType& quad, size_t max_sz ) {

    using point_container  = typename QuadType::point_container;
    using weight_container = typename QuadType::weight_container;
    using point_type  = typename point_container::value_type;
    using weight_type = typename weight_container::value_type;


    using quad_aos_type = std::pair<point_type,weight_type>;

    // Unroll quadrature into AoS
    const size_t npts = quad.npts();
    std::vector<quad_aos_type> q(npts);
    for( size_t i = 0; i < npts; ++i ) {
      q[i] = { quad.points()[i], quad.weights()[i] };
    }

    // Start the clock
    auto part_st = std::chrono::high_resolution_clock::now();


    // Perform initial box partition into 27 cubes
    auto partition_its = detail::partition_box( 3, q.begin(), q.end() );
    auto nparts = partition_its.size()-1;
    std::vector<size_t> part_counts( nparts );
    for( auto ip = 0ul; ip < nparts; ++ip ) {
      part_counts[ip] = std::distance(partition_its[ip], partition_its[ip+1]);
    }

    // Refine spatial partition until all partitions contain less
    // than max_sz points
    auto check_max_batch = [=]( size_t bsz ){ return bsz > max_sz; };

    size_t refine_iter = 0;
    while( std::any_of( part_counts.begin(), part_counts.end(), 
           check_max_batch ) ) {  

  #if INTEGRATORXX_VERBOSITY_LEVEL > 0
      std::cout << "OCTREE REFINE ITER = " << std::setw(4) << refine_iter 
                << " NBATCH = " << std::setw(6) << nparts 
                << " MAX BATCH  = " << 
        *std::max_element(part_counts.begin(), part_counts.end() ) << std::endl;
  #endif
      refine_iter++;

      // Partition large boxes into smaller boxes
      decltype( partition_its ) new_part_its;

      for( auto ip = 0ul; ip < nparts; ++ip ) 
      if( check_max_batch( part_counts[ip] ) ) {

        auto batch_st = partition_its[ip];
        auto batch_en = partition_its[ip+1];

        // Determine tight bounding box
        auto [bbox_lo, bbox_up] = detail::get_box_bounds( batch_st, batch_en );
 
        // If the box contains the origin, parition into 27 smalled boxes,
        // else partition into 8 (to preserve symmetry around origin)
        auto npart_batch = 
          detail::point_in_box( bbox_lo, bbox_up, {0., 0., 0.} ) ? 3 : 2;

        // Partition points into each of the new boxes
        auto batch_part_its = detail::partition_box( npart_batch, 
          bbox_lo, bbox_up, batch_st, batch_en );

        // Insert the new partition points
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
      for( auto ip = 0ul; ip < nparts; ++ip ) {
        part_counts[ip] = std::distance(partition_its[ip], partition_its[ip+1]);
      }

    } // Refinement loop

  #if INTEGRATORXX_VERBOSITY_LEVEL > 0
    std::cout << std::endl;
    std::cout << "OCTREE FINAL NBATCH = " <<  nparts << 
      " FINAL MAX BATCH = " << 
      *std::max_element(part_counts.begin(), part_counts.end() ) << std::endl;
  #endif

    // Stop the clock
    auto part_en = std::chrono::high_resolution_clock::now();
  
    double part_dur = 
      std::chrono::duration<double,std::milli>( part_en-part_st ).count();
  #if INTEGRATORXX_VERBOSITY_LEVEL > 0
    std::cout << "Octree Part Dur = " << part_dur << std::endl;
  #else
    (void)(part_dur); // Unused pascification
  #endif

    // Unroll back into quadrature SoA
    for( size_t i = 0; i < npts; ++i )
      std::tie( quad.points()[i], quad.weights()[i] ) = q[i];


    // Compute the partition points
    auto npart = partition_its.size();
    std::vector<size_t> partition_idx( partition_its.size() );
    for( size_t i = 0; i < npart; ++i )
      partition_idx[i] = std::distance( q.begin(), partition_its[i] );
    
    return partition_idx;
  
  }
};

}
