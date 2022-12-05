#pragma once
#include <climits>

namespace IntegratorXX {

template <size_t B, size_t N, typename WordType>
void inplace_hilbert_encode_transpose( WordType* X ) {
  WordType M=WordType(1)<<(B-1),P,Q,t;
  int i; // Inverse undo
  for( Q = M; Q > 1; Q >>= 1 ) { 
    P = Q - 1;
    for( i = 0; i < N; i++ ) {
      if( X[i] & Q ) X[0] ^= P;
      else { 
        t = (X[0]^X[i]) & P; X[0] ^= t; X[i] ^= t; 
      } 
    } // exchange
  }
  // Gray encode
  for( i = 1; i < N; i++ ) X[i] ^= X[i-1]; 

  t = 0;
  for( Q = M; Q > 1; Q >>= 1 ) {
    if( X[N-1] & Q ) t ^= Q-1; 
  }
  for( i = 0; i < N; i++ ) X[i] ^= t;
}

template <typename ResultType, size_t B, size_t N, typename WordType>
ResultType transpose_to_binary(const WordType* X) {
  static_assert(B   <= sizeof(WordType)*CHAR_BIT,   "WordType Insufficient for Encoding"  );
  static_assert(B*N <= sizeof(ResultType)*CHAR_BIT, "ResultTYpe Insufficient for Encoding");

  ResultType bin = 0;
  WordType mask = WordType(1) << (B-1);
  int shift = (N-1)*B;
  while(mask) {
    for(int i = 0; i < N; ++i) {
      bin |= ResultType(X[i] & mask) << (shift - i);
    }
    mask >>= 1;
    shift -= N-1;
  }

  return bin;
}

template <typename ResultType, size_t B, size_t N, typename WordType>
ResultType transpose_to_binary(const std::array<WordType,N>& X) {
  return transpose_to_binary<ResultType,B,N,WordType>(X.data());
}

template <typename ResultType, size_t B, size_t N, typename WordType>
ResultType hilbert_encode(const WordType* _X) {
  std::array<WordType,N> X; std::copy(_X,_X+N,X.data());
  inplace_hilbert_encode_transpose<B,N>(X.data());
  return transpose_to_binary<ResultType,B>(X);
}
template <typename ResultType, size_t B, size_t N, typename WordType>
ResultType hilbert_encode(const std::array<WordType,N>& _X) {
  std::array<WordType,N> X = _X;
  inplace_hilbert_encode_transpose<B,N>(X.data());
  return transpose_to_binary<ResultType,B>(X);
}



template <typename WordType, size_t B, size_t N, typename CoordContainer>
std::vector<uint64_t> points_to_hilbert(
  const std::array<double,3>& box_lo, double box_dim,
  const CoordContainer& points) {

  const size_t npts = points.size();
  std::vector<uint64_t> hilbert_indices(npts);

  const double spacing = box_dim/double((WordType(1) << B) - 1);
  for(size_t i = 0; i < npts; ++i) {
    const auto& pt = points[i];
    const std::array<WordType,3> coords = {
      std::round((pt[0] - box_lo[0]) / spacing),
      std::round((pt[1] - box_lo[1]) / spacing),
      std::round((pt[2] - box_lo[2]) / spacing)
    };
    hilbert_indices[i] = hilbert_encode<uint64_t,B>(coords);
  }

  return hilbert_indices;
}



template <typename WordType, size_t B, size_t N, typename CoordContainer>
std::vector<uint64_t> points_to_hilbert(const CoordContainer& points) {

  const size_t npts = points.size();

  auto [box_lo, box_hi] = 
    detail::get_box_bounds_points(points.begin(),points.end());

  const std::array<double,3> box_center = {
    0.5 * (box_hi[0] + box_lo[0]),
    0.5 * (box_hi[1] + box_lo[1]),
    0.5 * (box_hi[2] + box_lo[2])
  };

  const std::array<double,3> box_dims = {
    box_hi[0] - box_lo[0],
    box_hi[1] - box_lo[1],
    box_hi[2] - box_lo[2]
  };

  const double max_box_dim = *std::max_element(box_dims.begin(),box_dims.end());
  const std::array<double,3> box_scal_factors = {
    max_box_dim / box_dims[0],
    max_box_dim / box_dims[1],
    max_box_dim / box_dims[2]
  };
  
  box_lo[0] = box_scal_factors[0] * (box_lo[0] - box_center[0]) + box_center[0];
  box_lo[1] = box_scal_factors[1] * (box_lo[1] - box_center[1]) + box_center[1];
  box_lo[2] = box_scal_factors[2] * (box_lo[2] - box_center[2]) + box_center[2];
#if 0
  box_hi[0] = box_scal_factors[0] * (box_hi[0] - box_center[0]) + box_center[0];
  box_hi[1] = box_scal_factors[1] * (box_hi[1] - box_center[1]) + box_center[1];
  box_hi[2] = box_scal_factors[2] * (box_hi[2] - box_center[2]) + box_center[2];
  

  const double spacing = max_box_dim/double((WordType(1) << B) - 1);
  for(size_t i = 0; i < npts; ++i) {
    const auto& pt = points[i];
    const std::array<WordType,3> coords = {
      std::round((pt[0] - box_lo[0]) / spacing),
      std::round((pt[1] - box_lo[1]) / spacing),
      std::round((pt[2] - box_lo[2]) / spacing)
    };
    hilbert_indices[i] = hilbert_encode<uint64_t,B>(coords);
  }
  return hilbert_indices;
#else
  return points_to_hilbert<WordType,B,N>(box_lo, max_box_dim, points);
#endif

}










template <size_t NBITS>
struct HilbertGridPartitioner {

  static_assert(3*NBITS < 64,"Too Many Encoding Bits");

  using index_type      = size_t;
  using index_container = std::vector<index_type>;

  template <typename QuadType>
  static index_container generate_batches( QuadType& quad, size_t max_sz ) {

    using point_container  = typename QuadType::point_container;
    using weight_container = typename QuadType::weight_container;
    using point_type  = typename point_container::value_type;
    using weight_type = typename weight_container::value_type;

    // Compute Hilbert indices and sort grid index
    const size_t npts = quad.npts();
    auto hilbert_indices = points_to_hilbert<uint32_t,NBITS,3>(quad.points());

    std::vector<size_t> sort_indices(npts);
    std::iota(sort_indices.begin(), sort_indices.end(), 0);
    std::sort(sort_indices.begin(), sort_indices.end(),
       [&](size_t i, size_t j) { 
         return hilbert_indices[i] < hilbert_indices[j]; 
       }
    );

    std::vector<point_type> sorted_points(npts);
    for(size_t i = 0; i < npts; ++i) {
      sorted_points[i] = quad.points()[sort_indices[i]];
    }
    std::vector<weight_type> sorted_weights(npts);
    for(size_t i = 0; i < npts; ++i) {
      sorted_weights[i] = quad.weights()[sort_indices[i]];
    }

    quad.points()  = std::move(sorted_points);
    quad.weights() = std::move(sorted_weights);

    const size_t npart = (npts/max_sz) + !!(npts % max_sz) + 1;
    std::vector<size_t> partition_idx(npart);
    std::iota(partition_idx.begin(), partition_idx.end(), 0);
    std::transform(partition_idx.begin(),partition_idx.end(),partition_idx.begin(),
      [=](auto x){ return x * max_sz; });
    partition_idx[npart-1] = npts;
    return partition_idx;

  }
};


}
