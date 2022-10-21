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
  for( i = 1; i < N; i++ ) X[i] ^= X[i-1]; t = 0;
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


}
