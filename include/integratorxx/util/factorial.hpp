#pragma once

namespace IntegratorXX {

namespace detail {

constexpr static std::array<size_t, 20> factorials = {
  1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
  39916800, 479001600, 6227020800, 87178291200, 1307674368000,
  20922789888000, 355687428096000, 6402373705728000, 121645100408832000
};

constexpr static std::array<size_t,20> double_factorials = {
  0, 1, 1, 3, 6, 15, 38, 105, 306, 945, 3063, 10395,
  36766, 135135, 514731, 2027025, 8235700, 34459425,
  148242610, 654729075
};



}

template <typename T>
T factorial(unsigned n) {
  if( n < 20 ) {
    return detail::factorials[n];
  } else {
    return n * factorial<T>(n-1);
  }
}

template <typename T>
T double_factorial(int n) {
  if(n == -3)      return -1;
  else if(n == -1) return  1;
  else if( n < 20 ) {
    return detail::double_factorials[size_t(n)];
  } else {
    return n * double_factorial<T>(n-2);
  }
}

inline double pow_2(int n) {
  return 1ul << n;
}

// GAMMA(n/2)
template <typename T>
T half_integral_tgamma(int n) {
  assert(n%2);
  if(n == 1) return std::sqrt(M_PI);
  if(n == 3) return 0.5 * std::sqrt(M_PI);
  return double_factorial<T>(n-2) / pow_2(n-3) / M_2_SQRTPI;
}

template <typename T>
T integral_tgamma(T n) {
  std::cout << n << std::endl;
  return factorial<T>(n-1);
}


}
