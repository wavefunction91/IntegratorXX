#pragma once
#include <cmath>
#include <complex>
#include <random>
#include <version>
#include "quad_matcher.hpp"

// The C++17 mathematical special functions ([sf.cmath]) are an optional part of
// the standard library. libstdc++ and the MSVC STL provide them; libc++ does
// not, so std::assoc_legendre is unavailable on macOS and on any clang build
// against libc++. Detect it with the standard feature-test macro rather than by
// sniffing the standard library, so this keeps working if libc++ implements
// [sf.cmath] later.
#if defined(__cpp_lib_math_special_functions) && \
    __cpp_lib_math_special_functions >= 201603L
#  define IXX_TEST_HAS_STD_ASSOC_LEGENDRE 1
#endif

namespace detail {
  constexpr size_t factorial(size_t n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }
}
struct AssociatedLegendre {
  /**
   *  Associated Legendre polynomial P_l^m(x), evaluated with the standard
   *  upward recurrence. This is compiled on every platform, whether or not
   *  std::assoc_legendre exists, so that the parity test in assoc_legendre.cxx
   *  can check it against the standard library wherever that is available.
   *
   *  Follows the [sf.cmath] convention, P_l^m(x) = (1-x^2)^(m/2) d^m/dx^m
   *  P_l(x), i.e. WITHOUT the Condon-Shortley phase (-1)^m. This is the same
   *  convention std::assoc_legendre uses, which is what makes the two
   *  interchangeable.
   */
  static inline double evaluate_fallback(int l, int m, double x) {
    // Preconditions of std::assoc_legendre; callers pass std::abs(m) and a cos
    if( m < 0 or m > l or std::abs(x) > 1.0 ) return 0.0;

    // P_m^m(x) = (2m-1)!! (1-x^2)^(m/2), accumulated as a product of the exact
    // integer coefficients (2j-1) and sqrt(1-x^2) to avoid forming (2m-1)!! and
    // the fractional power separately
    const double s = std::sqrt( (1.0 - x) * (1.0 + x) );
    double p_mm = 1.0;
    for( int j = 1; j <= m; ++j ) p_mm *= (2 * j - 1) * s;
    if( l == m ) return p_mm;

    // P_{m+1}^m(x) = x (2m+1) P_m^m(x)
    double p_lm_m1 = p_mm;
    double p_lm    = x * (2 * m + 1) * p_mm;

    // (l-m) P_l^m(x) = x (2l-1) P_{l-1}^m(x) - (l+m-1) P_{l-2}^m(x)
    for( int ll = m + 2; ll <= l; ++ll ) {
      const double p_lm_m2 = p_lm_m1;
      p_lm_m1 = p_lm;
      p_lm = ( x * (2 * ll - 1) * p_lm_m1 - (ll + m - 1) * p_lm_m2 ) / (ll - m);
    }

    return p_lm;
  }

  static inline double evaluate(int l, int m, double x) {
#ifdef IXX_TEST_HAS_STD_ASSOC_LEGENDRE
    return std::assoc_legendre(l, m, x);
#else
    return evaluate_fallback(l, m, x);
#endif
  }
};

struct SphericalHarmonic {
  static auto evaluate(int l, int m, double theta, double phi) {
    double fac_ratio  = static_cast<double>(detail::factorial(l-m));
           fac_ratio /= detail::factorial(l+m);
    double prefactor  = (2.0 * l + 1.0) / (4.0 * M_PI);
           prefactor *= fac_ratio;
    prefactor = std::sqrt(prefactor);
    if( m < 0 ) {
      prefactor *= std::pow(1, std::abs(m)) * fac_ratio;
    }

    return prefactor * 
           AssociatedLegendre::evaluate(l, std::abs(m), std::cos(theta)) *
           std::complex<double>( std::cos(m*phi), std::sin(m*phi) );
  }

  static auto evaluate(int l, int m, double x, double y, double z) {
    const auto r     = std::sqrt( x*x + y*y + z*z );
    const auto theta = std::acos( z / r );
    const auto phi   = std::atan2( y, x );

    return evaluate(l, m, theta, phi);
  }

  static auto evaluate(int l, int m, std::array<double,3> x) {
    return evaluate(l, m, x[0], x[1], x[2]);
  }
};

struct MagnitudeSquaredSphericalHarmonic {
  template <typename... Args>
  static double evaluate( Args&&... args ) {
    auto eval = std::abs(SphericalHarmonic::evaluate(args...));
    return eval * eval;
  }
};

struct RadialGaussian {
  static double evaluate(double r) {
    return r*r * std::exp(-r*r);
  }
  static double evaluate(double x, double y, double z) {
    const auto r = std::hypot(x,y,z);
    return evaluate(r);
  }
  static double evaluate(std::array<double,3> x) {
    return evaluate(x[0], x[1], x[2]);
  }
};



struct Polynomial {

  static double evaluate( std::vector<double> c, double x ) {
    const auto N = c.size();
    double res = c[0];
    for(int i = 1; i < N; ++i) {
      res = x * res + c[i];
    }
    return res;
  }

};

template <typename WeightFunctor>
struct WeightedPolynomial {

  static double evaluate( std::vector<double> c, double x ) {
    return Polynomial::evaluate(c,x) * WeightFunctor::evaluate(x);
  }

};

struct ChebyshevT1WeightFunction {
  static double evaluate(double x){ 
    return 1./std::sqrt(1.0 - x*x);
  }
};
struct ChebyshevT2WeightFunction {
  static double evaluate(double x){ 
    return std::sqrt(1.0 - x*x);
  }
};
struct ChebyshevT3WeightFunction {
  static double evaluate(double x){ 
    return std::sqrt(x/(1.0 - x));
  }
};


template <typename TestFunction, typename QuadType, typename... PreArgs>
void test_quadrature(std::string msg, const QuadType& quad, double ref, double e, PreArgs&&... args) {
  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  double res = 0.0;
  for(auto i = 0; i < quad.npts(); ++i) {
    res += wgt[i] * TestFunction::evaluate(args..., pts[i]);
  }
  
  //standard_matcher(mes, res, ref, e);
  //printf("diff = %.6e\n", std::abs(ref - res));
  REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, e));
}

template <typename QuadType>
void test_angular_quadrature(std::string msg, const QuadType& quad, int maxL, double e) {

  for( auto l = 1; l < maxL; ++l )
  for( auto m = 0; m <= l; ++m ) {
    auto loc_msg = msg + "(L,M) = (" + std::to_string(l) + "," + std::to_string(m) + ")";
    test_quadrature<MagnitudeSquaredSphericalHarmonic>(loc_msg, quad, 1.0, e, l, m); 
  }

}


template <typename QuadType, typename TestFunction>
void test_random_polynomial(std::string qname, int min_order, int max_order,
  std::function<int(int)> max_poly_order_functor,
  std::function<double(const std::vector<double>&)> ref_functor,
  double e) {

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(int order = min_order; order < max_order; order++) {

    QuadType quad(order);
    REQUIRE(std::is_sorted(quad.points().begin(), quad.points().end()));

    const auto p_max = max_poly_order_functor(order);
    REQUIRE(p_max > 2);
    for(int p = 2; p < p_max; ++p) {
      // Generate a random polynomial
      std::vector<double> c(p); 
      std::generate(c.begin(), c.end(), rand_gen);

      // Evaluate reference value
      const auto ref = ref_functor(c);

      // Test Quadrature Value
      const std::string msg = qname + 
        " Order = " + std::to_string(order) +
        " PolyOrder = " + std::to_string(p-1);
      test_quadrature<TestFunction>(msg, quad, ref, e, c);
    }

  }

}
