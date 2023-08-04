#include "catch2/catch_all.hpp"
#include <integratorxx/quadratures/gausslegendre.hpp>
#include <integratorxx/quadratures/gausslobatto.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/treutlerahlrichs.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/gausschebyshev1.hpp>
#include <integratorxx/quadratures/gausschebyshev2.hpp>
#include <integratorxx/quadratures/gausschebyshev2modified.hpp>
#include <integratorxx/quadratures/gausschebyshev3.hpp>
#include <integratorxx/quadratures/ahrens_beylkin.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>

#include "quad_matcher.hpp"
#include "test_functions.hpp"

using namespace IntegratorXX;

inline constexpr double inf = std::numeric_limits<double>::infinity();
inline constexpr double eps = std::numeric_limits<double>::epsilon();

template <typename T = double>
constexpr T gaussian( T alpha, T c, T x ) {
  return std::exp( -alpha * (x-c) * (x-c) );
}

template <typename T = double>
constexpr T gaussian( T x ) {
  return gaussian( 1., 0., x );
}

template <typename T = double>
constexpr T ref_gaussian_int( T alpha, T c, T a, T b ) {

  const auto low = (a == inf) ? T(-1.) : std::erf( std::sqrt(alpha) * (c - a) );
  const auto hgh = (b == inf) ? T(-1.) : std::erf( std::sqrt(alpha) * (c - b) );

  return 0.5 * std::sqrt( M_PI / alpha ) * (low - hgh);
}

template <typename T = double>
constexpr T ref_gaussian_int( T a, T b ) {

  return ref_gaussian_int( 1., 0., a, b );

}


template <typename T>
auto chebyshev_T(int n, T x) {
  return std::cos( n * std::acos(x) );
}






TEST_CASE( "Gauss-Legendre Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    std::vector<double> cp(p+1, 0.0); 
    for(int i = 0; i < p; ++i) {
      cp[i] = c[i] / (p-i);
    }
    return Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussLegendre<double,double>;
  test_random_polynomial<quad_type, Polynomial>("Gauss-Legendre", 10, 14, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Lobatto Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    std::vector<double> cp(p+1, 0.0); 
    for(int i = 0; i < p; ++i) {
      cp[i] = c[i] / (p-i);
    }
    return Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussLobatto<double,double>;
  test_random_polynomial<quad_type, Polynomial>("Gauss-Lobatto", 10, 14, 
    [](int o){ return 2*o-1; }, // Max order 2N-3
    ref_value, 1e-12 );


}

TEST_CASE( "Gauss-Chebyshev T1 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T1 evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
    for(int i = 0; i < p; ++i) {
      int k = p - i - 1;
      if(k > 0)
        ref += c[i] * (std::sqrt(M_PI) / k) * (std::pow(-1,k)+1) * 
               std::tgamma((k+1)/2.0) / std::tgamma(k/2.0);
      else
        ref += c[i] * M_PI;
    }
    return ref;
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussChebyshev1<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT1WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T1)", 10, 100, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Chebyshev T2 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T2 evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
    for(int i = 0; i < p; ++i) {
      int k = p - i - 1;
      if(k > 0)
        ref += c[i] * (std::sqrt(M_PI) / 4.0) * (std::pow(-1,k)+1) * 
               std::tgamma((k+1)/2.0) / std::tgamma(k/2.0 + 2);
      else
        ref += c[i] * M_PI/2.0;
    }
    return ref;
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussChebyshev2<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT2WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T2)", 10, 100, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Chebyshev T3 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T3 evaluated over [0,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
      for(int i = 0; i < p; ++i) {
        int k = p - i - 1;
        ref += c[i] * std::sqrt(M_PI) * std::tgamma(k+1.5) / std::tgamma(k+2); 
      }
    return ref;
  };

  // Test Quadrature for Correctness
  // TODO: Code breaks down for large orders here
  using quad_type = IntegratorXX::GaussChebyshev3<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT3WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T2)", 10, 50, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}


TEST_CASE( "Euler-Maclaurin Quadratures", "[1d-quad]" ) {
  IntegratorXX::MurrayHandyLaming<double,double> quad(150);
  const auto msg = "Euler-Maclaurin N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Treutler-Ahlrichs Quadratures", "[1d-quad]" ) {
  IntegratorXX::TreutlerAhlrichs<double,double> quad(150);
  const auto msg = "Treutler-Ahlrichs N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Knowles Quadratures", "[1d-quad]" ) {
  IntegratorXX::MuraKnowles<double,double> quad(350);
  const auto msg = "Mura-Knowles N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Lebedev-Laikov", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::LebedevLaikov<double> quad( nPts );
    for( auto l = 1; l < 10; ++l )
    for( auto m = 0; m <= l; ++m ) {
      const auto msg = "Lebedev-Laikov N = " + std::to_string(quad.npts());
      test_quadrature<MagnitudeSquaredSphericalHarmonic>(msg, quad, 1.0, 1e-10, l, m); 
    }

  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}

TEST_CASE( "Ahrens-Beylkin", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::AhrensBeylkin<double> quad( nPts );

    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    for( auto l = 1; l < 10; ++l )
    for( auto m = 0; m <= l; ++m ) {
      const auto msg = "Ahrens-Beylkin N = " + std::to_string(quad.npts());
      test_quadrature<MagnitudeSquaredSphericalHarmonic>(msg, quad, 1.0, 1e-10, l, m); 
    }

  };

  test_fn(312);
  test_fn(792);
  test_fn(972);
}
