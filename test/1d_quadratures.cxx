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

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(unsigned order=10;order<14;order++) {
    IntegratorXX::GaussLegendre<double,double> quad( order );
    for(int p = 2; p < 2*order+1; ++p) {
      std::vector<double> c(p); std::generate(c.begin(), c.end(), rand_gen);
      std::vector<double> cp(p+1, 0.0); 
      for(int i = 0; i < p; ++i) {
        cp[i] = c[i] / (p-i);
      }
      const auto ref = 
        Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
      const auto msg = "Gauss-Legendre Order = " + std::to_string(order) + " PolyOrder = " + std::to_string(p-1);
      test_quadrature<Polynomial>(msg, quad, ref, 1e-12, c);
    }
  }

}

TEST_CASE( "Gauss-Lobatto Quadratures", "[1d-quad]" ) {

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(unsigned order=10;order<14;order++) {
    IntegratorXX::GaussLobatto<double,double> quad( order );
    auto f = quad.points().front();
    auto e = quad.points().back();
    REQUIRE_THAT(f, Catch::Matchers::WithinAbs(-1.0, eps));
    REQUIRE_THAT(e, Catch::Matchers::WithinAbs(1.0, eps));

    for(int p = 2; p < 2*order-1; ++p) {
      std::vector<double> c(p); std::generate(c.begin(), c.end(), rand_gen);
      std::vector<double> cp(p+1, 0.0); 
      for(int i = 0; i < p; ++i) {
        cp[i] = c[i] / (p-i);
      }
      const auto ref = 
        Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
      const auto msg = "Gauss-Lobatto Order = " + std::to_string(order) + " PolyOrder = " + std::to_string(p-1);
      test_quadrature<Polynomial>(msg, quad, ref, 1e-12, c);
    }
  }
}

TEST_CASE( "Gauss-Chebyshev T1 Quadratures", "[1d-quad]" ) {

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(unsigned order=10; order<11; order++) {
    IntegratorXX::GaussChebyshev1<double,double> quad( order );
    REQUIRE(std::is_sorted(quad.points().begin(), quad.points().end()));
    for(int p = 2; p < 2*order+1; ++p) {
      std::vector<double> c(p); std::generate(c.begin(), c.end(), rand_gen);

      double ref = 0.0;
      for(int i = 0; i < p; ++i) {
        int k = p - i - 1;
        if(k > 0)
          ref += c[i] * (std::sqrt(M_PI) / k) * std::tgamma((k+1)/2.0) / std::tgamma(k/2.0) * (std::pow(-1,k)+1);
        else
          ref += c[i] * M_PI;
      }

      const auto msg = "Gauss-Chebyshev T1 Order = " + std::to_string(order) + " PolyOrder = " + std::to_string(p-1);
      test_quadrature<WeightedPolynomial<ChebyshevT1WeightFunction>>(msg, quad, ref, 1e-12, c);
    }
  }

}

TEST_CASE( "Gauss-Chebyshev T2 Quadratures", "[1d-quad]" ) {

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(unsigned order=10; order<100; order++) {
    IntegratorXX::GaussChebyshev2<double,double> quad( order );
    REQUIRE(std::is_sorted(quad.points().begin(), quad.points().end()));
    if(order%2) REQUIRE(quad.points()[order/2] == 0.0); // FP equality b/c data should be explicitly set
    for(int p = 2; p < 2*order+1; ++p) {
      std::vector<double> c(p); std::generate(c.begin(), c.end(), rand_gen);

      double ref = 0.0;
      for(int i = 0; i < p; ++i) {
        int k = p - i - 1;
        if(k > 0)
          ref += c[i] * (std::sqrt(M_PI) / 4.0) * std::tgamma((k+1)/2.0) / std::tgamma(k/2.0 + 2) * (std::pow(-1,k)+1);
        else
          ref += c[i] * M_PI/2.0;
      }

      const auto msg = "Gauss-Chebyshev T2 Order = " + std::to_string(order) + " PolyOrder = " + std::to_string(p-1);
      test_quadrature<WeightedPolynomial<ChebyshevT2WeightFunction>>(msg, quad, ref, 1e-12, c);
    }
  }

}

TEST_CASE( "Gauss-Chebyshev T3 Quadratures", "[1d-quad]" ) {

  std::default_random_engine gen;
  std::uniform_real_distribution<> dist(-1.,1.);
  auto rand_gen = [&]{ return dist(gen); };

  for(unsigned order=10; order<50; order++) { // Code breaks down for large orders here
    IntegratorXX::GaussChebyshev3<double,double> quad( order );
    REQUIRE(std::is_sorted(quad.points().begin(), quad.points().end()));
    for(int p = 2; p < 2*order-1; ++p) {
      std::vector<double> c(p); std::generate(c.begin(), c.end(), rand_gen);

      double ref = 0.0;
      for(int i = 0; i < p; ++i) {
        int k = p - i - 1;
        ref += c[i] * std::sqrt(M_PI) * std::tgamma(k+3.0/2.0) / std::tgamma(k + 2); 
      }

      const auto msg = "Gauss-Chebyshev T3 Order = " + std::to_string(order) + " PolyOrder = " + std::to_string(p-1);
      test_quadrature<WeightedPolynomial<ChebyshevT3WeightFunction>>(msg, quad, ref, 1e-12, c);
    }
  } 

}

#if 0
TEST_CASE( "Gauss-Chebyshev Quadratures", "[1d-quad]") {

  auto integrate = [&](auto& quad) {
    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    // Check that nodes are in increasing value
    REQUIRE(std::is_sorted(pts.begin(), pts.end()));

    //auto f = [=]( double x ){ return std::cyl_bessel_j(6,std::abs(x)); };
    auto f = [=]( double x ){ 
      return chebyshev_T(6,x) * chebyshev_T(7,x)  / std::sqrt(1 - x*x); 
    };

    double res = 0.;
    for( auto i = 0; i < quad.npts(); ++i ) {
      res += wgt[i] * f(pts[i]);
    }

    const auto ref = 0;
    const auto msg = "Gauss-Chebyshev N = " + std::to_string(quad.npts());
    REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));
  };

  SECTION("First Kind") {
    IntegratorXX::GaussChebyshev1<double, double> quad_even(500);
    integrate(quad_even);
    IntegratorXX::GaussChebyshev1<double, double> quad_odd(501);
    integrate(quad_odd);
  }
  SECTION("Second Kind") {
    IntegratorXX::GaussChebyshev2<double, double> quad_even(500);
    integrate(quad_even);
    IntegratorXX::GaussChebyshev2<double, double> quad_odd(501);
    integrate(quad_odd);
  }
  SECTION("Second Kind (Modified)") {
    IntegratorXX::GaussChebyshev2Modified<double, double> quad_even(500);
    integrate(quad_even);
    IntegratorXX::GaussChebyshev2Modified<double, double> quad_odd(501);
    integrate(quad_odd);
  }
  SECTION("Third Kind") {
    IntegratorXX::GaussChebyshev3<double, double> quad(500);
    integrate(quad);
  }
}
#endif

TEST_CASE( "Euler-Maclaurin Quadratures", "[1d-quad]" ) {

  IntegratorXX::MurrayHandyLaming<double,double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return x*x*gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i )
    res += wgt[i] * f(pts[i]);

  const auto ref = std::sqrt(M_PI)/4;
  const auto msg = "Euler-Maclaurin N = " + std::to_string(quad.npts());
  REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));

}

TEST_CASE( "Treutler-Ahlrichs Quadratures", "[1d-quad]" ) {

  IntegratorXX::TreutlerAhlrichs<double,double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return x*x*gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i ) {
    res += wgt[i] * f(pts[i]);
  }

  const auto ref = std::sqrt(M_PI)/4;
  const auto msg = "Treutler-Ahlrichs N = " + std::to_string(quad.npts());
  REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));

}

TEST_CASE( "Knowles Quadratures", "[1d-quad]" ) {

  IntegratorXX::MuraKnowles<double,double> quad(350);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return x*x*gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i ) 
    res += wgt[i] * f(pts[i]);
  

  const auto ref = std::sqrt(M_PI)/4;
  const auto msg = "Mura-Knowles N = " + std::to_string(quad.npts());
  REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));

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
