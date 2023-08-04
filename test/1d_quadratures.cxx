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

#include "quad_matcher.hpp"
#include "test_functions.hpp"

#if 0
// NR
template <typename T, typename IntT>
std::enable_if_t< std::is_arithmetic_v<T>, T>
assoc_legendre( const IntT l, const IntT m, const T x ) {

  if( m < 0 or m > l or std::abs(x) > 1. )
    throw std::runtime_error("Bad args to assoc_legendre");

  T pmm = 1.;
  if( m > 0 ) {

    const T somx2 = std::sqrt( (1. - x)*(1. + x) );
    T fact = 1.;
    for( IntT i = 0; i < m; ++i ) {
      pmm *= fact*somx2;
      fact += 2.;
    }

  }

  if( l == m ) return pmm;
  else {
    T pmmp1 = x * (2*m + 1) * pmm;

    if( l == (m+1) ) return pmmp1;
    else {

      T pll;
      for( IntT ll = m+1; ll < l; ++ll ) {
        pll = (x*(2*ll+1)*pmmp1-(ll+m)*pmm)/(ll-m+1);
        pmm = pmmp1;
        pmmp1 = pll;
      }

      return pll;
    }
  }

}
#endif


using namespace IntegratorXX;

inline constexpr double inf = std::numeric_limits<double>::infinity();

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


#if 0
template <typename T = double>
T ref_bessel_int( int n, T a, T b ) {
  return std::pow(2.0,-n) * (std::pow(-1.0,n) +1) * boost::math::hypergeometric_pFq( {(n+1.0)/2.0}, {(n+3.0)/2.0, n+1.0}, -0.25) / std::tgamma(n+2);
};
#endif 

template <typename T>
auto chebyshev_T(int n, T x) {
  return std::cos( n * std::acos(x) );
}


constexpr size_t factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


#if 0
template <typename T>
constexpr std::complex<T> spherical_harmonics( int64_t l, int64_t m, T theta, T phi ) {


  auto prefactor = std::sqrt( (2.*l + 1) / (4.*M_PI) * T(factorial(l-m))/T(factorial(l+m)) );
  if( m < 0 )
    prefactor *= std::pow(-1,std::abs(m)) * T(factorial(l-m))/T(factorial(l+m));


  return prefactor * assoc_legendre( l, std::abs(m), std::cos(theta) ) *
         std::complex<T>( std::cos(m*phi), std::sin(m*phi) );

}
#endif

template <typename T>
constexpr std::complex<T> spherical_harmonics( int64_t l, int64_t m, T x, T y, T z ) {

  const auto r     = std::sqrt( x*x + y*y + z*z );
  const auto theta = std::acos( z / r );
  const auto phi   = std::atan2( y, x );

  //return spherical_harmonics( l, m, theta, phi  );
  return SphericalHarmonic::evaluate(l, m, theta, phi);

}


template <typename T, typename... Args>
constexpr T real_spherical_harmonics( int64_t l, int64_t m, Args&&... args ) {

  const T phase = std::sqrt(2.) * ((m % 2) ? 1. : -1.);

  if( m == 0 )
    return std::real(spherical_harmonics( l, m, std::forward<Args>(args)... ));
  else if( m < 0 )
    return phase * std::imag(spherical_harmonics( l, std::abs(m), std::forward<Args>(args)... ));
  else
    return phase * std::real(spherical_harmonics( l, m, std::forward<Args>(args)... ));

}





const double x_tolerance = 10 * std::numeric_limits<double>::epsilon();
const double w_tolerance = 10 * std::numeric_limits<double>::epsilon();

TEST_CASE( "Gauss-Legendre Quadratures", "[1d-quad]" ) {

  for(unsigned order=10;order<14;order++) {
    std::ostringstream oss;
    oss << "order " << order;
    SECTION( oss.str()) {

      IntegratorXX::GaussLegendre<double,double> quad( order );

      const auto& pts = quad.points();
      const auto& wgt = quad.weights();

      auto f = [=]( double x ){ return gaussian(x); };

      double res = 0.;
      for( auto i = 0; i < quad.npts(); ++i ) {
        res += wgt[i] * f(pts[i]);
      }

      const auto ref = ref_gaussian_int(-1.,1.);
      const auto msg = "Gauss-Legendre N = " + std::to_string(quad.npts());
      REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));
    }
  }
}

TEST_CASE( "Gauss-Lobatto Quadratures", "[1d-quad]" ) {

  for(unsigned order=10;order<14;order++) {
    std::ostringstream oss;
    oss << "order " << order;
    SECTION( oss.str()) {

      IntegratorXX::GaussLobatto<double,double> quad( order );

      const auto& pts = quad.points();
      const auto& wgt = quad.weights();

      auto f = [=]( double x ){ return gaussian(x); };

      double res = 0.;
      for( auto i = 0; i < quad.npts(); ++i ) {
        res += wgt[i] * f(pts[i]);
      }

      const auto ref = ref_gaussian_int(-1.,1.);
      const auto msg = "Gauss-Lobatto N = " + std::to_string(quad.npts());
      REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, 1e-10));
    }
  }
}

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

      auto f = [=]( decltype(pts[0]) x ){
        return spherical_harmonics(l,m,x[0],x[1],x[2]);
      };

      std::complex<double> res = 0.;
      for( auto i = 0; i < quad.npts(); ++i )
        res += wgt[i] * f(pts[i])* std::conj(f(pts[i]));

      const auto msg = "Ahrens-Beylkin N = " + std::to_string(quad.npts());
      REQUIRE_THAT( std::real(res), IntegratorXX::Matchers::WithinAbs(msg, 1.0, 1e-10) );
      REQUIRE_THAT( std::imag(res), IntegratorXX::Matchers::WithinAbs(msg, 0.0, 1e-10) );

    }

  };

  test_fn(312);
  test_fn(792);
  test_fn(972);
}
