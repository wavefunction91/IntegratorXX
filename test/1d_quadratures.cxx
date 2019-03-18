#include "catch2/catch.hpp"
#include <integratorxx/quadrature.hpp>
#include <cmath>
#include <complex>

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


constexpr size_t factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename T>
constexpr std::complex<T> spherical_harmonics( int64_t l, int64_t m, T theta, T phi ) {


  auto prefactor = std::sqrt( (2.*l + 1) / (4.*M_PI) * T(factorial(l-m))/T(factorial(l+m)) );
  if( m < 0 )
    prefactor *= std::pow(-1,std::abs(m)) * T(factorial(l-m))/T(factorial(l+m));


  return prefactor * std::assoc_legendre( l, std::abs(m), std::cos(theta) ) * 
         std::complex<T>( std::cos(m*phi), std::sin(m*phi) );

}

template <typename T>
constexpr std::complex<T> spherical_harmonics( int64_t l, int64_t m, T x, T y, T z ) {

  const auto r = std::sqrt( x*x + y*y + z*z );

  return spherical_harmonics( l, m, std::acos( z / r ), std::atan2( y, x ) );

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


TEST_CASE( "Gauss-Legendre Quadratures", "[1d-quad]" ) {

  constexpr unsigned order = 10;
  GaussLegendre<double> quad( order, -1, 1 );

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < pts.size(); ++i )
    res += wgt[i] * f(pts[i]);

  CHECK( res == Approx(ref_gaussian_int(-1.,1.)) );

}

TEST_CASE( "Euler-Maclaurin Quadratures", "[1d-quad]" ) {

  EulerMaclaurin<double> quad( 150 );

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < pts.size(); ++i )
    res += wgt[i] * f(pts[i]);

  CHECK( res == Approx(ref_gaussian_int(0.,inf)) );

}

TEST_CASE( "Lebedev-Laikov", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {
  
    Lebedev<double> quad( nPts );
  
    const auto& pts = quad.points();
    const auto& wgt = quad.weights();
  
    for( auto l = 1; l < 10; ++l ) 
    for( auto m = 0; m <= l; ++m ) {
  
      auto f = [=]( decltype(pts[0]) x ){ 
        return spherical_harmonics(l,m,x[0],x[1],x[2]); 
      };
  
      std::complex<double> res = 0.;
      for( auto i = 0; i < pts.size(); ++i )
        res += wgt[i] * f(pts[i])* std::conj(f(pts[i]));
  
      CHECK( 4.*M_PI*std::real(res) == Approx(1.) );
  
    }

  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}
