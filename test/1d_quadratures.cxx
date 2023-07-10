#include "catch2/catch_all.hpp"
#include <integratorxx/quadratures/gausslegendre.hpp>
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


  return prefactor * assoc_legendre( l, std::abs(m), std::cos(theta) ) *
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

      CHECK( res == Catch::Approx(ref_gaussian_int(-1.,1.)) );
    }
  }
}

TEST_CASE( "Gauss-Chebyshev Quadratures", "[1d-quad]") {

  constexpr unsigned order = 200;
  auto integrate = [&](auto& quad) {
    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    // Check that nodes are in increasing value
    for(auto i = 1; i < quad.npts(); ++i) {
      if(pts[i] <= pts[i - 1]) {
        std::ostringstream oss;
        oss << "Quadrature points are not in increasing order: " << pts[i-1] << " " << pts[i] << "!\n";
	throw std::runtime_error(oss.str());
      }
    }

    auto f = [=]( double x ){ return gaussian(x); };

    double res = 0.;
    for( auto i = 0; i < quad.npts(); ++i ) {
      res += wgt[i] * f(pts[i]);
    }

    CHECK( res == Catch::Approx(ref_gaussian_int(-1.,1.)) );
  };

  SECTION("First Kind") {
    IntegratorXX::GaussChebyshev1<double, double> quad(order);
    integrate(quad);
  }
  SECTION("Second Kind") {
    IntegratorXX::GaussChebyshev2<double, double> quad(order);
    integrate(quad);
  }
  SECTION("Second Kind (Modified)") {
    IntegratorXX::GaussChebyshev2Modified<double, double> quad(order);
    integrate(quad);
  }
  SECTION("Third Kind") {
    IntegratorXX::GaussChebyshev3<double, double> quad(order);
    integrate(quad);
  }
}

TEST_CASE( "Euler-Maclaurin Quadratures", "[1d-quad]" ) {

  IntegratorXX::MurrayHandyLaming<double,double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i )
    res += wgt[i] * f(pts[i]);

  CHECK( res == Catch::Approx(ref_gaussian_int(0.,inf)) );

}

TEST_CASE( "Treutler-Ahlrichs Quadratures", "[1d-quad]" ) {

  IntegratorXX::TreutlerAhlrichs<double,double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i ) {
    res += wgt[i] * f(pts[i]);
  }

  CHECK( res == Catch::Approx(ref_gaussian_int(0.,inf)) );

}

TEST_CASE( "Knowles Quadratures", "[1d-quad]" ) {

  IntegratorXX::MuraKnowles<double,double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=]( double x ){ return gaussian(x); };

  double res = 0.;
  for( auto i = 0; i < quad.npts(); ++i ) 
    res += wgt[i] * f(pts[i]);
  

  CHECK( res == Catch::Approx(ref_gaussian_int(0.,inf)) );

}

TEST_CASE( "Lebedev-Laikov", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::LebedevLaikov<double> quad( nPts );

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
  
      CHECK( std::real(res) == Catch::Approx(1.) );
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

      CHECK( std::real(res) == Catch::Approx(1.) );

    }

  };

  test_fn(312);
  test_fn(792);
  test_fn(972);
}
