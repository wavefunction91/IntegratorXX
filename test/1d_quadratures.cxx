#include "catch2/catch_all.hpp"
#include <integratorxx/quadratures/all.hpp>
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







// Compose one primitive rule with every radial transform. The pairing is only
// legal if the primitive's traits declare bound_inclusive, so instantiating
// these is itself the regression guard for #116.
template <typename Primitive>
void compose_with_all_transforms() {
  constexpr size_t npts = 16;
  CHECK( IntegratorXX::RadialTransformQuadrature<Primitive,
           IntegratorXX::BeckeRadialTraits>(npts).npts() == npts );
  CHECK( IntegratorXX::RadialTransformQuadrature<Primitive,
           IntegratorXX::MuraKnowlesRadialTraits>(npts).npts() == npts );
  CHECK( IntegratorXX::RadialTransformQuadrature<Primitive,
           IntegratorXX::MurrayHandyLamingRadialTraits<2>>(npts).npts() == npts );
  CHECK( IntegratorXX::RadialTransformQuadrature<Primitive,
           IntegratorXX::TreutlerAhlrichsRadialTraits>(npts).npts() == npts );
}

TEST_CASE( "Gauss-Legendre Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const int p = static_cast<int>(c.size());
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
    const int p = static_cast<int>(c.size());
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
    const int p = static_cast<int>(c.size());
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
    const int p = static_cast<int>(c.size());
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
    const int p = static_cast<int>(c.size());
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


TEST_CASE( "Gauss-Chebyshev extremal weights at large N", "[1d-quad]" ) {

  // The extremal weight of each Chebyshev rule is sin(t) at an angle that the
  // generators reach as pi - t. Rounding t to the nearest double is an
  // O(ulp(pi)) absolute error there, which is an O(N eps) *relative* error in
  // the weight: ~1e-13 at N=800, growing with N. Evaluating the sine at the
  // complement, which is exact in the loop index, keeps it at O(eps).
  //
  // The reference is the Taylor series of sin about 0. At these angles
  // (h < 2e-3) the first omitted term is O(h^7) ~ 1e-22 relative, so the
  // series is good to rounding and is independent of std::sin.
  auto sin_series = [](double h) {
    const double h2 = h * h;
    return h * (1.0 - h2 / 6.0 * (1.0 - h2 / 20.0 * (1.0 - h2 / 42.0)));
  };

  // An absolute 10*eps check cannot see this: the weights are themselves
  // O(1/N^2) here, so the error hides below the tolerance.
  const double rel_tol = 1e-14;
  const size_t npts = 800;

  SECTION("T1") {
    // w_0 = (pi/N) sin(h), h = pi/2N, at the node nearest -1
    const double h = M_PI / (2.0 * npts);
    const double ref = (M_PI / npts) * sin_series(h);
    const GaussChebyshev1<double,double> q(npts);
    const auto& wgts = q.weights();
    CHECK_THAT( wgts.front(), Catch::Matchers::WithinRel(ref, rel_tol) );
    CHECK_THAT( wgts.back(),  Catch::Matchers::WithinRel(ref, rel_tol) );
  }

  SECTION("T2") {
    // w_0 = h sin(h), h = pi/(N+1)
    const double h = M_PI / (npts + 1);
    const double ref = h * sin_series(h);
    const GaussChebyshev2<double,double> q(npts);
    const auto& wgts = q.weights();
    CHECK_THAT( wgts.front(), Catch::Matchers::WithinRel(ref, rel_tol) );
    CHECK_THAT( wgts.back(),  Catch::Matchers::WithinRel(ref, rel_tol) );
  }

  SECTION("T3") {
    // w = h sin(2 h) at the node nearest 1, and h sin(h) at the one nearest 0
    const double h = M_PI / (2.0 * npts + 1);
    const GaussChebyshev3<double,double> q(npts);
    const auto& wgts = q.weights();
    CHECK_THAT( wgts.front(), Catch::Matchers::WithinRel(h * sin_series(2.0*h), rel_tol) );
    CHECK_THAT( wgts.back(),  Catch::Matchers::WithinRel(h * sin_series(h),     rel_tol) );
  }

}

TEST_CASE( "Euler-Maclaurin Quadratures by Murray, Handy, and Laming", "[1d-quad]" ) {
  IntegratorXX::MurrayHandyLaming<double,double> quad(150);
  const auto msg = "Euler-Maclaurin N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Treutler-Ahlrichs Quadratures", "[1d-quad]" ) {
  IntegratorXX::TreutlerAhlrichs<double,double> quad(150);
  const auto msg = "Treutler-Ahlrichs N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Mura-Knowles Quadratures", "[1d-quad]" ) {
  IntegratorXX::MuraKnowles<double,double> quad(350);
  const auto msg = "Mura-Knowles N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Becke Quadratures", "[1d-quad]" ) {
  IntegratorXX::Becke<double,double> quad(350);
  const auto msg = "Becke N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}


TEST_CASE( "Radial transforms over the primitive rules", "[1d-quad]" ) {

  // The named radial grids are all fixed aliases over GaussChebyshev2, so
  // nothing above pairs a transform with any other primitive by hand.

  SECTION("Gauss-Lobatto x Becke") {
    // Lobatto is the only bound_inclusive primitive, so this is also the
    // endpoint-dropping path: the base rule is asked for npts+2 nodes and the
    // two that map to r = 0 and r = infinity are skipped.
    IntegratorXX::RadialTransformQuadrature<
      IntegratorXX::GaussLobatto<double,double>,
      IntegratorXX::BeckeRadialTraits> quad(350);
    REQUIRE( quad.npts() == 350 );
    const auto msg = "Gauss-Lobatto x Becke N = " + std::to_string(quad.npts());
    test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
  }

  SECTION("Gauss-Legendre x Becke") {
    IntegratorXX::RadialTransformQuadrature<
      IntegratorXX::GaussLegendre<double,double>,
      IntegratorXX::BeckeRadialTraits> quad(350);
    REQUIRE( quad.npts() == 350 );
    const auto msg = "Gauss-Legendre x Becke N = " + std::to_string(quad.npts());
    test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
  }

  SECTION("Every primitive composes with every transform") {
    // All 28 pairings. Accuracy is not the point here: several of these are
    // not sensible rules, but all of them must build.
    compose_with_all_transforms<IntegratorXX::GaussChebyshev1<double,double>>();
    compose_with_all_transforms<IntegratorXX::GaussChebyshev2<double,double>>();
    compose_with_all_transforms<IntegratorXX::GaussChebyshev2Modified<double,double>>();
    compose_with_all_transforms<IntegratorXX::GaussChebyshev3<double,double>>();
    compose_with_all_transforms<IntegratorXX::GaussLegendre<double,double>>();
    compose_with_all_transforms<IntegratorXX::GaussLobatto<double,double>>();
    compose_with_all_transforms<IntegratorXX::UniformTrapezoid<double,double>>();
  }

}

TEST_CASE( "Lebedev-Laikov", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {
    IntegratorXX::LebedevLaikov<double> quad( nPts );
    const auto msg = "Lebedev-Laikov N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);
  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}

TEST_CASE( "Ahrens-Beylkin", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::AhrensBeylkin<double> quad( nPts );
    const auto msg = "Ahrens-Beylkin N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(312);
  test_fn(792);
  test_fn(972);
}

TEST_CASE( "Womersley", "[1d-quad]" ) {

  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::Womersley<double> quad( nPts );
    const auto msg = "Womersley N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(314);
  test_fn(801);
  test_fn(969);
}

TEST_CASE( "Delley", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::Delley<double> quad( nPts );
    const auto msg = "Delley N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}
