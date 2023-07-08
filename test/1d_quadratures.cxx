#include <cmath>
#include <complex>
#include <integratorxx/quadratures/gausscheby1.hpp>
#include <integratorxx/quadratures/gausscheby2.hpp>
#include <integratorxx/quadratures/gausscheby2_mod.hpp>
#include <integratorxx/quadratures/gausscheby3.hpp>
#include <integratorxx/quadratures/gausslegendre.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/quadratures/mhl.hpp>
#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/treutlerahlrichs.hpp>

#include "catch2/catch_all.hpp"

// NR
template <typename T, typename IntT>
std::enable_if_t<std::is_arithmetic_v<T>, T> assoc_legendre(const IntT l,
                                                            const IntT m,
                                                            const T x) {
  if(m < 0 or m > l or std::abs(x) > 1.)
    throw std::runtime_error("Bad args to assoc_legendre");

  T pmm = 1.;
  if(m > 0) {
    const T somx2 = std::sqrt((1. - x) * (1. + x));
    T fact = 1.;
    for(IntT i = 0; i < m; ++i) {
      pmm *= fact * somx2;
      fact += 2.;
    }
  }

  if(l == m)
    return pmm;
  else {
    T pmmp1 = x * (2 * m + 1) * pmm;

    if(l == (m + 1))
      return pmmp1;
    else {
      T pll;
      for(IntT ll = m + 1; ll < l; ++ll) {
        pll = (x * (2 * ll + 1) * pmmp1 - (ll + m) * pmm) / (ll - m + 1);
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
constexpr T gaussian(T alpha, T c, T x) {
  return std::exp(-alpha * (x - c) * (x - c));
}

template <typename T = double>
constexpr T gaussian(T x) {
  return gaussian(1., 0., x);
}

template <typename T = double>
constexpr T ref_gaussian_int(T alpha, T c, T a, T b) {
  const auto low = (a == inf) ? T(-1.) : std::erf(std::sqrt(alpha) * (c - a));
  const auto hgh = (b == inf) ? T(-1.) : std::erf(std::sqrt(alpha) * (c - b));

  return 0.5 * std::sqrt(M_PI / alpha) * (low - hgh);
}

template <typename T = double>
constexpr T ref_gaussian_int(T a, T b) {
  return ref_gaussian_int(1., 0., a, b);
}

constexpr size_t factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename T>
constexpr std::complex<T> spherical_harmonics(int64_t l, int64_t m, T theta,
                                              T phi) {
  auto prefactor = std::sqrt((2. * l + 1) / (4. * M_PI) * T(factorial(l - m)) /
                             T(factorial(l + m)));
  if(m < 0)
    prefactor *=
        std::pow(-1, std::abs(m)) * T(factorial(l - m)) / T(factorial(l + m));

  return prefactor * assoc_legendre(l, std::abs(m), std::cos(theta)) *
         std::complex<T>(std::cos(m * phi), std::sin(m * phi));
}

template <typename T>
constexpr std::complex<T> spherical_harmonics(int64_t l, int64_t m, T x, T y,
                                              T z) {
  const auto r = std::sqrt(x * x + y * y + z * z);

  return spherical_harmonics(l, m, std::acos(z / r), std::atan2(y, x));
}

template <typename T, typename... Args>
constexpr T real_spherical_harmonics(int64_t l, int64_t m, Args&&... args) {
  const T phase = std::sqrt(2.) * ((m % 2) ? 1. : -1.);

  if(m == 0)
    return std::real(spherical_harmonics(l, m, std::forward<Args>(args)...));
  else if(m < 0)
    return phase * std::imag(spherical_harmonics(l, std::abs(m),
                                                 std::forward<Args>(args)...));
  else
    return phase *
           std::real(spherical_harmonics(l, m, std::forward<Args>(args)...));
}

TEST_CASE("Gauss-Legendre Quadratures", "[1d-quad]") {
  constexpr unsigned order = 10;

  SECTION("untransformed bounds") {
    IntegratorXX::GaussLegendre<double, double> quad(order, -1, 1);

    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    auto f = [=](double x) { return gaussian(x); };

    double res = 0.;
    for(auto i = 0; i < quad.npts(); ++i) res += wgt[i] * f(pts[i]);

    CHECK(res == Catch::Approx(ref_gaussian_int(-1., 1.)));
  }

  SECTION("transformed bounds") {
    const double lo = 0.;
    const double up = 4.;
    IntegratorXX::GaussLegendre<double, double> quad(100, lo, up);

    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    auto f = [=](double x) { return gaussian(x); };

    double res = 0.;
    for(auto i = 0; i < quad.npts(); ++i) res += wgt[i] * f(pts[i]);

    CHECK(res == Catch::Approx(ref_gaussian_int(lo, up)));
  }
}

TEST_CASE("Euler-Maclaurin Quadratures", "[1d-quad]") {
  IntegratorXX::MurrayHandyLaming<double, double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=](double x) { return gaussian(x); };

  double res = 0.;
  for(auto i = 0; i < quad.npts(); ++i) res += wgt[i] * f(pts[i]);

  CHECK(res == Catch::Approx(ref_gaussian_int(0., inf)));
}

TEST_CASE("Ahlrichs Quadratures", "[1d-quad]") {
  IntegratorXX::TreutlerAhlrichs<double, double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=](double x) { return gaussian(x); };

  double res = 0.;
  for(auto i = 0; i < quad.npts(); ++i) {
    // std::cout << wgt[i] << ", " << pts[i] << std::endl;
    res += wgt[i] * f(pts[i]);
  }

  CHECK(res == Catch::Approx(ref_gaussian_int(0., inf)));
}

TEST_CASE("Knowles Quadratures", "[1d-quad]") {
  IntegratorXX::MuraKnowles<double, double> quad(150);

  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  auto f = [=](double x) { return gaussian(x); };

  double res = 0.;
  for(auto i = 0; i < quad.npts(); ++i) res += wgt[i] * f(pts[i]);

  CHECK(res == Catch::Approx(ref_gaussian_int(0., inf)));
}

TEST_CASE("Lebedev-Laikov", "[1d-quad]") {
  auto test_fn = [&](size_t nPts) {
    IntegratorXX::LebedevLaikov<double> quad(nPts);

    const auto& pts = quad.points();
    const auto& wgt = quad.weights();

    for(auto l = 1; l < 10; ++l)
      for(auto m = 0; m <= l; ++m) {
        auto f = [=](decltype(pts[0]) x) {
          return spherical_harmonics(l, m, x[0], x[1], x[2]);
        };

        std::complex<double> res = 0.;
        for(auto i = 0; i < quad.npts(); ++i)
          res += wgt[i] * f(pts[i]) * std::conj(f(pts[i]));

        CHECK(4. * M_PI * std::real(res) == Catch::Approx(1.));
      }
  };

  test_fn(302);
  test_fn(770);
  test_fn(974);
}
