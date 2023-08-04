#pragma once
#include <cmath>
#include "quad_matcher.hpp"


namespace detail {
  constexpr size_t factorial(size_t n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }
}
struct AssociatedLegendre {
  static inline double evaluate(int l, int m, double x) {
    return std::assoc_legendre(l, m, x);
  }
};

struct SphericalHarmonic {
  static auto evaluate(int l, int m, double theta, double phi) {
    double fac_ratio  = detail::factorial(l-m);
           fac_ratio /= detail::factorial(l+m);
    double prefactor  = (2.0 * l + 1.0) / (4.0 * M_PI);
           prefactor *= fac_ratio;
    prefactor = std::sqrt(prefactor);
    if( m < 0 ) {
      prefactor *= std::pow(1, std::abs(m)) * fac_ratio;
    }

    return prefactor * 
           std::assoc_legendre(l, std::abs(m), std::cos(theta)) *
           std::complex( std::cos(m*phi), std::sin(m*phi) );
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


template <typename TestFunction, typename QuadType, typename... PreArgs, 
  typename ResType = std::invoke_result_t<PreArgs..., typename QuadType::point_type>
>
void test_quadrature(std::string msg, const QuadType& quad, const ResType& ref, double e, PreArgs&&... args) {
  const auto& pts = quad.points();
  const auto& wgt = quad.weights();

  ResType res = 0.0;
  for(auto i = 0; i < quad.npts(); ++i) {
    res += wgt[i] * TestFunction::evaluate(args..., pts[i]);
  }
  
  //standard_matcher(mes, res, ref, e);
  //printf("diff = %.6e\n", std::abs(ref - res));
  REQUIRE_THAT(res, IntegratorXX::Matchers::WithinAbs(msg, ref, e));
}
