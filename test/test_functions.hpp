#pragma once
#include <cmath>
#include <random>
#include "quad_matcher.hpp"
#include <integratorxx/util/factorial.hpp>

namespace detail {
}
struct AssociatedLegendre {
  static inline double evaluate(int l, int m, double x) {
    return std::assoc_legendre(l, m, x);
  }
};

struct SphericalHarmonic {
  static auto evaluate(int l, int m, double theta, double phi) {
    double fac_ratio  = IntegratorXX::factorial(l-m);
           fac_ratio /= IntegratorXX::factorial(l+m);
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
