#include <integratorxx/util/lambert_w.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>

namespace IntegratorXX {
namespace lmg {

inline double r_upper_obj(int m, double alpha, double r) {
  const double am_term = (m + 1.0) / 2.0;
  const double g_term  = std::tgamma((m + 3.0) / 2.0);
  const double x = alpha * r * r;
  return g_term * std::pow(x, am_term) * std::exp(-x);
}

// Equation 19 in the LMG paper
inline double r_upper(int m, double alpha, double prec) {
  // P = G * (ALPHA * R^2)^((M+1)/2) * EXP(-ALPHA * R^2)
  // P = G  * X^L * EXP(-X):   X = ALPHA * R^2, L = (M+1)/2
  // X = -L * LAMBERT_W( - (P/G)^(1/L) / L )
  // R = SQRT(X / ALPHA)
  const double am_term = (m + 1.0) / 2.0;
  const double g_term  = std::tgamma((m + 3.0) / 2.0);
  const double arg = std::pow(prec/g_term, 1.0 / am_term) / am_term;
  const double wval = lambert_wm1(-arg); // W_(-1) is the larger value here
  const double x = -am_term * wval;
  const double r = std::sqrt(x / alpha);
  return r;
}

inline double r_lower(int m, double alpha, double prec) {
  // Magic numbers from the LMG paper below Eq. 25
  double Dm;
  switch(m) {
    case -2:
      Dm = 9.1;
      break;
    case 0:
      Dm = 1.9;
      break;
    case 2:
      Dm = -1.0;
      break;
    case 4:
      Dm = -2.3;
      break;
    default:
      Dm = NAN;
  }

  return std::exp(1.0 / (m + 3.0) * (Dm - std::log(1.0 / prec))) /
	 std::sqrt(alpha);
}

}
}
