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
      Dm = 4.0; // Not in the paper, taken from soruce provided by Roland Lindh
  }

  return std::exp(1.0 / (m + 3.0) * (Dm - std::log(1.0 / prec))) /
	 std::sqrt(alpha);
}

inline double step_size(int m, double prec) {

  // Recast Eqs 17/18 into the form
  // R == C * x^L * EXP[-PI/2 * X] with X == PI/h
  // X = - 2*L / PI * LAMBERT_W( -PI/(2*L) * (R/C)^(1/L) )
  double C = 4 * M_SQRT2; // This is the constant in Eq 17
  double L;
  if( m == 0 ) {
    L = 1;
  } else {
    L = m / 2.0 + 1.0; // Eq 17 -> Eq 18
    C = C * std::tgamma(1.5) / std::tgamma((m + 3.0) / 2.0); 
  }

  const double L_FAC = M_PI_2 / L;
  const double X = -(1.0 / L_FAC) * lambert_wm1( -L_FAC * std::pow(prec/C, 1/L));
  return M_PI / X;

}

}
}
