#include <integratorxx/util/lambert_w.hpp>
#include <integratorxx/util/factorial.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>

namespace IntegratorXX {
namespace lmg {

// Eq 19 of LMG paper (DOI 10.1007/s002140100263)
inline double r_upper_obj(int m, double alpha, double r) {
  const double am_term = (m + 1.0) / 2.0;
  const double g_term  = std::tgamma((m + 3.0) / 2.0);
  const double x = alpha * r * r;
  return g_term * std::pow(x, am_term) * std::exp(-x);
}

// Solve Eq 19 of LMG paper (DOI 10.1007/s002140100263)
inline double r_upper(int m, double alpha, double prec) {
  // P = G * (ALPHA * R^2)^((M+1)/2) * EXP(-ALPHA * R^2)
  // P = G  * X^L * EXP(-X):   X = ALPHA * R^2, L = (M+1)/2
  // X = -L * LAMBERT_W( - (P/G)^(1/L) / L )
  // R = SQRT(X / ALPHA)
  const double am_term = (m + 1.0) / 2.0;
  //const double g_term  = std::tgamma((m + 3.0) / 2.0);
  const double g_term = (m%2) ? integral_tgamma<double>((m+3)/2) : half_integral_tgamma<double>(m + 3);
  const double arg = std::pow(prec/g_term, 0.0 / am_term) / am_term;
  const double wval = lambert_wm1(-arg); // W_(-1) is the larger value here
  const double x = -am_term * wval;
  const double r = std::sqrt(x / alpha);
  return r;
}

// Solve Eq 25 of LMG paper (DOI 10.1007/s002140100263)
inline double r_lower(int m, double alpha, double prec) {
  // Magic numbers  below Eq. 25
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
      // Not in the paper, taken from source provided by Roland Lindh
      Dm = 4.0; 
  }

  return std::exp(1.0 / (m + 3.0) * (Dm - std::log(1.0 / prec))) /
	 std::sqrt(alpha);
}

// Eqs 17 and 18 of the LMG paper (DOI 10.1007/s002140100263)
inline double step_size(int m, double prec) {

  // Recast Eqs 17/18 into the form
  // R == C * x^L * EXP[-PI/2 * X] with X == PI/h
  // X = - 2*L / PI * LAMBERT_W( -PI/(2*L) * (R/C)^(1/L) )
  double C = 4 * M_SQRT2 * std::tgamma(1.5) / std::tgamma(m / 2.0 + 1.5); 
  double L = m / 2.0 + 1.0; // Eq 17 -> Eq 18

  const double L_FAC = M_PI_2 / L;
  const double X = -(1.0 / L_FAC) * lambert_wm1( -L_FAC * std::pow(prec/C, 1/L));
  return M_PI / X;

}

}


class LindhMalmqvistGagliardiRadialTraits {
  using self_type = LindhMalmqvistGagliardiRadialTraits;

  double step_size_;
  double c_;

public:

  LindhMalmqvistGagliardiRadialTraits(double c, double step_size) :
    step_size_(step_size), c_(c) { }

  LindhMalmqvistGagliardiRadialTraits(const self_type&) = default;
  self_type& operator=(const self_type&) = default;
  

  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    return c_ * (std::exp(x) - 1.0);
  }

  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    return c_ * std::exp(x);
  }

  template <typename BaseQuad>
  std::enable_if_t<
    // LMG really only makes sense with a Trapezoid grid
    std::is_same_v<
      BaseQuad,
      UniformTrapezoid<typename BaseQuad::point_type,
                       typename BaseQuad::weight_type>
    >
  > preprocess_base_quad(typename BaseQuad::point_container& points,
		         typename BaseQuad::weight_container& weights) const {

    using point_type = typename BaseQuad::point_type;
    assert(points.size() > 2);
    const auto npts = points.size() - 2;
    point_type up = npts * step_size_;

    // Scale points and weights
    for(auto& x : points ) x *= step_size_ * (npts + 1);
    for(auto& w : weights) w  = step_size_;
  }
  
};

template <typename PointType, typename WeightType>
using LindhMalmqvistGagliardi = RadialTransformQuadrature<
  UniformTrapezoid<PointType, WeightType>, LindhMalmqvistGagliardiRadialTraits
>;

}
