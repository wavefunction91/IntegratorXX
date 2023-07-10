#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/gausslegendre.hpp>

namespace IntegratorXX {

template <typename PointType, typename WeightType>
class GaussLobatto : public Quadrature<GaussLobatto<PointType, WeightType>> {
  using base_type = Quadrature<GaussLobatto<PointType, WeightType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  GaussLobatto(size_t npts) : base_type(npts) {}

  GaussLobatto(const GaussLobatto&) = default;
  GaussLobatto(GaussLobatto&&) noexcept = default;
};

template <typename PointType, typename WeightType>
struct quadrature_traits<GaussLobatto<PointType, WeightType>> {
  using point_type = PointType;
  using weight_type = WeightType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts) {
    point_container points(npts);
    weight_container weights(npts);

    // Absolute precision for the nodes
    const auto eps = std::numeric_limits<double>::epsilon();

    // 2/(n(n-1)) appears in many expressions
    weight_type two_ov_nnm1 = 2.0 / (npts * (npts - 1.0));

    // Since the rules are symmetric around the origin, we only need
    // to compute one half of the points
    const size_t mid = (npts+1) / 2;
    for(size_t idx = 1; idx < mid; ++idx) {
      // Initial guess
      const point_type i = npts - 1 - idx;
      point_type z = cos (i * M_PI / ( npts - 1.0));

      // Old value of root
      point_type z_old = -2.0;
      // Values of P_{n-1}(x), dP_{n-1}/dx, P_{n-2}(x)
      point_type p_nm1, p_nm2;

      // Solve the root to eps absolute precision
      bool converged = false;
      const int maxit = 100;
      for(int it = 0; it < maxit; ++it) {
        // Evaluate the Legendre polynomial at z and its derivative
        std::tie(p_nm1, std::ignore, p_nm2) = IntegratorXX::eval_Pn(z, npts - 1);

        // Newton update for the root. This equation might look
        // peculiar, but it is correct: you can derive it by using
        // Newton's method to search for the roots of f(x) = (x^2-1)
        // P_{n}'(x), and using the differential equation satisfied by
        // the Legendre polynomials to find that f'(x) = n(n+1)
        // P_{n}(x).
        z_old = z;
        z -= (z*p_nm1 - p_nm2) / (npts*p_nm1);
        // Convergence check
        if(std::abs(z - z_old) <= eps) {
          converged = true;
          break;
        }
      }  // end while

      if(not converged) {
        throw std::runtime_error(
            "Gauss-Lobatto Newton Iterations Failed to Converge");
      }

      // Quadrature node is z
      point_type pt = z;
      // Quadrature weight is 2 / [n (n-1) (1-x^2) P_{n-1}(x)^2]
      weight_type wgt = two_ov_nnm1 / (p_nm1 * p_nm1);

      // Store the symmetric points
      points[idx] = pt;
      weights[idx] = wgt;

      points[npts - 1 - idx] = -pt;
      weights[npts - 1 - idx] = wgt;
    }  // Loop over points

    // Nodes and weights of the end points
    points[0] = -1.0;
    weights[0] = two_ov_nnm1;
    points[npts - 1] = 1.0;
    weights[npts - 1] = two_ov_nnm1;

    if(npts%2==1) {
      // Rule with even number of points has a node at the origin
      points[npts/2]=0.0;
    }

    return std::make_tuple(points, weights);
  }
};

}  // namespace IntegratorXX
