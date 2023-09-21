#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/util/copy_grid.hpp>
#include <integratorxx/quadratures/s2/delley/delley_grids.hpp>
#include <vector>

namespace IntegratorXX {

/**
 *  @brief Implementation of Delley's reoptimized Lebedev quadratures on the
 * sphere.
 *
 *  Reference:
 *  B. Delley, J. Comput. Chem. 17, 1152 (1996)
 *  doi:10.1002/(SICI)1096-987X(19960715)17:9<1152::AID-JCC7>3.0.CO;2-R
 *
 *  @tparam PointContainer  Quadrature points
 *  @tparam WeightContainer Quadrature weights
 */

namespace detail::delley {

inline static int64_t npts_by_algebraic_order(int64_t order) {
  switch(order) {
    case 5:
      return 14;
    case 7:
      return 26;
    case 11:
      return 50;
    case 17:
      return 110;
    case 23:
      return 194;
    case 29:
      return 302;
    case 35:
      return 434;
    case 41:
      return 590;
    case 47:
      return 770;
    case 53:
      return 974;
    case 59:
      return 1202;
    case 65:
      return 1454;
    case 71:
      return 1730;
    case 77:
      return 2030;
    case 83:
      return 2354;
    case 89:
      return 2702;
    case 95:
      return 3074;
    case 101:
      return 3470;
    default:
      return -1;
  }
}

inline static int64_t algebraic_order_by_npts(int64_t npts) {
  switch(npts) {
    case 14:
      return 5;
    case 26:
      return 7;
    case 50:
      return 11;
    case 110:
      return 17;
    case 194:
      return 23;
    case 302:
      return 29;
    case 434:
      return 35;
    case 590:
      return 41;
    case 770:
      return 47;
    case 974:
      return 53;
    case 1202:
      return 59;
    case 1454:
      return 65;
    case 1730:
      return 71;
    case 2030:
      return 77;
    case 2354:
      return 83;
    case 2702:
      return 89;
    case 3074:
      return 95;
    case 3470:
      return 101;
    default:
      return -1;
  }
}

inline static int64_t next_algebraic_order(int64_t order) {
  if(order <= 5)
    return 5;
  else if(order <= 7)
    return 7;
  else if(order <= 11)
    return 11;
  else if(order <= 17)
    return 17;
  else if(order <= 23)
    return 23;
  else if(order <= 29)
    return 29;
  else if(order <= 35)
    return 35;
  else if(order <= 41)
    return 41;
  else if(order <= 47)
    return 47;
  else if(order <= 53)
    return 53;
  else if(order <= 59)
    return 59;
  else if(order <= 65)
    return 65;
  else if(order <= 71)
    return 71;
  else if(order <= 77)
    return 77;
  else if(order <= 83)
    return 83;
  else if(order <= 89)
    return 89;
  else if(order <= 95)
    return 95;
  else
    return 101;
}
}  // namespace detail::delley

template <typename RealType>
class Delley : public Quadrature<Delley<RealType>> {
  using base_type = Quadrature<Delley<RealType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  Delley(size_t npts) : base_type(npts) {}

  Delley(const Delley &) = default;
  Delley(Delley &&) noexcept = default;
};

template <typename RealType>
struct quadrature_traits<Delley<RealType>> {
  using point_type = cartesian_pt_t<RealType>;
  using weight_type = RealType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts) {
    point_container points(npts);
    weight_container weights(npts);

    using namespace DelleyGrids;

    if(npts == 14)
      detail::copy_grid<delley_14<RealType>>(points, weights);
    else if(npts == 26)
      detail::copy_grid<delley_26<RealType>>(points, weights);
    else if(npts == 50)
      detail::copy_grid<delley_50<RealType>>(points, weights);
    else if(npts == 110)
      detail::copy_grid<delley_110<RealType>>(points, weights);
    else if(npts == 194)
      detail::copy_grid<delley_194<RealType>>(points, weights);
    else if(npts == 302)
      detail::copy_grid<delley_302<RealType>>(points, weights);
    else if(npts == 434)
      detail::copy_grid<delley_434<RealType>>(points, weights);
    else if(npts == 590)
      detail::copy_grid<delley_590<RealType>>(points, weights);
    else if(npts == 770)
      detail::copy_grid<delley_770<RealType>>(points, weights);
    else if(npts == 974)
      detail::copy_grid<delley_974<RealType>>(points, weights);
    else if(npts == 1202)
      detail::copy_grid<delley_1202<RealType>>(points, weights);
    else if(npts == 1454)
      detail::copy_grid<delley_1454<RealType>>(points, weights);
    else if(npts == 1730)
      detail::copy_grid<delley_1730<RealType>>(points, weights);
    else if(npts == 2030)
      detail::copy_grid<delley_2030<RealType>>(points, weights);
    else if(npts == 2354)
      detail::copy_grid<delley_2354<RealType>>(points, weights);
    else if(npts == 2702)
      detail::copy_grid<delley_2702<RealType>>(points, weights);
    else if(npts == 3074)
      detail::copy_grid<delley_3074<RealType>>(points, weights);
    else if(npts == 3470)
      detail::copy_grid<delley_3470<RealType>>(points, weights);

    // Pretabulated weights are missing 4 pi
    for(auto i=0; i < npts; i++)
      weights[i] *= 4.0*M_PI;

    return std::make_tuple(points, weights);
  }
};
}  // namespace IntegratorXX
