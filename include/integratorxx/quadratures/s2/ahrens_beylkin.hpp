#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/s2/ahrens_beylkin/ahrens_beylkin_grids.hpp>
#include <integratorxx/util/copy_grid.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of Ahrens-Beylkin quadrature on the sphere.
 *
 *  Reference:
 *  Proc. R. Soc. A (2009) 465, 3103.
 *  DOI: https://doi.org/10.1098/rspa.2009.0104
 *
 *  @tparam PointContainer  Quadrature points
 *  @tparam WeightContainer Quadrature weights
 */

namespace detail::ahrensbeylkin {

inline static int64_t npts_by_algebraic_order(int64_t order) {
  switch (order) {
  case 14:
    return 72;
  case 19:
    return 132;
  case 23:
    return 192;
  case 29:
    return 312;
  case 32:
    return 372;
  case 35:
    return 432;
  case 37:
    return 492;
  case 39:
    return 552;
  case 41:
    return 612;
  case 44:
    return 672;
  case 45:
    return 732;
  case 47:
    return 792;
  case 51:
    return 912;
  case 53:
    return 972;
  case 54:
    return 1032;
  case 56:
    return 1092;
  case 60:
    return 1272;
  case 62:
    return 1332;
  case 63:
    return 1392;
  case 66:
    return 1512;
  case 67:
    return 1572;
  case 69:
    return 1632;
  case 70:
    return 1692;
  case 72:
    return 1812;
  case 75:
    return 1932;
  case 76:
    return 1992;
  case 78:
    return 2112;
  case 79:
    return 2172;
  case 80:
    return 2232;
  case 81:
    return 2292;
  case 84:
    return 2412;
  case 85:
    return 2472;
  case 89:
    return 2712;
  case 90:
    return 2772;
  case 91:
    return 2832;
  case 94:
    return 3012;
  case 98:
    return 3312;
  case 99:
    return 3372;
  case 100:
    return 3432;
  case 101:
    return 3492;
  case 103:
    return 3612;
  case 105:
    return 3792;
  case 111:
    return 4212;
  case 113:
    return 4332;
  case 116:
    return 4572;
  case 117:
    return 4692;
  case 121:
    return 4992;
  case 124:
    return 5232;
  case 127:
    return 5472;
  case 128:
    return 5592;
  case 130:
    return 5772;
  case 135:
    return 6192;
  case 136:
    return 6312;
  case 143:
    return 6912;
  case 149:
    return 7512;
  case 210:
    return 15012;
  default:
    return -1;
  }
}

inline static int64_t algebraic_order_by_npts(int64_t npts) {
  switch (npts) {
  case 72:
    return 14;
  case 132:
    return 19;
  case 192:
    return 23;
  case 312:
    return 29;
  case 372:
    return 32;
  case 432:
    return 35;
  case 492:
    return 37;
  case 552:
    return 39;
  case 612:
    return 41;
  case 672:
    return 44;
  case 732:
    return 45;
  case 792:
    return 47;
  case 912:
    return 51;
  case 972:
    return 53;
  case 1032:
    return 54;
  case 1092:
    return 56;
  case 1272:
    return 60;
  case 1332:
    return 62;
  case 1392:
    return 63;
  case 1512:
    return 66;
  case 1572:
    return 67;
  case 1632:
    return 69;
  case 1692:
    return 70;
  case 1812:
    return 72;
  case 1932:
    return 75;
  case 1992:
    return 76;
  case 2112:
    return 78;
  case 2172:
    return 79;
  case 2232:
    return 80;
  case 2292:
    return 81;
  case 2412:
    return 84;
  case 2472:
    return 85;
  case 2712:
    return 89;
  case 2772:
    return 90;
  case 2832:
    return 91;
  case 3012:
    return 94;
  case 3312:
    return 98;
  case 3372:
    return 99;
  case 3432:
    return 100;
  case 3492:
    return 101;
  case 3612:
    return 103;
  case 3792:
    return 105;
  case 4212:
    return 111;
  case 4332:
    return 113;
  case 4572:
    return 116;
  case 4692:
    return 117;
  case 4992:
    return 121;
  case 5232:
    return 124;
  case 5472:
    return 127;
  case 5592:
    return 128;
  case 5772:
    return 130;
  case 6192:
    return 135;
  case 6312:
    return 136;
  case 6912:
    return 143;
  case 7512:
    return 149;
  case 15012:
    return 210;
  default:
    return -1;
  }
}

inline static int64_t next_algebraic_order(int64_t order) {

  if (order <= 14)
    return 14;
  else if (order <= 19)
    return 19;
  else if (order <= 23)
    return 23;
  else if (order <= 29)
    return 29;
  else if (order <= 32)
    return 32;
  else if (order <= 35)
    return 35;
  else if (order <= 37)
    return 37;
  else if (order <= 39)
    return 39;
  else if (order <= 41)
    return 41;
  else if (order <= 44)
    return 44;
  else if (order <= 45)
    return 45;
  else if (order <= 47)
    return 47;
  else if (order <= 51)
    return 51;
  else if (order <= 53)
    return 53;
  else if (order <= 54)
    return 54;
  else if (order <= 56)
    return 56;
  else if (order <= 60)
    return 60;
  else if (order <= 62)
    return 62;
  else if (order <= 63)
    return 63;
  else if (order <= 66)
    return 66;
  else if (order <= 67)
    return 67;
  else if (order <= 69)
    return 69;
  else if (order <= 70)
    return 70;
  else if (order <= 72)
    return 72;
  else if (order <= 75)
    return 75;
  else if (order <= 76)
    return 76;
  else if (order <= 78)
    return 78;
  else if (order <= 79)
    return 79;
  else if (order <= 80)
    return 80;
  else if (order <= 81)
    return 81;
  else if (order <= 84)
    return 84;
  else if (order <= 85)
    return 85;
  else if (order <= 89)
    return 89;
  else if (order <= 90)
    return 90;
  else if (order <= 91)
    return 91;
  else if (order <= 94)
    return 94;
  else if (order <= 98)
    return 98;
  else if (order <= 99)
    return 99;
  else if (order <= 100)
    return 100;
  else if (order <= 101)
    return 101;
  else if (order <= 103)
    return 103;
  else if (order <= 105)
    return 105;
  else if (order <= 111)
    return 111;
  else if (order <= 113)
    return 113;
  else if (order <= 116)
    return 116;
  else if (order <= 117)
    return 117;
  else if (order <= 121)
    return 121;
  else if (order <= 124)
    return 124;
  else if (order <= 127)
    return 127;
  else if (order <= 128)
    return 128;
  else if (order <= 130)
    return 130;
  else if (order <= 135)
    return 135;
  else if (order <= 136)
    return 136;
  else if (order <= 143)
    return 143;
  else if (order <= 149)
    return 149;
  else
    return 210;
}
} // namespace detail::ahrensbeylkin

template <typename RealType>
class AhrensBeylkin : public Quadrature<AhrensBeylkin<RealType>> {

  using base_type = Quadrature<AhrensBeylkin<RealType>>;

public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  AhrensBeylkin(size_t npts) : base_type(npts) {}

  AhrensBeylkin(const AhrensBeylkin &) = default;
  AhrensBeylkin(AhrensBeylkin &&) noexcept = default;
};

template <typename RealType> struct quadrature_traits<AhrensBeylkin<RealType>> {

  using point_type = cartesian_pt_t<RealType>;
  using weight_type = RealType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container>
  generate(size_t npts) {

    point_container points(npts);
    weight_container weights(npts);

    using namespace AhrensBeylkinGrids;

    if (npts == 72)
      detail::copy_grid<ahrens_beylkin_72<RealType>>(points, weights);
    else if (npts == 132)
      detail::copy_grid<ahrens_beylkin_132<RealType>>(points, weights);
    else if (npts == 192)
      detail::copy_grid<ahrens_beylkin_192<RealType>>(points, weights);
    else if (npts == 312)
      detail::copy_grid<ahrens_beylkin_312<RealType>>(points, weights);
    else if (npts == 372)
      detail::copy_grid<ahrens_beylkin_372<RealType>>(points, weights);
    else if (npts == 432)
      detail::copy_grid<ahrens_beylkin_432<RealType>>(points, weights);
    else if (npts == 492)
      detail::copy_grid<ahrens_beylkin_492<RealType>>(points, weights);
    else if (npts == 552)
      detail::copy_grid<ahrens_beylkin_552<RealType>>(points, weights);
    else if (npts == 612)
      detail::copy_grid<ahrens_beylkin_612<RealType>>(points, weights);
    else if (npts == 672)
      detail::copy_grid<ahrens_beylkin_672<RealType>>(points, weights);
    else if (npts == 732)
      detail::copy_grid<ahrens_beylkin_732<RealType>>(points, weights);
    else if (npts == 792)
      detail::copy_grid<ahrens_beylkin_792<RealType>>(points, weights);
    else if (npts == 912)
      detail::copy_grid<ahrens_beylkin_912<RealType>>(points, weights);
    else if (npts == 972)
      detail::copy_grid<ahrens_beylkin_972<RealType>>(points, weights);
    else if (npts == 1032)
      detail::copy_grid<ahrens_beylkin_1032<RealType>>(points, weights);
    else if (npts == 1092)
      detail::copy_grid<ahrens_beylkin_1092<RealType>>(points, weights);
    else if (npts == 1272)
      detail::copy_grid<ahrens_beylkin_1272<RealType>>(points, weights);
    else if (npts == 1332)
      detail::copy_grid<ahrens_beylkin_1332<RealType>>(points, weights);
    else if (npts == 1392)
      detail::copy_grid<ahrens_beylkin_1392<RealType>>(points, weights);
    else if (npts == 1512)
      detail::copy_grid<ahrens_beylkin_1512<RealType>>(points, weights);
    else if (npts == 1572)
      detail::copy_grid<ahrens_beylkin_1572<RealType>>(points, weights);
    else if (npts == 1632)
      detail::copy_grid<ahrens_beylkin_1632<RealType>>(points, weights);
    else if (npts == 1692)
      detail::copy_grid<ahrens_beylkin_1692<RealType>>(points, weights);
    else if (npts == 1812)
      detail::copy_grid<ahrens_beylkin_1812<RealType>>(points, weights);
    else if (npts == 1932)
      detail::copy_grid<ahrens_beylkin_1932<RealType>>(points, weights);
    else if (npts == 1992)
      detail::copy_grid<ahrens_beylkin_1992<RealType>>(points, weights);
    else if (npts == 2112)
      detail::copy_grid<ahrens_beylkin_2112<RealType>>(points, weights);
    else if (npts == 2172)
      detail::copy_grid<ahrens_beylkin_2172<RealType>>(points, weights);
    else if (npts == 2232)
      detail::copy_grid<ahrens_beylkin_2232<RealType>>(points, weights);
    else if (npts == 2292)
      detail::copy_grid<ahrens_beylkin_2292<RealType>>(points, weights);
    else if (npts == 2412)
      detail::copy_grid<ahrens_beylkin_2412<RealType>>(points, weights);
    else if (npts == 2472)
      detail::copy_grid<ahrens_beylkin_2472<RealType>>(points, weights);
    else if (npts == 2712)
      detail::copy_grid<ahrens_beylkin_2712<RealType>>(points, weights);
    else if (npts == 2772)
      detail::copy_grid<ahrens_beylkin_2772<RealType>>(points, weights);
    else if (npts == 2832)
      detail::copy_grid<ahrens_beylkin_2832<RealType>>(points, weights);
    else if (npts == 3012)
      detail::copy_grid<ahrens_beylkin_3012<RealType>>(points, weights);
    else if (npts == 3312)
      detail::copy_grid<ahrens_beylkin_3312<RealType>>(points, weights);
    else if (npts == 3372)
      detail::copy_grid<ahrens_beylkin_3372<RealType>>(points, weights);
    else if (npts == 3432)
      detail::copy_grid<ahrens_beylkin_3432<RealType>>(points, weights);
    else if (npts == 3492)
      detail::copy_grid<ahrens_beylkin_3492<RealType>>(points, weights);
    else if (npts == 3612)
      detail::copy_grid<ahrens_beylkin_3612<RealType>>(points, weights);
    else if (npts == 3792)
      detail::copy_grid<ahrens_beylkin_3792<RealType>>(points, weights);
    else if (npts == 4212)
      detail::copy_grid<ahrens_beylkin_4212<RealType>>(points, weights);
    else if (npts == 4332)
      detail::copy_grid<ahrens_beylkin_4332<RealType>>(points, weights);
    else if (npts == 4572)
      detail::copy_grid<ahrens_beylkin_4572<RealType>>(points, weights);
    else if (npts == 4692)
      detail::copy_grid<ahrens_beylkin_4692<RealType>>(points, weights);
    else if (npts == 4992)
      detail::copy_grid<ahrens_beylkin_4992<RealType>>(points, weights);
    else if (npts == 5232)
      detail::copy_grid<ahrens_beylkin_5232<RealType>>(points, weights);
    else if (npts == 5472)
      detail::copy_grid<ahrens_beylkin_5472<RealType>>(points, weights);
    else if (npts == 5592)
      detail::copy_grid<ahrens_beylkin_5592<RealType>>(points, weights);
    else if (npts == 5772)
      detail::copy_grid<ahrens_beylkin_5772<RealType>>(points, weights);
    else if (npts == 6192)
      detail::copy_grid<ahrens_beylkin_6192<RealType>>(points, weights);
    else if (npts == 6312)
      detail::copy_grid<ahrens_beylkin_6312<RealType>>(points, weights);
    else if (npts == 6912)
      detail::copy_grid<ahrens_beylkin_6912<RealType>>(points, weights);
    else if (npts == 7512)
      detail::copy_grid<ahrens_beylkin_7512<RealType>>(points, weights);
    else if (npts == 15012)
      detail::copy_grid<ahrens_beylkin_15012<RealType>>(points, weights);

    return std::make_tuple(points, weights);
  }
};
} // namespace IntegratorXX
