#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/util/copy_grid.hpp>
#include <integratorxx/util/create_array.hpp>
#include <integratorxx/quadratures/s2/womersley/womersley_grids.hpp>
#include <vector>

namespace IntegratorXX {

/**
 *  @brief Implementation of Womersley's spherical t design quadrature on the
 * sphere.
 *
 *  References:
 *  Geometriae Dedicata (1977) 6, 363–388.
 *  doi:10.1007/bf03187604
 *
 *  Geometry and Combinatorics (1991) 68-93.
 *  doi:10.1016/B978-0-12-189420-7.50013-X
 *
 *  R. S. Womersley, Efficient Spherical Designs with Good Geometric
 *  Properties, in "Contemporary Computational Mathematics - A
 *  Celebration of the 80th Birthday of Ian Sloan", Springer
 *  2018. doi:10.1007/978-3-319-72456-0_57
 *
 *  @tparam PointContainer  Quadrature points
 *  @tparam WeightContainer Quadrature weights
 */
template <typename RealType>
class Womersley : public Quadrature<Womersley<RealType>> {
  using base_type = Quadrature<Womersley<RealType>>;

 public:
  using point_type = typename base_type::point_type;
  using weight_type = typename base_type::weight_type;
  using point_container = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  Womersley(size_t npts) : base_type(npts) {}

  Womersley(const Womersley &) = default;
  Womersley(Womersley &&) noexcept = default;
};

template <typename RealType>
struct quadrature_traits<Womersley<RealType>> {
  using point_type = cartesian_pt_t<RealType>;
  using weight_type = RealType;

  using point_container = std::vector<point_type>;
  using weight_container = std::vector<weight_type>;

  inline static std::tuple<point_container, weight_container> generate(
      size_t npts) {
    point_container points(npts);
    weight_container weights(npts);

    using namespace WomersleyGrids;

    if(npts == 3)
      detail::copy_grid<womersley_3<RealType>>(points, weights);
    else if(npts == 6)
      detail::copy_grid<womersley_6<RealType>>(points, weights);
    else if(npts == 8)
      detail::copy_grid<womersley_8<RealType>>(points, weights);
    else if(npts == 14)
      detail::copy_grid<womersley_14<RealType>>(points, weights);
    else if(npts == 18)
      detail::copy_grid<womersley_18<RealType>>(points, weights);
    else if(npts == 26)
      detail::copy_grid<womersley_26<RealType>>(points, weights);
    else if(npts == 32)
      detail::copy_grid<womersley_32<RealType>>(points, weights);
    else if(npts == 42)
      detail::copy_grid<womersley_42<RealType>>(points, weights);
    else if(npts == 50)
      detail::copy_grid<womersley_50<RealType>>(points, weights);
    else if(npts == 62)
      detail::copy_grid<womersley_62<RealType>>(points, weights);
    else if(npts == 72)
      detail::copy_grid<womersley_72<RealType>>(points, weights);
    else if(npts == 86)
      detail::copy_grid<womersley_86<RealType>>(points, weights);
    else if(npts == 98)
      detail::copy_grid<womersley_98<RealType>>(points, weights);
    else if(npts == 114)
      detail::copy_grid<womersley_114<RealType>>(points, weights);
    else if(npts == 128)
      detail::copy_grid<womersley_128<RealType>>(points, weights);
    else if(npts == 146)
      detail::copy_grid<womersley_146<RealType>>(points, weights);
    else if(npts == 163)
      detail::copy_grid<womersley_163<RealType>>(points, weights);
    else if(npts == 182)
      detail::copy_grid<womersley_182<RealType>>(points, weights);
    else if(npts == 201)
      detail::copy_grid<womersley_201<RealType>>(points, weights);
    else if(npts == 222)
      detail::copy_grid<womersley_222<RealType>>(points, weights);
    else if(npts == 243)
      detail::copy_grid<womersley_243<RealType>>(points, weights);
    else if(npts == 266)
      detail::copy_grid<womersley_266<RealType>>(points, weights);
    else if(npts == 289)
      detail::copy_grid<womersley_289<RealType>>(points, weights);
    else if(npts == 314)
      detail::copy_grid<womersley_314<RealType>>(points, weights);
    else if(npts == 339)
      detail::copy_grid<womersley_339<RealType>>(points, weights);
    else if(npts == 366)
      detail::copy_grid<womersley_366<RealType>>(points, weights);
    else if(npts == 393)
      detail::copy_grid<womersley_393<RealType>>(points, weights);
    else if(npts == 422)
      detail::copy_grid<womersley_422<RealType>>(points, weights);
    else if(npts == 451)
      detail::copy_grid<womersley_451<RealType>>(points, weights);
    else if(npts == 482)
      detail::copy_grid<womersley_482<RealType>>(points, weights);
    else if(npts == 513)
      detail::copy_grid<womersley_513<RealType>>(points, weights);
    else if(npts == 546)
      detail::copy_grid<womersley_546<RealType>>(points, weights);
    else if(npts == 579)
      detail::copy_grid<womersley_579<RealType>>(points, weights);
    else if(npts == 614)
      detail::copy_grid<womersley_614<RealType>>(points, weights);
    else if(npts == 649)
      detail::copy_grid<womersley_649<RealType>>(points, weights);
    else if(npts == 686)
      detail::copy_grid<womersley_686<RealType>>(points, weights);
    else if(npts == 723)
      detail::copy_grid<womersley_723<RealType>>(points, weights);
    else if(npts == 762)
      detail::copy_grid<womersley_762<RealType>>(points, weights);
    else if(npts == 801)
      detail::copy_grid<womersley_801<RealType>>(points, weights);
    else if(npts == 842)
      detail::copy_grid<womersley_842<RealType>>(points, weights);
    else if(npts == 883)
      detail::copy_grid<womersley_883<RealType>>(points, weights);
    else if(npts == 926)
      detail::copy_grid<womersley_926<RealType>>(points, weights);
    else if(npts == 969)
      detail::copy_grid<womersley_969<RealType>>(points, weights);
    else if(npts == 1014)
      detail::copy_grid<womersley_1014<RealType>>(points, weights);
    else if(npts == 1059)
      detail::copy_grid<womersley_1059<RealType>>(points, weights);
    else if(npts == 1106)
      detail::copy_grid<womersley_1106<RealType>>(points, weights);
    else if(npts == 1153)
      detail::copy_grid<womersley_1153<RealType>>(points, weights);
    else if(npts == 1202)
      detail::copy_grid<womersley_1202<RealType>>(points, weights);
    else if(npts == 1251)
      detail::copy_grid<womersley_1251<RealType>>(points, weights);
    else if(npts == 1302)
      detail::copy_grid<womersley_1302<RealType>>(points, weights);
    else if(npts == 1353)
      detail::copy_grid<womersley_1353<RealType>>(points, weights);
    else if(npts == 1406)
      detail::copy_grid<womersley_1406<RealType>>(points, weights);
    else if(npts == 1459)
      detail::copy_grid<womersley_1459<RealType>>(points, weights);
    else if(npts == 1514)
      detail::copy_grid<womersley_1514<RealType>>(points, weights);
    else if(npts == 1569)
      detail::copy_grid<womersley_1569<RealType>>(points, weights);
    else if(npts == 1626)
      detail::copy_grid<womersley_1626<RealType>>(points, weights);
    else if(npts == 1683)
      detail::copy_grid<womersley_1683<RealType>>(points, weights);
    else if(npts == 1742)
      detail::copy_grid<womersley_1742<RealType>>(points, weights);
    else if(npts == 1801)
      detail::copy_grid<womersley_1801<RealType>>(points, weights);
    else if(npts == 1862)
      detail::copy_grid<womersley_1862<RealType>>(points, weights);
    else if(npts == 1923)
      detail::copy_grid<womersley_1923<RealType>>(points, weights);
    else if(npts == 1986)
      detail::copy_grid<womersley_1986<RealType>>(points, weights);
    else if(npts == 2049)
      detail::copy_grid<womersley_2049<RealType>>(points, weights);
    else if(npts == 2114)
      detail::copy_grid<womersley_2114<RealType>>(points, weights);
    else if(npts == 2179)
      detail::copy_grid<womersley_2179<RealType>>(points, weights);
    else if(npts == 2246)
      detail::copy_grid<womersley_2246<RealType>>(points, weights);
    else if(npts == 2313)
      detail::copy_grid<womersley_2313<RealType>>(points, weights);
    else if(npts == 2382)
      detail::copy_grid<womersley_2382<RealType>>(points, weights);
    else if(npts == 2451)
      detail::copy_grid<womersley_2451<RealType>>(points, weights);
    else if(npts == 2522)
      detail::copy_grid<womersley_2522<RealType>>(points, weights);
    else if(npts == 2593)
      detail::copy_grid<womersley_2593<RealType>>(points, weights);
    else if(npts == 2666)
      detail::copy_grid<womersley_2666<RealType>>(points, weights);
    else if(npts == 2739)
      detail::copy_grid<womersley_2739<RealType>>(points, weights);
    else if(npts == 2814)
      detail::copy_grid<womersley_2814<RealType>>(points, weights);
    else if(npts == 2889)
      detail::copy_grid<womersley_2889<RealType>>(points, weights);
    else if(npts == 2966)
      detail::copy_grid<womersley_2966<RealType>>(points, weights);
    else if(npts == 3043)
      detail::copy_grid<womersley_3043<RealType>>(points, weights);
    else if(npts == 3122)
      detail::copy_grid<womersley_3122<RealType>>(points, weights);
    else if(npts == 3201)
      detail::copy_grid<womersley_3201<RealType>>(points, weights);
    else if(npts == 3282)
      detail::copy_grid<womersley_3282<RealType>>(points, weights);
    else if(npts == 3363)
      detail::copy_grid<womersley_3363<RealType>>(points, weights);
    else if(npts == 3446)
      detail::copy_grid<womersley_3446<RealType>>(points, weights);
    else if(npts == 3529)
      detail::copy_grid<womersley_3529<RealType>>(points, weights);
    else if(npts == 3614)
      detail::copy_grid<womersley_3614<RealType>>(points, weights);
    else if(npts == 3699)
      detail::copy_grid<womersley_3699<RealType>>(points, weights);
    else if(npts == 3786)
      detail::copy_grid<womersley_3786<RealType>>(points, weights);
    else if(npts == 3873)
      detail::copy_grid<womersley_3873<RealType>>(points, weights);
    else if(npts == 3962)
      detail::copy_grid<womersley_3962<RealType>>(points, weights);
    else if(npts == 4051)
      detail::copy_grid<womersley_4051<RealType>>(points, weights);
    else if(npts == 4142)
      detail::copy_grid<womersley_4142<RealType>>(points, weights);
    else if(npts == 4233)
      detail::copy_grid<womersley_4233<RealType>>(points, weights);
    else if(npts == 4326)
      detail::copy_grid<womersley_4326<RealType>>(points, weights);
    else if(npts == 4419)
      detail::copy_grid<womersley_4419<RealType>>(points, weights);
    else if(npts == 4514)
      detail::copy_grid<womersley_4514<RealType>>(points, weights);
    else if(npts == 4609)
      detail::copy_grid<womersley_4609<RealType>>(points, weights);
    else if(npts == 4706)
      detail::copy_grid<womersley_4706<RealType>>(points, weights);
    else if(npts == 4803)
      detail::copy_grid<womersley_4803<RealType>>(points, weights);
    else if(npts == 4902)
      detail::copy_grid<womersley_4902<RealType>>(points, weights);
    else if(npts == 5001)
      detail::copy_grid<womersley_5001<RealType>>(points, weights);
    else if(npts == 5102)
      detail::copy_grid<womersley_5102<RealType>>(points, weights);
    else if(npts == 5203)
      detail::copy_grid<womersley_5203<RealType>>(points, weights);
    else if(npts == 5306)
      detail::copy_grid<womersley_5306<RealType>>(points, weights);
    else if(npts == 5409)
      detail::copy_grid<womersley_5409<RealType>>(points, weights);
    else if(npts == 5514)
      detail::copy_grid<womersley_5514<RealType>>(points, weights);
    else if(npts == 5619)
      detail::copy_grid<womersley_5619<RealType>>(points, weights);
    else if(npts == 5726)
      detail::copy_grid<womersley_5726<RealType>>(points, weights);
    else if(npts == 5833)
      detail::copy_grid<womersley_5833<RealType>>(points, weights);
    else if(npts == 5942)
      detail::copy_grid<womersley_5942<RealType>>(points, weights);
    else if(npts == 6051)
      detail::copy_grid<womersley_6051<RealType>>(points, weights);
    else if(npts == 6162)
      detail::copy_grid<womersley_6162<RealType>>(points, weights);
    else if(npts == 6273)
      detail::copy_grid<womersley_6273<RealType>>(points, weights);
    else if(npts == 6386)
      detail::copy_grid<womersley_6386<RealType>>(points, weights);
    else if(npts == 6499)
      detail::copy_grid<womersley_6499<RealType>>(points, weights);
    else if(npts == 6614)
      detail::copy_grid<womersley_6614<RealType>>(points, weights);
    else if(npts == 6729)
      detail::copy_grid<womersley_6729<RealType>>(points, weights);
    else if(npts == 6846)
      detail::copy_grid<womersley_6846<RealType>>(points, weights);
    else if(npts == 6963)
      detail::copy_grid<womersley_6963<RealType>>(points, weights);
    else if(npts == 7082)
      detail::copy_grid<womersley_7082<RealType>>(points, weights);
    else if(npts == 7201)
      detail::copy_grid<womersley_7201<RealType>>(points, weights);
    else if(npts == 7322)
      detail::copy_grid<womersley_7322<RealType>>(points, weights);
    else if(npts == 7443)
      detail::copy_grid<womersley_7443<RealType>>(points, weights);
    else if(npts == 7566)
      detail::copy_grid<womersley_7566<RealType>>(points, weights);
    else if(npts == 7689)
      detail::copy_grid<womersley_7689<RealType>>(points, weights);
    else if(npts == 7814)
      detail::copy_grid<womersley_7814<RealType>>(points, weights);
    else if(npts == 7939)
      detail::copy_grid<womersley_7939<RealType>>(points, weights);
    return std::make_tuple(points, weights);
  }

inline static int64_t npts_by_algebraic_order(int64_t order) {
  switch(order) {
    case 1:
      return 3;
    case 2:
      return 6;
    case 3:
      return 8;
    case 4:
      return 14;
    case 5:
      return 18;
    case 6:
      return 26;
    case 7:
      return 32;
    case 8:
      return 42;
    case 9:
      return 50;
    case 10:
      return 62;
    case 11:
      return 72;
    case 12:
      return 86;
    case 13:
      return 98;
    case 14:
      return 114;
    case 15:
      return 128;
    case 16:
      return 146;
    case 17:
      return 163;
    case 18:
      return 182;
    case 19:
      return 201;
    case 20:
      return 222;
    case 21:
      return 243;
    case 22:
      return 266;
    case 23:
      return 289;
    case 24:
      return 314;
    case 25:
      return 339;
    case 26:
      return 366;
    case 27:
      return 393;
    case 28:
      return 422;
    case 29:
      return 451;
    case 30:
      return 482;
    case 31:
      return 513;
    case 32:
      return 546;
    case 33:
      return 579;
    case 34:
      return 614;
    case 35:
      return 649;
    case 36:
      return 686;
    case 37:
      return 723;
    case 38:
      return 762;
    case 39:
      return 801;
    case 40:
      return 842;
    case 41:
      return 883;
    case 42:
      return 926;
    case 43:
      return 969;
    case 44:
      return 1014;
    case 45:
      return 1059;
    case 46:
      return 1106;
    case 47:
      return 1153;
    case 48:
      return 1202;
    case 49:
      return 1251;
    case 50:
      return 1302;
    case 51:
      return 1353;
    case 52:
      return 1406;
    case 53:
      return 1459;
    case 54:
      return 1514;
    case 55:
      return 1569;
    case 56:
      return 1626;
    case 57:
      return 1683;
    case 58:
      return 1742;
    case 59:
      return 1801;
    case 60:
      return 1862;
    case 61:
      return 1923;
    case 62:
      return 1986;
    case 63:
      return 2049;
    case 64:
      return 2114;
    case 65:
      return 2179;
    case 66:
      return 2246;
    case 67:
      return 2313;
    case 68:
      return 2382;
    case 69:
      return 2451;
    case 70:
      return 2522;
    case 71:
      return 2593;
    case 72:
      return 2666;
    case 73:
      return 2739;
    case 74:
      return 2814;
    case 75:
      return 2889;
    case 76:
      return 2966;
    case 77:
      return 3043;
    case 78:
      return 3122;
    case 79:
      return 3201;
    case 80:
      return 3282;
    case 81:
      return 3363;
    case 82:
      return 3446;
    case 83:
      return 3529;
    case 84:
      return 3614;
    case 85:
      return 3699;
    case 86:
      return 3786;
    case 87:
      return 3873;
    case 88:
      return 3962;
    case 89:
      return 4051;
    case 90:
      return 4142;
    case 91:
      return 4233;
    case 92:
      return 4326;
    case 93:
      return 4419;
    case 94:
      return 4514;
    case 95:
      return 4609;
    case 96:
      return 4706;
    case 97:
      return 4803;
    case 98:
      return 4902;
    case 99:
      return 5001;
    case 100:
      return 5102;
    case 101:
      return 5203;
    case 102:
      return 5306;
    case 103:
      return 5409;
    case 104:
      return 5514;
    case 105:
      return 5619;
    case 106:
      return 5726;
    case 107:
      return 5833;
    case 108:
      return 5942;
    case 109:
      return 6051;
    case 110:
      return 6162;
    case 111:
      return 6273;
    case 112:
      return 6386;
    case 113:
      return 6499;
    case 114:
      return 6614;
    case 115:
      return 6729;
    case 116:
      return 6846;
    case 117:
      return 6963;
    case 118:
      return 7082;
    case 119:
      return 7201;
    case 120:
      return 7322;
    case 121:
      return 7443;
    case 122:
      return 7566;
    case 123:
      return 7689;
    case 124:
      return 7814;
    case 125:
      return 7939;
    default:
      return -1;
  }
}

inline static int64_t algebraic_order_by_npts(int64_t npts) {
  switch(npts) {
    case 3:
      return 1;
    case 6:
      return 2;
    case 8:
      return 3;
    case 14:
      return 4;
    case 18:
      return 5;
    case 26:
      return 6;
    case 32:
      return 7;
    case 42:
      return 8;
    case 50:
      return 9;
    case 62:
      return 10;
    case 72:
      return 11;
    case 86:
      return 12;
    case 98:
      return 13;
    case 114:
      return 14;
    case 128:
      return 15;
    case 146:
      return 16;
    case 163:
      return 17;
    case 182:
      return 18;
    case 201:
      return 19;
    case 222:
      return 20;
    case 243:
      return 21;
    case 266:
      return 22;
    case 289:
      return 23;
    case 314:
      return 24;
    case 339:
      return 25;
    case 366:
      return 26;
    case 393:
      return 27;
    case 422:
      return 28;
    case 451:
      return 29;
    case 482:
      return 30;
    case 513:
      return 31;
    case 546:
      return 32;
    case 579:
      return 33;
    case 614:
      return 34;
    case 649:
      return 35;
    case 686:
      return 36;
    case 723:
      return 37;
    case 762:
      return 38;
    case 801:
      return 39;
    case 842:
      return 40;
    case 883:
      return 41;
    case 926:
      return 42;
    case 969:
      return 43;
    case 1014:
      return 44;
    case 1059:
      return 45;
    case 1106:
      return 46;
    case 1153:
      return 47;
    case 1202:
      return 48;
    case 1251:
      return 49;
    case 1302:
      return 50;
    case 1353:
      return 51;
    case 1406:
      return 52;
    case 1459:
      return 53;
    case 1514:
      return 54;
    case 1569:
      return 55;
    case 1626:
      return 56;
    case 1683:
      return 57;
    case 1742:
      return 58;
    case 1801:
      return 59;
    case 1862:
      return 60;
    case 1923:
      return 61;
    case 1986:
      return 62;
    case 2049:
      return 63;
    case 2114:
      return 64;
    case 2179:
      return 65;
    case 2246:
      return 66;
    case 2313:
      return 67;
    case 2382:
      return 68;
    case 2451:
      return 69;
    case 2522:
      return 70;
    case 2593:
      return 71;
    case 2666:
      return 72;
    case 2739:
      return 73;
    case 2814:
      return 74;
    case 2889:
      return 75;
    case 2966:
      return 76;
    case 3043:
      return 77;
    case 3122:
      return 78;
    case 3201:
      return 79;
    case 3282:
      return 80;
    case 3363:
      return 81;
    case 3446:
      return 82;
    case 3529:
      return 83;
    case 3614:
      return 84;
    case 3699:
      return 85;
    case 3786:
      return 86;
    case 3873:
      return 87;
    case 3962:
      return 88;
    case 4051:
      return 89;
    case 4142:
      return 90;
    case 4233:
      return 91;
    case 4326:
      return 92;
    case 4419:
      return 93;
    case 4514:
      return 94;
    case 4609:
      return 95;
    case 4706:
      return 96;
    case 4803:
      return 97;
    case 4902:
      return 98;
    case 5001:
      return 99;
    case 5102:
      return 100;
    case 5203:
      return 101;
    case 5306:
      return 102;
    case 5409:
      return 103;
    case 5514:
      return 104;
    case 5619:
      return 105;
    case 5726:
      return 106;
    case 5833:
      return 107;
    case 5942:
      return 108;
    case 6051:
      return 109;
    case 6162:
      return 110;
    case 6273:
      return 111;
    case 6386:
      return 112;
    case 6499:
      return 113;
    case 6614:
      return 114;
    case 6729:
      return 115;
    case 6846:
      return 116;
    case 6963:
      return 117;
    case 7082:
      return 118;
    case 7201:
      return 119;
    case 7322:
      return 120;
    case 7443:
      return 121;
    case 7566:
      return 122;
    case 7689:
      return 123;
    case 7814:
      return 124;
    case 7939:
      return 125;
    default:
      return -1;
  }
}

inline static int64_t next_algebraic_order(int64_t order) {
  if(order <= 1)
    return 1;
  else if(order <= 2)
    return 2;
  else if(order <= 3)
    return 3;
  else if(order <= 4)
    return 4;
  else if(order <= 5)
    return 5;
  else if(order <= 6)
    return 6;
  else if(order <= 7)
    return 7;
  else if(order <= 8)
    return 8;
  else if(order <= 9)
    return 9;
  else if(order <= 10)
    return 10;
  else if(order <= 11)
    return 11;
  else if(order <= 12)
    return 12;
  else if(order <= 13)
    return 13;
  else if(order <= 14)
    return 14;
  else if(order <= 15)
    return 15;
  else if(order <= 16)
    return 16;
  else if(order <= 17)
    return 17;
  else if(order <= 18)
    return 18;
  else if(order <= 19)
    return 19;
  else if(order <= 20)
    return 20;
  else if(order <= 21)
    return 21;
  else if(order <= 22)
    return 22;
  else if(order <= 23)
    return 23;
  else if(order <= 24)
    return 24;
  else if(order <= 25)
    return 25;
  else if(order <= 26)
    return 26;
  else if(order <= 27)
    return 27;
  else if(order <= 28)
    return 28;
  else if(order <= 29)
    return 29;
  else if(order <= 30)
    return 30;
  else if(order <= 31)
    return 31;
  else if(order <= 32)
    return 32;
  else if(order <= 33)
    return 33;
  else if(order <= 34)
    return 34;
  else if(order <= 35)
    return 35;
  else if(order <= 36)
    return 36;
  else if(order <= 37)
    return 37;
  else if(order <= 38)
    return 38;
  else if(order <= 39)
    return 39;
  else if(order <= 40)
    return 40;
  else if(order <= 41)
    return 41;
  else if(order <= 42)
    return 42;
  else if(order <= 43)
    return 43;
  else if(order <= 44)
    return 44;
  else if(order <= 45)
    return 45;
  else if(order <= 46)
    return 46;
  else if(order <= 47)
    return 47;
  else if(order <= 48)
    return 48;
  else if(order <= 49)
    return 49;
  else if(order <= 50)
    return 50;
  else if(order <= 51)
    return 51;
  else if(order <= 52)
    return 52;
  else if(order <= 53)
    return 53;
  else if(order <= 54)
    return 54;
  else if(order <= 55)
    return 55;
  else if(order <= 56)
    return 56;
  else if(order <= 57)
    return 57;
  else if(order <= 58)
    return 58;
  else if(order <= 59)
    return 59;
  else if(order <= 60)
    return 60;
  else if(order <= 61)
    return 61;
  else if(order <= 62)
    return 62;
  else if(order <= 63)
    return 63;
  else if(order <= 64)
    return 64;
  else if(order <= 65)
    return 65;
  else if(order <= 66)
    return 66;
  else if(order <= 67)
    return 67;
  else if(order <= 68)
    return 68;
  else if(order <= 69)
    return 69;
  else if(order <= 70)
    return 70;
  else if(order <= 71)
    return 71;
  else if(order <= 72)
    return 72;
  else if(order <= 73)
    return 73;
  else if(order <= 74)
    return 74;
  else if(order <= 75)
    return 75;
  else if(order <= 76)
    return 76;
  else if(order <= 77)
    return 77;
  else if(order <= 78)
    return 78;
  else if(order <= 79)
    return 79;
  else if(order <= 80)
    return 80;
  else if(order <= 81)
    return 81;
  else if(order <= 82)
    return 82;
  else if(order <= 83)
    return 83;
  else if(order <= 84)
    return 84;
  else if(order <= 85)
    return 85;
  else if(order <= 86)
    return 86;
  else if(order <= 87)
    return 87;
  else if(order <= 88)
    return 88;
  else if(order <= 89)
    return 89;
  else if(order <= 90)
    return 90;
  else if(order <= 91)
    return 91;
  else if(order <= 92)
    return 92;
  else if(order <= 93)
    return 93;
  else if(order <= 94)
    return 94;
  else if(order <= 95)
    return 95;
  else if(order <= 96)
    return 96;
  else if(order <= 97)
    return 97;
  else if(order <= 98)
    return 98;
  else if(order <= 99)
    return 99;
  else if(order <= 100)
    return 100;
  else if(order <= 101)
    return 101;
  else if(order <= 102)
    return 102;
  else if(order <= 103)
    return 103;
  else if(order <= 104)
    return 104;
  else if(order <= 105)
    return 105;
  else if(order <= 106)
    return 106;
  else if(order <= 107)
    return 107;
  else if(order <= 108)
    return 108;
  else if(order <= 109)
    return 109;
  else if(order <= 110)
    return 110;
  else if(order <= 111)
    return 111;
  else if(order <= 112)
    return 112;
  else if(order <= 113)
    return 113;
  else if(order <= 114)
    return 114;
  else if(order <= 115)
    return 115;
  else if(order <= 116)
    return 116;
  else if(order <= 117)
    return 117;
  else if(order <= 118)
    return 118;
  else if(order <= 119)
    return 119;
  else if(order <= 120)
    return 120;
  else if(order <= 121)
    return 121;
  else if(order <= 122)
    return 122;
  else if(order <= 123)
    return 123;
  else if(order <= 124)
    return 124;
  else
    return 125;
}
};
namespace detail {

template <typename QuadType>
static constexpr bool is_womersley_v = std::is_same_v<
  QuadType, 
  Womersley<typename QuadType::weight_type>
>;

}
}  // namespace IntegratorXX
