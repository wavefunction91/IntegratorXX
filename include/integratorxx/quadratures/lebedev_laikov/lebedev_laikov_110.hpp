#pragma once

namespace IntegratorXX {
namespace LebedevLaikovGrids {

/**
 *  \brief Lebedev-Laikov Quadrature specification for Order = 110
 * 
 */
template <typename T>
struct lebedev_laikov_110 {

  static constexpr std::array<cartesian_pt_t<T>,110> points = {
      1.000000000000000e+00,      0.000000000000000e+00,      0.000000000000000e+00,
     -1.000000000000000e+00,      0.000000000000000e+00,      0.000000000000000e+00,
      0.000000000000000e+00,      1.000000000000000e+00,      0.000000000000000e+00,
      0.000000000000000e+00,     -1.000000000000000e+00,      0.000000000000000e+00,
      0.000000000000000e+00,      0.000000000000000e+00,      1.000000000000000e+00,
      0.000000000000000e+00,      0.000000000000000e+00,     -1.000000000000000e+00,
      5.773502691896257e-01,      5.773502691896257e-01,      5.773502691896257e-01,
     -5.773502691896257e-01,      5.773502691896257e-01,      5.773502691896257e-01,
      5.773502691896257e-01,     -5.773502691896257e-01,      5.773502691896257e-01,
     -5.773502691896257e-01,     -5.773502691896257e-01,      5.773502691896257e-01,
      5.773502691896257e-01,      5.773502691896257e-01,     -5.773502691896257e-01,
     -5.773502691896257e-01,      5.773502691896257e-01,     -5.773502691896257e-01,
      5.773502691896257e-01,     -5.773502691896257e-01,     -5.773502691896257e-01,
     -5.773502691896257e-01,     -5.773502691896257e-01,     -5.773502691896257e-01,
      1.851156353447362e-01,      1.851156353447362e-01,      9.651240350865941e-01,
     -1.851156353447362e-01,      1.851156353447362e-01,      9.651240350865941e-01,
      1.851156353447362e-01,     -1.851156353447362e-01,      9.651240350865941e-01,
     -1.851156353447362e-01,     -1.851156353447362e-01,      9.651240350865941e-01,
      1.851156353447362e-01,      1.851156353447362e-01,     -9.651240350865941e-01,
     -1.851156353447362e-01,      1.851156353447362e-01,     -9.651240350865941e-01,
      1.851156353447362e-01,     -1.851156353447362e-01,     -9.651240350865941e-01,
     -1.851156353447362e-01,     -1.851156353447362e-01,     -9.651240350865941e-01,
      1.851156353447362e-01,      9.651240350865941e-01,      1.851156353447362e-01,
     -1.851156353447362e-01,      9.651240350865941e-01,      1.851156353447362e-01,
      1.851156353447362e-01,     -9.651240350865941e-01,      1.851156353447362e-01,
     -1.851156353447362e-01,     -9.651240350865941e-01,      1.851156353447362e-01,
      1.851156353447362e-01,      9.651240350865941e-01,     -1.851156353447362e-01,
     -1.851156353447362e-01,      9.651240350865941e-01,     -1.851156353447362e-01,
      1.851156353447362e-01,     -9.651240350865941e-01,     -1.851156353447362e-01,
     -1.851156353447362e-01,     -9.651240350865941e-01,     -1.851156353447362e-01,
      9.651240350865941e-01,      1.851156353447362e-01,      1.851156353447362e-01,
     -9.651240350865941e-01,      1.851156353447362e-01,      1.851156353447362e-01,
      9.651240350865941e-01,     -1.851156353447362e-01,      1.851156353447362e-01,
     -9.651240350865941e-01,     -1.851156353447362e-01,      1.851156353447362e-01,
      9.651240350865941e-01,      1.851156353447362e-01,     -1.851156353447362e-01,
     -9.651240350865941e-01,      1.851156353447362e-01,     -1.851156353447362e-01,
      9.651240350865941e-01,     -1.851156353447362e-01,     -1.851156353447362e-01,
     -9.651240350865941e-01,     -1.851156353447362e-01,     -1.851156353447362e-01,
      6.904210483822922e-01,      6.904210483822922e-01,      2.159572918458484e-01,
     -6.904210483822922e-01,      6.904210483822922e-01,      2.159572918458484e-01,
      6.904210483822922e-01,     -6.904210483822922e-01,      2.159572918458484e-01,
     -6.904210483822922e-01,     -6.904210483822922e-01,      2.159572918458484e-01,
      6.904210483822922e-01,      6.904210483822922e-01,     -2.159572918458484e-01,
     -6.904210483822922e-01,      6.904210483822922e-01,     -2.159572918458484e-01,
      6.904210483822922e-01,     -6.904210483822922e-01,     -2.159572918458484e-01,
     -6.904210483822922e-01,     -6.904210483822922e-01,     -2.159572918458484e-01,
      6.904210483822922e-01,      2.159572918458484e-01,      6.904210483822922e-01,
     -6.904210483822922e-01,      2.159572918458484e-01,      6.904210483822922e-01,
      6.904210483822922e-01,     -2.159572918458484e-01,      6.904210483822922e-01,
     -6.904210483822922e-01,     -2.159572918458484e-01,      6.904210483822922e-01,
      6.904210483822922e-01,      2.159572918458484e-01,     -6.904210483822922e-01,
     -6.904210483822922e-01,      2.159572918458484e-01,     -6.904210483822922e-01,
      6.904210483822922e-01,     -2.159572918458484e-01,     -6.904210483822922e-01,
     -6.904210483822922e-01,     -2.159572918458484e-01,     -6.904210483822922e-01,
      2.159572918458484e-01,      6.904210483822922e-01,      6.904210483822922e-01,
     -2.159572918458484e-01,      6.904210483822922e-01,      6.904210483822922e-01,
      2.159572918458484e-01,     -6.904210483822922e-01,      6.904210483822922e-01,
     -2.159572918458484e-01,     -6.904210483822922e-01,      6.904210483822922e-01,
      2.159572918458484e-01,      6.904210483822922e-01,     -6.904210483822922e-01,
     -2.159572918458484e-01,      6.904210483822922e-01,     -6.904210483822922e-01,
      2.159572918458484e-01,     -6.904210483822922e-01,     -6.904210483822922e-01,
     -2.159572918458484e-01,     -6.904210483822922e-01,     -6.904210483822922e-01,
      3.956894730559419e-01,      3.956894730559419e-01,      8.287699812525923e-01,
     -3.956894730559419e-01,      3.956894730559419e-01,      8.287699812525923e-01,
      3.956894730559419e-01,     -3.956894730559419e-01,      8.287699812525923e-01,
     -3.956894730559419e-01,     -3.956894730559419e-01,      8.287699812525923e-01,
      3.956894730559419e-01,      3.956894730559419e-01,     -8.287699812525923e-01,
     -3.956894730559419e-01,      3.956894730559419e-01,     -8.287699812525923e-01,
      3.956894730559419e-01,     -3.956894730559419e-01,     -8.287699812525923e-01,
     -3.956894730559419e-01,     -3.956894730559419e-01,     -8.287699812525923e-01,
      3.956894730559419e-01,      8.287699812525923e-01,      3.956894730559419e-01,
     -3.956894730559419e-01,      8.287699812525923e-01,      3.956894730559419e-01,
      3.956894730559419e-01,     -8.287699812525923e-01,      3.956894730559419e-01,
     -3.956894730559419e-01,     -8.287699812525923e-01,      3.956894730559419e-01,
      3.956894730559419e-01,      8.287699812525923e-01,     -3.956894730559419e-01,
     -3.956894730559419e-01,      8.287699812525923e-01,     -3.956894730559419e-01,
      3.956894730559419e-01,     -8.287699812525923e-01,     -3.956894730559419e-01,
     -3.956894730559419e-01,     -8.287699812525923e-01,     -3.956894730559419e-01,
      8.287699812525923e-01,      3.956894730559419e-01,      3.956894730559419e-01,
     -8.287699812525923e-01,      3.956894730559419e-01,      3.956894730559419e-01,
      8.287699812525923e-01,     -3.956894730559419e-01,      3.956894730559419e-01,
     -8.287699812525923e-01,     -3.956894730559419e-01,      3.956894730559419e-01,
      8.287699812525923e-01,      3.956894730559419e-01,     -3.956894730559419e-01,
     -8.287699812525923e-01,      3.956894730559419e-01,     -3.956894730559419e-01,
      8.287699812525923e-01,     -3.956894730559419e-01,     -3.956894730559419e-01,
     -8.287699812525923e-01,     -3.956894730559419e-01,     -3.956894730559419e-01,
      4.783690288121502e-01,      8.781589106040661e-01,      0.000000000000000e+00,
     -4.783690288121502e-01,      8.781589106040661e-01,      0.000000000000000e+00,
      4.783690288121502e-01,     -8.781589106040661e-01,      0.000000000000000e+00,
     -4.783690288121502e-01,     -8.781589106040661e-01,      0.000000000000000e+00,
      8.781589106040661e-01,      4.783690288121502e-01,      0.000000000000000e+00,
     -8.781589106040661e-01,      4.783690288121502e-01,      0.000000000000000e+00,
      8.781589106040661e-01,     -4.783690288121502e-01,      0.000000000000000e+00,
     -8.781589106040661e-01,     -4.783690288121502e-01,      0.000000000000000e+00,
      4.783690288121502e-01,      0.000000000000000e+00,      8.781589106040661e-01,
     -4.783690288121502e-01,      0.000000000000000e+00,      8.781589106040661e-01,
      4.783690288121502e-01,      0.000000000000000e+00,     -8.781589106040661e-01,
     -4.783690288121502e-01,      0.000000000000000e+00,     -8.781589106040661e-01,
      8.781589106040661e-01,      0.000000000000000e+00,      4.783690288121502e-01,
     -8.781589106040661e-01,      0.000000000000000e+00,      4.783690288121502e-01,
      8.781589106040661e-01,      0.000000000000000e+00,     -4.783690288121502e-01,
     -8.781589106040661e-01,      0.000000000000000e+00,     -4.783690288121502e-01,
      0.000000000000000e+00,      4.783690288121502e-01,      8.781589106040661e-01,
      0.000000000000000e+00,     -4.783690288121502e-01,      8.781589106040661e-01,
      0.000000000000000e+00,      4.783690288121502e-01,     -8.781589106040661e-01,
      0.000000000000000e+00,     -4.783690288121502e-01,     -8.781589106040661e-01,
      0.000000000000000e+00,      8.781589106040661e-01,      4.783690288121502e-01,
      0.000000000000000e+00,     -8.781589106040661e-01,      4.783690288121502e-01,
      0.000000000000000e+00,      8.781589106040661e-01,     -4.783690288121502e-01,
      0.000000000000000e+00,     -8.781589106040661e-01,     -4.783690288121502e-01
  };


  static constexpr std::array<T,110> weights = {
        3.828270494937162e-03,
        3.828270494937162e-03,
        3.828270494937162e-03,
        3.828270494937162e-03,
        3.828270494937162e-03,
        3.828270494937162e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        9.793737512487513e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        8.211737283191111e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.942814891178103e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.595471336070962e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03,
        9.694996361663029e-03
  };
};
}
}
