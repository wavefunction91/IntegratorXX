#include "catch2/catch_all.hpp"
#include <integratorxx/quadratures/all.hpp>
#include <integratorxx/util/lambert_w.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>

#include "quad_matcher.hpp"
#include "test_functions.hpp"

using namespace IntegratorXX;

inline constexpr double inf = std::numeric_limits<double>::infinity();
inline constexpr double eps = std::numeric_limits<double>::epsilon();

template <typename T = double>
constexpr T gaussian( T alpha, T c, T x ) {
  return std::exp( -alpha * (x-c) * (x-c) );
}

template <typename T = double>
constexpr T gaussian( T x ) {
  return gaussian( 1., 0., x );
}

template <typename T = double>
constexpr T ref_gaussian_int( T alpha, T c, T a, T b ) {

  const auto low = (a == inf) ? T(-1.) : std::erf( std::sqrt(alpha) * (c - a) );
  const auto hgh = (b == inf) ? T(-1.) : std::erf( std::sqrt(alpha) * (c - b) );

  return 0.5 * std::sqrt( M_PI / alpha ) * (low - hgh);
}

template <typename T = double>
constexpr T ref_gaussian_int( T a, T b ) {

  return ref_gaussian_int( 1., 0., a, b );

}


template <typename T>
auto chebyshev_T(int n, T x) {
  return std::cos( n * std::acos(x) );
}






TEST_CASE( "Gauss-Legendre Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    std::vector<double> cp(p+1, 0.0); 
    for(int i = 0; i < p; ++i) {
      cp[i] = c[i] / (p-i);
    }
    return Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussLegendre<double,double>;
  test_random_polynomial<quad_type, Polynomial>("Gauss-Legendre", 10, 14, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Lobatto Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    std::vector<double> cp(p+1, 0.0); 
    for(int i = 0; i < p; ++i) {
      cp[i] = c[i] / (p-i);
    }
    return Polynomial::evaluate(cp, 1.0) - Polynomial::evaluate(cp, -1.0);
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussLobatto<double,double>;
  test_random_polynomial<quad_type, Polynomial>("Gauss-Lobatto", 10, 14, 
    [](int o){ return 2*o-1; }, // Max order 2N-3
    ref_value, 1e-12 );


}

TEST_CASE( "Gauss-Chebyshev T1 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T1 evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
    for(int i = 0; i < p; ++i) {
      int k = p - i - 1;
      if(k > 0)
        ref += c[i] * (std::sqrt(M_PI) / k) * (std::pow(-1,k)+1) * 
               std::tgamma((k+1)/2.0) / std::tgamma(k/2.0);
      else
        ref += c[i] * M_PI;
    }
    return ref;
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussChebyshev1<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT1WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T1)", 10, 100, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Chebyshev T2 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T2 evaluated over [-1,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
    for(int i = 0; i < p; ++i) {
      int k = p - i - 1;
      if(k > 0)
        ref += c[i] * (std::sqrt(M_PI) / 4.0) * (std::pow(-1,k)+1) * 
               std::tgamma((k+1)/2.0) / std::tgamma(k/2.0 + 2);
      else
        ref += c[i] * M_PI/2.0;
    }
    return ref;
  };

  // Test Quadrature for Correctness
  using quad_type = IntegratorXX::GaussChebyshev2<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT2WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T2)", 10, 100, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}

TEST_CASE( "Gauss-Chebyshev T3 Quadratures", "[1d-quad]" ) {

  // Reference integral for polynomial * T3 evaluated over [0,1]
  auto ref_value = [](const std::vector<double>& c) {
    const auto p = c.size();
    double ref = 0.0;
      for(int i = 0; i < p; ++i) {
        int k = p - i - 1;
        ref += c[i] * std::sqrt(M_PI) * std::tgamma(k+1.5) / std::tgamma(k+2); 
      }
    return ref;
  };

  // Test Quadrature for Correctness
  // TODO: Code breaks down for large orders here
  using quad_type = IntegratorXX::GaussChebyshev3<double,double>;
  using func_type = WeightedPolynomial<ChebyshevT3WeightFunction>;
  test_random_polynomial<quad_type, func_type>("Gauss-Chebyshev (T2)", 10, 50, 
    [](int o){ return 2*o+1; }, // Max order 2N-1
    ref_value, 1e-12 );

}


TEST_CASE( "Euler-Maclaurin Quadratures by Murray, Handy, and Laming", "[1d-quad]" ) {
  IntegratorXX::MurrayHandyLaming<double,double> quad(150);
  const auto msg = "Euler-Maclaurin N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Treutler-Ahlrichs Quadratures", "[1d-quad]" ) {
  IntegratorXX::TreutlerAhlrichs<double,double> quad(150);
  const auto msg = "Treutler-Ahlrichs N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Mura-Knowles Quadratures", "[1d-quad]" ) {
  IntegratorXX::MuraKnowles<double,double> quad(350);
  const auto msg = "Mura-Knowles N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Becke Quadratures", "[1d-quad]" ) {
  IntegratorXX::Becke<double,double> quad(350);
  const auto msg = "Becke N = " + std::to_string(quad.npts());
  test_quadrature<RadialGaussian>(msg, quad, std::sqrt(M_PI)/4, 1e-10);
}

TEST_CASE( "Lebedev-Laikov", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {
    IntegratorXX::LebedevLaikov<double> quad( nPts );
    const auto msg = "Lebedev-Laikov N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);
  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}

TEST_CASE( "Ahrens-Beylkin", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::AhrensBeylkin<double> quad( nPts );
    const auto msg = "Ahrens-Beylkin N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(312);
  test_fn(792);
  test_fn(972);
}

TEST_CASE( "Womersley", "[1d-quad]" ) {

  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::Womersley<double> quad( nPts );
    const auto msg = "Womersley N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(314);
  test_fn(801);
  test_fn(969);
}

TEST_CASE( "Delley", "[1d-quad]" ) {


  auto test_fn = [&]( size_t nPts ) {

    IntegratorXX::Delley<double> quad( nPts );
    const auto msg = "Delley N = " + std::to_string(quad.npts());
    test_angular_quadrature(msg, quad, 10, 1e-10);

  };

  test_fn(302);
  test_fn(770);
  test_fn(974);

}

TEST_CASE("LMG", "[1d-quad]") {
  using namespace IntegratorXX;
  SECTION("Lambert W") {
    size_t n = 50;

    SECTION("-1/e") {
      REQUIRE(lambert_w0(-1.0 / std::exp(1.0)) == -1.0);
    }

    double dx;
    double x_min;
    std::vector<double> w_ref;
    int branch = 0;
    SECTION("W0 (-1/e,0)") {
      dx    = (1.0/std::exp(1.0)) / (n-1);
      x_min = -1.0/std::exp(1.0); 
      w_ref = { -1.0000000000000000e+00, -8.1043436396882251e-01, -7.3839124833180692e-01,
                -6.8534285732547906e-01, -6.4201643964687316e-01, -6.0484100186256395e-01,
                -5.7199479483565563e-01, -5.4240102237003140e-01, -5.1536144833012332e-01,
                -4.9039331229166899e-01, -4.6714668763808043e-01, -4.4535855202428648e-01,
                -4.2482542402809581e-01, -4.0538614978331228e-01, -3.8691058975718701e-01,
                -3.6929191010441786e-01, -3.5244116871063264e-01, -3.3628341359547681e-01,
                -3.2075480804994166e-01, -3.0580047091961676e-01, -2.9137282629002875e-01,
                -2.7743032325972444e-01, -2.6393642934736894e-01, -2.5085882941845555e-01,
                -2.3816878116899459e-01, -2.2584059140279072e-01, -2.1385118659837130e-01,
                -2.0217975786260850e-01, -1.9080746514340260e-01, -1.7971718907514189e-01,
                -1.6889332142980482e-01, -1.5832158709811098e-01, -1.4798889200596840e-01,
                -1.3788319250617834e-01, -1.2799338266282390e-01, -1.1830919653009347e-01,
                -1.0882112306525034e-01, -9.9520331741533477e-02, -9.0398607266594919e-02,
                -8.1448292084908017e-02, -7.2662235562987412e-02, -6.4033748935368034e-02,
                -5.5556565235647053e-02, -4.7224803557117209e-02, -3.9032937086775930e-02,
                -3.0975764438825409e-02, -2.3048383882382367e-02, -1.5246170115525802e-02,
                -7.5647532860519534e-03, 6.4184768611141850e-17 };
    }
    SECTION("W0 (0,e)") {
      dx    = std::exp(1.0) / (n-1);
      x_min = 0.0; 
      w_ref = { 0.0000000000000000e+00, 5.2630934056544072e-02, 1.0035620951623317e-01,
                1.4409241373767120e-01, 1.8451244238478223e-01, 2.2212638256855327e-01,
                2.5733105242307630e-01, 2.9044182654506334e-01, 3.2171389544393414e-01,
                3.5135693697900772e-01, 3.7954552465616309e-01, 4.0642668865609727e-01,
                4.3212552273462468e-01, 4.5674941777919742e-01, 4.8039130984974648e-01,
                5.0313220781385692e-01, 5.2504318560584629e-01, 5.4618697067229260e-01,
                5.6661922372849582e-01, 5.8638957965403482e-01, 6.0554250149726696e-01,
                6.2411798675730124e-01, 6.4215215580659113e-01, 6.5967774546581859e-01,
                6.7672452563807517e-01, 6.9331965306433241e-01, 7.0948797333649716e-01,
                7.2525228005705800e-01, 7.4063353829291223e-01, 7.5565107811039622e-01,
                7.7032276290732105e-01, 7.8466513640828828e-01, 7.9869355151124433e-01,
                8.1242228362814539e-01, 8.2586463072191718e-01, 8.3903300188356378e-01,
                8.5193899600034140e-01, 8.6459347182518476e-01, 8.7700661055879869e-01,
                8.8918797189091481e-01, 9.0114654430979890e-01, 9.1289079037410314e-01,
                9.2442868754455021e-01, 9.3576776509144632e-01, 9.4691513752504819e-01,
                9.5787753493721217e-01, 9.6866133059280313e-01, 9.7927256606663837e-01,
                9.8971697418509541e-01, 9.9999999999999989e-01 };
    }

    SECTION("W0 (e,inf)") {
      dx    = 1;
      x_min = std::exp(1.0); 
      w_ref = { 1.0000000000000000e+00, 1.1626015113014869e+00, 1.2938344571456706e+00,
                1.4042003708950586e+00, 1.4996204192818812e+00, 1.5837783960587868e+00,
                1.6591292586935935e+00, 1.7273945838205316e+00, 1.7898301426340337e+00,
                1.8473811371144686e+00, 1.9007774734056948e+00, 1.9505949530993991e+00,
                1.9972960780278193e+00, 2.0412581329768691e+00, 2.0827930401857122e+00,
                2.1221617266541868e+00, 2.1595847338181962e+00, 2.1952501935697217e+00,
                2.2293199201838059e+00, 2.2619341295913213e+00, 2.2932151421500375e+00,
                2.3232703215087871e+00, 2.3521944316956782e+00, 2.3800715457358486e+00,
                2.4069766047089463e+00, 2.4329767015585579e+00, 2.4581321461297474e+00,
                2.4824973548124829e+00, 2.5061215984360530e+00, 2.5290496347488083e+00,
                2.5513222462703027e+00, 2.5729767000540669e+00, 2.5940471426162146e+00,
                2.6145649407274134e+00, 2.6345589767577975e+00, 2.6540559056765503e+00,
                2.6730803795436775e+00, 2.6916552443184649e+00, 2.7098017129924283e+00,
                2.7275395183923634e+00, 2.7448870484592360e+00, 2.7618614663662542e+00,
                2.7784788174751269e+00, 2.7947541248280854e+00, 2.8107014746227152e+00,
                2.8263340929076008e+00, 2.8416644145615426e+00, 2.8567041454717463e+00,
                2.8714643187019151e+00, 2.8859553453357312e+00 };
    }
    SECTION("W1 (-1/e,-0.25)") {
      branch = -1;
      dx    = (1.0/std::exp(1.0) - 0.25) / (n-1);
      x_min = -1.0/std::exp(1.0); 
      w_ref = { -1.0000000000000000e+00, -1.1189650446807282e+00, -1.1711582725750418e+00,
                -1.2124855589488057e+00, -1.2482413372948526e+00, -1.2804726429337259e+00,
                -1.3102276924123728e+00, -1.3381283510412803e+00, -1.3645793574591196e+00,
                -1.3898616085396616e+00, -1.4141795078524806e+00, -1.4376873219364974e+00,
                -1.4605049067902260e+00, -1.4827276158806526e+00, -1.5044328217548120e+00,
                -1.5256843663192678e+00, -1.5465356909502610e+00, -1.5670320955171273e+00,
                -1.5872124053596470e+00, -1.6071102254331719e+00, -1.6267549000718433e+00,
                -1.6461722586522056e+00, -1.6653852027928755e+00, -1.6844141744163048e+00,
                -1.7032775329676604e+00, -1.7219918624772559e+00, -1.7405722238112493e+00,
                -1.7590323636442364e+00, -1.7773848889282362e+00, -1.7956414136079653e+00,
                -1.8138126828282584e+00, -1.8319086787493766e+00, -1.8499387112277788e+00,
                -1.8679114959619136e+00, -1.8858352221933619e+00, -1.9037176116563133e+00,
                -1.9215659701557612e+00, -1.9393872329070991e+00, -1.9571880045721319e+00,
                -1.9749745947677455e+00, -1.9927530496951602e+00, -2.0105291804333940e+00,
                -2.0283085883553316e+00, -2.0460966880547553e+00, -2.0638987281149239e+00,
                -2.0817198100013492e+00, -2.0995649053215297e+00, -2.1174388716610637e+00,
                -2.1353464671775746e+00, -2.1532923641103534e+00 };
    }
    SECTION("W1 (-0.25,0)") {
      branch = -1;
      dx    = (0.25-1e-10) / (n-1);
      x_min = -0.25; 
      w_ref = { -2.1532923641103499e+00, -2.1914997961914482e+00, -2.2299431488704058e+00,
                -2.2686646700353417e+00, -2.3077062240124904e+00, -2.3471096678770036e+00,
                -2.3869172040331206e+00, -2.4271717185362274e+00, -2.4679171134091620e+00,
                -2.5091986404347137e+00, -2.5510632435016558e+00, -2.5935599164917043e+00,
                -2.6367400838934105e+00, -2.6806580118111714e+00, -2.7253712578170455e+00,
                -2.7709411692030885e+00, -2.8174334406865955e+00, -2.8649187445794793e+00,
                -2.9134734489670580e+00, -2.9631804427033619e+00, -3.0141300902276456e+00,
                -3.0664213446240929e+00, -3.1201630543731218e+00, -3.1754755084163060e+00,
                -3.2324922762280690e+00, -3.2913624156132215e+00, -3.3522531424379425e+00,
                -3.4153530856305858e+00, -3.4808762907481663e+00, -3.5490671909401166e+00,
                -3.6202068424308829e+00, -3.6946208337432305e+00, -3.7726894411668441e+00,
                -3.8548608453419511e+00, -3.9416685911800489e+00, -4.0337550432658453e+00,
                -4.1319034964882038e+00, -4.2370830903073911e+00, -4.3505132002644684e+00,
                -4.4737584327019935e+00, -4.6088735602501680e+00, -4.7586337059829287e+00,
                -4.9269181859532187e+00, -5.1193905946137361e+00, -5.3448010688929983e+00,
                -5.6177518509657292e+00, -5.9654942798986532e+00, -6.4488728297255040e+00,
                -7.2605734538864919e+00, -2.6295240572605621e+01 };
    }

    double x = x_min;
    for(auto i = 0; i < w_ref.size(); ++i) {
      auto w = branch == 0 ? lambert_w0(x) : lambert_wm1(x);
      REQUIRE_THAT(w, Catch::Matchers::WithinAbs(w_ref[i], 1e-15));
      //printf("%.16e %.16e %.4e\n", x, lambert_wm1(x),
      //  lambert_wm1(x) - boost::math::lambert_wm1(x));
      x = x + dx;
    }

  }

  SECTION("LMG Radial Bounds") {
    // aug-cc-pVTZ Hydrogen Data
    SECTION("H aug-cc-pVTZ") {
      const double alpha_min_0 = 0.02526;
      const double alpha_min_1 = 0.10200;
      const double alpha_min_2 = 0.24700;

      const double alpha_max_0 = 33.87;
      //const double alpha_max_1 = 1.407;
      const double alpha_max_2 = 1.057;

      auto r_lower_0 = lmg::r_lower(0, alpha_max_0, 1e-12);
      auto r_lower_2 = lmg::r_lower(2, alpha_max_2, 1e-12);

      auto r_upper_0 = lmg::r_upper(0, alpha_min_0, 1e-12);
      auto r_upper_1 = lmg::r_upper(1, alpha_min_1, 1e-12);
      auto r_upper_2 = lmg::r_upper(2, alpha_min_2, 1e-12);


      REQUIRE_THAT(r_lower_0, Catch::Matchers::WithinAbs(3.2370214224347012e-05, 1e-15));
      REQUIRE_THAT(r_lower_2, Catch::Matchers::WithinAbs(3.1703237193762383e-03, 1e-15));

      REQUIRE_THAT(r_upper_0, Catch::Matchers::WithinAbs(3.3998088412766542e+01, 1e-15));
      REQUIRE_THAT(r_upper_1, Catch::Matchers::WithinAbs(1.7452224096146878e+01, 1e-15));
      REQUIRE_THAT(r_upper_2, Catch::Matchers::WithinAbs(1.1588086549342464e+01, 1e-15));
    }

    // aug-cc-pVTZ Carbon Data
    SECTION("C aug-cc-pVTZ") {
      const double alpha_min_0 = 3.3423e-02;
      const double alpha_min_1 = 3.3258e-02;
      const double alpha_min_2 = 7.6815e-02;
      const double alpha_min_3 = 1.7636e-01;;
      const double alpha_min_4 = 5.1329e-01;

      const double alpha_max_0 = 6.177194e+06;
      const double alpha_max_2 = 8.857680e+01;
      const double alpha_max_4 = 1.445900e+00;

      auto r_lower_0 = lmg::r_lower(0, alpha_max_0, 1e-12);
      auto r_lower_2 = lmg::r_lower(2, alpha_max_2, 1e-12);
      auto r_lower_4 = lmg::r_lower(4, alpha_max_4, 1e-12);

      auto r_upper_0 = lmg::r_upper(0, alpha_min_0, 1e-12);
      auto r_upper_1 = lmg::r_upper(1, alpha_min_1, 1e-12);
      auto r_upper_2 = lmg::r_upper(2, alpha_min_2, 1e-12);
      auto r_upper_3 = lmg::r_upper(3, alpha_min_3, 1e-12);
      auto r_upper_4 = lmg::r_upper(4, alpha_min_4, 1e-12);

      REQUIRE_THAT(r_lower_0, Catch::Matchers::WithinAbs(7.5797965945427875e-08,1e-15));
      REQUIRE_THAT(r_lower_2, Catch::Matchers::WithinAbs(3.4632282095817583e-04,1e-15));
      REQUIRE_THAT(r_lower_4, Catch::Matchers::WithinAbs(1.1559748847781975e-02,1e-15));

      REQUIRE_THAT(r_upper_0, Catch::Matchers::WithinAbs(2.9556190534662171e+01,1e-15));
      REQUIRE_THAT(r_upper_1, Catch::Matchers::WithinAbs(3.0563480016561439e+01,1e-15));
      REQUIRE_THAT(r_upper_2, Catch::Matchers::WithinAbs(2.0779600292386970e+01,1e-15));
      REQUIRE_THAT(r_upper_3, Catch::Matchers::WithinAbs(1.4179981794582680e+01,1e-15));
      REQUIRE_THAT(r_upper_4, Catch::Matchers::WithinAbs(8.5952200832939525e+00,1e-15));
    }
  }


  SECTION("LMG Step Size") {
    REQUIRE_THAT(lmg::step_size(0, 1e-12), Catch::Matchers::WithinAbs(1.5235497564176387e-01,1e-15));
    REQUIRE_THAT(lmg::step_size(1, 1e-12), Catch::Matchers::WithinAbs(1.4579053103691331e-01,1e-15));
    REQUIRE_THAT(lmg::step_size(2, 1e-12), Catch::Matchers::WithinAbs(1.4028895246473957e-01,1e-15));
    REQUIRE_THAT(lmg::step_size(3, 1e-12), Catch::Matchers::WithinAbs(1.3554180272204153e-01,1e-15));
    REQUIRE_THAT(lmg::step_size(4, 1e-12), Catch::Matchers::WithinAbs(1.3136472924710046e-01,1e-15));
  }

}
