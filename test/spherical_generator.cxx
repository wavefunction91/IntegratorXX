#include "catch2/catch_all.hpp"
#include <iostream>

#include <integratorxx/quadratures/radial.hpp>
#include <integratorxx/quadratures/s2.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>
#include <integratorxx/generators/spherical_factory.hpp>
#include <integratorxx/generators/radial_factory.hpp>

using bk_type  = IntegratorXX::Becke<double,double>;
using mk_type  = IntegratorXX::MuraKnowles<double,double>;
using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
using ta_type  = IntegratorXX::TreutlerAhlrichs<double,double>;

using ah_type = IntegratorXX::AhrensBeylkin<double>;
using de_type = IntegratorXX::Delley<double>;
using ll_type = IntegratorXX::LebedevLaikov<double>;
using wo_type = IntegratorXX::Womersley<double>;

std::ostream& operator<<( std::ostream& out, IntegratorXX::PruningRegion p) {

  using namespace IntegratorXX;
  out << "ANGULAR QUAD: [" << p.idx_st << "," << p.idx_en << ") ";
  switch(p.angular_quad) {
    case AngularQuad::AhrensBeylkin:  out << "AH"; break;
    case AngularQuad::Delley:         out << "DE"; break;
    case AngularQuad::LebedevLaikov:  out << "LL"; break;
    case AngularQuad::Womersley:      out << "WO"; break;
  }
  out << "  SIZE = " << p.angular_size;
  
  return out;
}
std::ostream& operator<<( std::ostream& out, 
  IntegratorXX::PrunedSphericalGridSpecification p) {


  using namespace IntegratorXX;
  out << "RADIAL QUAD: ";
  switch(p.radial_quad) {
    case RadialQuad::Becke:             out << "Becke"; break;
    case RadialQuad::MuraKnowles:       out << "MK";    break;
    case RadialQuad::MurrayHandyLaming: out << "MHL";   break;
    case RadialQuad::TreutlerAhlrichs:  out << "TA";    break;
  } 

  out << std::endl;
  for(auto pr : p.pruning_regions) {
    out << "  " << pr << std::endl;
  }

  
  return out;
}

TEST_CASE( "String Getter", "[sph-gen]" ) {
  using namespace IntegratorXX;
  SECTION("Radial") {
    REQUIRE(radial_from_type<bk_type>()  == radial_from_string("Becke"));
    REQUIRE(radial_from_type<mk_type>()  == radial_from_string("MuraKnowles"));
    REQUIRE(radial_from_type<mk_type>()  == radial_from_string("MK"));
    REQUIRE(radial_from_type<mhl_type>() == radial_from_string("MurrayHandyLaming"));
    REQUIRE(radial_from_type<mhl_type>() == radial_from_string("MHL"));
    REQUIRE(radial_from_type<ta_type>() == radial_from_string("TreutlerAhlrichs"));
    REQUIRE(radial_from_type<ta_type>() == radial_from_string("TA"));
  }

  SECTION("Angular") {
    REQUIRE(angular_from_type<ah_type>()  == angular_from_string("AhrensBeylkin"));
    REQUIRE(angular_from_type<ah_type>()  == angular_from_string("AB"));
    REQUIRE(angular_from_type<de_type>()  == angular_from_string("Delley"));
    REQUIRE(angular_from_type<ll_type>()  == angular_from_string("LebedevLaikov"));
    REQUIRE(angular_from_type<ll_type>()  == angular_from_string("Lebedev"));
    REQUIRE(angular_from_type<ll_type>()  == angular_from_string("LL"));
    REQUIRE(angular_from_type<wo_type>()  == angular_from_string("Womersley"));
  }
}

TEST_CASE( "Pruning Schemes", "[sph-gen]" ) {

  using namespace IntegratorXX;

  auto rad_traits = make_radial_traits(RadialQuad::MuraKnowles, 99, 1.0);

  SECTION("Ahrens-Beylkin") {

    //UnprunedSphericalGridSpecification unp = make_unpruned_spec(
    //  AngularQuad::AhrensBeylkin, 372,
    //  RadialQuad::MuraKnowles, 99, 1.0
    //);
    UnprunedSphericalGridSpecification unp(
      RadialQuad::MuraKnowles, *rad_traits,
      AngularQuad::AhrensBeylkin, 372
    );

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 25, AngularQuad::AhrensBeylkin, 72},
          {25, 50, AngularQuad::AhrensBeylkin, 312},
          {50, 99, AngularQuad::AhrensBeylkin, 372},
        }
      );
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 34, AngularQuad::AhrensBeylkin, 72},
          {34, 50, AngularQuad::AhrensBeylkin, 72},
          {50, 99, AngularQuad::AhrensBeylkin, 372},
        }
      );
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Delley") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, *rad_traits,
      AngularQuad::Delley, 302
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles,  *rad_traits,
        std::vector<PruningRegion>{
          {0 , 25, AngularQuad::Delley, 26},
          {25, 50, AngularQuad::Delley, 194},
          {50, 99, AngularQuad::Delley, 302},
        }
      );
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 34, AngularQuad::Delley, 26},
          {34, 50, AngularQuad::Delley, 50},
          {50, 99, AngularQuad::Delley, 302},
        }
      );
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Lebedev-Laikov") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, *rad_traits,
      AngularQuad::LebedevLaikov, 302
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 25, AngularQuad::LebedevLaikov, 26},
          {25, 50, AngularQuad::LebedevLaikov, 194},
          {50, 99, AngularQuad::LebedevLaikov, 302},
        }
      );
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 34, AngularQuad::LebedevLaikov, 26},
          {34, 50, AngularQuad::LebedevLaikov, 50},
          {50, 99, AngularQuad::LebedevLaikov, 302},
        }
      );
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Womersley") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, *rad_traits,
      AngularQuad::Womersley, 339
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 25, AngularQuad::Womersley, 32},
          {25, 50, AngularQuad::Womersley, 201},
          {50, 99, AngularQuad::Womersley, 339},
        }
      );
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification(
        RadialQuad::MuraKnowles, *rad_traits,
        std::vector<PruningRegion>{
          {0 , 34, AngularQuad::Womersley, 32},
          {34, 50, AngularQuad::Womersley, 72},
          {50, 99, AngularQuad::Womersley, 339},
        }
      );
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

}

using radial_test_types = std::tuple<
  bk_type, mk_type, mhl_type, ta_type
>;

TEMPLATE_LIST_TEST_CASE("Radial Generator", "[sph-gen]", radial_test_types) {
  using namespace IntegratorXX;
  using radial_type = TestType;

  size_t npts = 10;
  radial_type rq(npts, 1.0);
  auto rad_spec = radial_from_type<radial_type>();
  auto rad_traits = make_radial_traits(rad_spec, npts, 1.0);
  auto rad_grid = RadialFactory::generate(rad_spec, *rad_traits);

  REQUIRE(rad_grid->npts() == npts);
  for(auto i = 0; i < npts; ++i) {
    auto pt = rad_grid->points()[i];
    auto pt_ref = rq.points()[i];
    REQUIRE_THAT(pt, Catch::Matchers::WithinAbs(pt_ref,1e-15));
  
    auto w = rad_grid->weights()[i];
    auto w_ref = rq.weights()[i];
    REQUIRE_THAT(w, Catch::Matchers::WithinAbs(w_ref,1e-15));
  }
}

using s2_test_types = std::tuple<
  ah_type, de_type, ll_type, wo_type
>;

TEMPLATE_LIST_TEST_CASE("S2 Generator", "[sph-gen]", s2_test_types) {
  using namespace IntegratorXX;
  using angular_type = TestType;
  using angular_traits = quadrature_traits<angular_type>;

  size_t npts = angular_traits::npts_by_algebraic_order(
    angular_traits::next_algebraic_order(1)); // Smallest possible angular grid

  angular_type aq(npts);
  
  auto ang_spec = angular_from_type<angular_type>();
  auto ang_grid = S2Factory::generate(ang_spec, npts);
  REQUIRE(ang_grid->npts() == npts);
  for(auto i = 0; i < npts; ++i) {
    auto pt = ang_grid->points()[i];
    auto pt_ref = aq.points()[i];
    REQUIRE_THAT(pt[0], Catch::Matchers::WithinAbs(pt_ref[0],1e-15));
    REQUIRE_THAT(pt[1], Catch::Matchers::WithinAbs(pt_ref[1],1e-15));
    REQUIRE_THAT(pt[2], Catch::Matchers::WithinAbs(pt_ref[2],1e-15));
  
    auto w = ang_grid->weights()[i];
    auto w_ref = aq.weights()[i];
    REQUIRE_THAT(w, Catch::Matchers::WithinAbs(w_ref,1e-15));
  }

}

using sph_test_types = std::tuple<
  std::tuple<bk_type, ah_type>,
  std::tuple<bk_type, de_type>,
  std::tuple<bk_type, ll_type>,
  std::tuple<bk_type, wo_type>,

  std::tuple<mk_type, ah_type>,
  std::tuple<mk_type, de_type>,
  std::tuple<mk_type, ll_type>,
  std::tuple<mk_type, wo_type>,

  std::tuple<mhl_type, ah_type>,
  std::tuple<mhl_type, de_type>,
  std::tuple<mhl_type, ll_type>,
  std::tuple<mhl_type, wo_type>,

  std::tuple<ta_type, ah_type>,
  std::tuple<ta_type, de_type>,
  std::tuple<ta_type, ll_type>,
  std::tuple<ta_type, wo_type>
>;



TEMPLATE_LIST_TEST_CASE("Unpruned", "[sph-gen]", sph_test_types) {
  using namespace IntegratorXX;
  using radial_type  = std::decay_t<decltype(std::get<0>(std::declval<TestType>()))>;
  using angular_type = std::decay_t<decltype(std::get<1>(std::declval<TestType>()))>;
  using angular_traits = quadrature_traits<angular_type>;

  using spherical_type = SphericalQuadrature<radial_type,angular_type>;

  size_t nrad = 10;
  size_t nang = angular_traits::npts_by_algebraic_order(
    angular_traits::next_algebraic_order(1)); // Smallest possible angular grid

  // Generate the quadrature manually
  radial_type rq(nrad, 1.0);
  angular_type aq(nang);
  spherical_type sph_ref(rq,aq);


  // Generate via runtime API
  auto rad_spec = radial_from_type<radial_type>();
  auto rad_traits = make_radial_traits(rad_spec, nrad, 1.0);
  UnprunedSphericalGridSpecification unp(
    rad_spec, *rad_traits,
    angular_from_type<angular_type>(), nang
  );

  auto sph = SphericalGridFactory::generate_grid(unp);

  // Check that they're the same
  REQUIRE(sph->npts() == sph_ref.npts());

  const auto npts = sph->npts();
  for(auto i = 0; i < npts; ++i) {
    auto pt = sph->points()[i];
    auto pt_ref = sph_ref.points()[i];
    REQUIRE_THAT(pt[0], Catch::Matchers::WithinAbs(pt_ref[0],1e-15));
    REQUIRE_THAT(pt[1], Catch::Matchers::WithinAbs(pt_ref[1],1e-15));
    REQUIRE_THAT(pt[2], Catch::Matchers::WithinAbs(pt_ref[2],1e-15));

    auto w = sph->weights()[i];
    auto w_ref = sph_ref.weights()[i];
    REQUIRE_THAT(w, Catch::Matchers::WithinAbs(w_ref,1e-15));
  }
  
}

TEMPLATE_LIST_TEST_CASE("Pruned", "[sph-gen]", sph_test_types) {
  using namespace IntegratorXX;
  using radial_type  = std::decay_t<decltype(std::get<0>(std::declval<TestType>()))>;
  using angular_type = std::decay_t<decltype(std::get<1>(std::declval<TestType>()))>;
  using angular_traits = quadrature_traits<angular_type>;

  using spherical_type = PrunedSphericalQuadrature<radial_type,angular_type>;

  size_t nrad = 99;
  size_t nang = angular_traits::npts_by_algebraic_order(
    angular_traits::next_algebraic_order(29)); // Smallest possible angular grid

  // Generate pruning scheme
  auto rad_spec = radial_from_type<radial_type>();
  auto rad_traits = make_radial_traits(rad_spec, nrad, 1.0);
  UnprunedSphericalGridSpecification unp(
    rad_spec, *rad_traits,
    angular_from_type<angular_type>(), nang
  );

  auto pruning_spec = robust_psi4_pruning_scheme(unp);


  // Generate the quadrature manually
  radial_type rq(nrad, 1.0);
  RadialGridPartition<angular_type> rgp;
  for(auto pr : pruning_spec.pruning_regions) {
    angular_type aq(pr.angular_size);
    rgp.add_quad(rq, pr.idx_st, aq);
  }
  rgp.finalize(rq);

  spherical_type sph_ref(rq, rgp);


  // Generate via runtime API

  auto sph = SphericalGridFactory::generate_grid(pruning_spec);

  // Check that they're the same
  REQUIRE(sph->npts() == sph_ref.npts());

  const auto npts = sph->npts();
  for(auto i = 0; i < npts; ++i) {
    auto pt = sph->points()[i];
    auto pt_ref = sph_ref.points()[i];
    REQUIRE_THAT(pt[0], Catch::Matchers::WithinAbs(pt_ref[0],1e-15));
    REQUIRE_THAT(pt[1], Catch::Matchers::WithinAbs(pt_ref[1],1e-15));
    REQUIRE_THAT(pt[2], Catch::Matchers::WithinAbs(pt_ref[2],1e-15));

    auto w = sph->weights()[i];
    auto w_ref = sph_ref.weights()[i];
    REQUIRE_THAT(w, Catch::Matchers::WithinAbs(w_ref,1e-15));
  }
  
}
