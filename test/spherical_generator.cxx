#include "catch2/catch_all.hpp"
#include <iostream>

#include <integratorxx/quadratures/radial.hpp>
#include <integratorxx/quadratures/s2.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>
#include <integratorxx/generators/spherical_factory.hpp>

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

TEST_CASE( "Pruning Schemes", "[sph-gen]" ) {

  using namespace IntegratorXX;

  SECTION("Ahrens-Beylkin") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, 99, 1.0,
      AngularQuad::AhrensBeylkin, 372
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 25, AngularQuad::AhrensBeylkin, 72},
          {25, 50, AngularQuad::AhrensBeylkin, 312},
          {50, 99, AngularQuad::AhrensBeylkin, 372},
        }
      };
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 34, AngularQuad::AhrensBeylkin, 72},
          {34, 50, AngularQuad::AhrensBeylkin, 72},
          {50, 99, AngularQuad::AhrensBeylkin, 372},
        }
      };
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Delley") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, 99, 1.0,
      AngularQuad::Delley, 302
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 25, AngularQuad::Delley, 26},
          {25, 50, AngularQuad::Delley, 194},
          {50, 99, AngularQuad::Delley, 302},
        }
      };
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 34, AngularQuad::Delley, 26},
          {34, 50, AngularQuad::Delley, 50},
          {50, 99, AngularQuad::Delley, 302},
        }
      };
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Lebedev-Laikov") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, 99, 1.0,
      AngularQuad::LebedevLaikov, 302
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 25, AngularQuad::LebedevLaikov, 26},
          {25, 50, AngularQuad::LebedevLaikov, 194},
          {50, 99, AngularQuad::LebedevLaikov, 302},
        }
      };
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 34, AngularQuad::LebedevLaikov, 26},
          {34, 50, AngularQuad::LebedevLaikov, 50},
          {50, 99, AngularQuad::LebedevLaikov, 302},
        }
      };
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

  SECTION("Womersley") {

    UnprunedSphericalGridSpecification unp{
      RadialQuad::MuraKnowles, 99, 1.0,
      AngularQuad::Womersley, 339
    };

    PrunedSphericalGridSpecification pruning_spec, ref_pruning_spec;

    SECTION("Robust") {
      pruning_spec = robust_psi4_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 25, AngularQuad::Womersley, 32},
          {25, 50, AngularQuad::Womersley, 201},
          {50, 99, AngularQuad::Womersley, 339},
        }
      };
    }

    SECTION("Treutler") {
      pruning_spec = treutler_pruning_scheme(unp);
      ref_pruning_spec = PrunedSphericalGridSpecification{
        RadialQuad::MuraKnowles, 99, 1.0,
        {
          {0 , 34, AngularQuad::Womersley, 32},
          {34, 50, AngularQuad::Womersley, 72},
          {50, 99, AngularQuad::Womersley, 339},
        }
      };
    }

    REQUIRE(pruning_spec == ref_pruning_spec);
  }

}

using bk_type  = IntegratorXX::Becke<double,double>;
using mk_type  = IntegratorXX::MuraKnowles<double,double>;
using mhl_type = IntegratorXX::MurrayHandyLaming<double,double>;
using ta_type  = IntegratorXX::TreutlerAhlrichs<double,double>;

using ah_type = IntegratorXX::AhrensBeylkin<double>;
using de_type = IntegratorXX::Delley<double>;
using ll_type = IntegratorXX::LebedevLaikov<double>;
using wo_type = IntegratorXX::Womersley<double>;

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
  UnprunedSphericalGridSpecification unp{
    rad_from_type<radial_type>(), nrad, 1.0,
    ang_from_type<angular_type>(), nang
  };

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
