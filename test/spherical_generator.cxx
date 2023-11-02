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
