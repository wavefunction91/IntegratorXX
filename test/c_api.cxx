#include "catch2/catch_all.hpp"
#include "quad_matcher.hpp"
#include <integratorxx/quadratures/all.hpp>
#include <integratorxx/c_api.h>


using namespace IntegratorXX;

TEST_CASE("C API") {

  intxx_quad_type quad;
  int error;

  SECTION("Invalid") {
    error = intxx_quad_init(&quad, -1);
    REQUIRE(error == INTXX_INVALID_QUAD);
    REQUIRE(quad.info == NULL);
    REQUIRE(quad.npoints == -1);
  }


  SECTION("Primitive") {

    const char* name;
    const int base_npts = 100;
    using base_quad_type = QuadratureBase<std::vector<double>, std::vector<double>>;
    std::unique_ptr<base_quad_type> base_quad = nullptr;

    SECTION("Uniform") {
      using quad_type = UniformTrapezoid<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_UNIFORM);
      name = "UNIFORM";
      base_quad = std::make_unique<quad_type>(base_npts, 0.0, 1.0);
    }

    SECTION("Gauss-Legendre") {
      using quad_type = GaussLegendre<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSLEG);
      name = "GAUSS_LEGENDRE";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Lobatto") {
      using quad_type = GaussLobatto<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSLOB);
      name = "GAUSS_LOBATTO";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 1") {
      using quad_type = GaussChebyshev1<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_1);
      name = "GAUSS_CHEBYSHEV_1";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 2") {
      using quad_type = GaussChebyshev2<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_2);
      name = "GAUSS_CHEBYSHEV_2";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 2MOD") {
      using quad_type = GaussChebyshev2Modified<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_2MOD);
      name = "GAUSS_CHEBYSHEV_2MOD";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 3") {
      using quad_type = GaussChebyshev3<double,double>;
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_3);
      name = "GAUSS_CHEBYSHEV_3";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad._state == NULL);
    REQUIRE(quad.info->ext_params.n == 0);
    REQUIRE(!strcmp(quad.info->name, name));

    // Get before set
    int npts;
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_INVALID_OUT);
    REQUIRE(npts == -1);

    // Set NPTS
    error = intxx_quad_set_npts(&quad, base_npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(quad.npoints == base_npts);

    // Get NPTS
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(npts == base_npts);

    // Check Quadrature Generation and Destruction
    if(base_quad) {
      intxx_generate_quad(&quad);
      REQUIRE(quad._state != NULL);

      // Check validity of the state
      auto state_as_quad = reinterpret_cast<base_quad_type*>(quad._state);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE_THAT(state_as_quad->points()[i], Matchers::WithinAbs(name, base_quad->points()[i], 1e-15));
        REQUIRE_THAT(state_as_quad->weights()[i], Matchers::WithinAbs(name, base_quad->weights()[i], 1e-15));
      }

      intxx_destroy_quad(&quad);
      REQUIRE(quad._state == NULL);
    }
  }


  SECTION("Radial") {
    const char* name;
    const int base_npts = 100;
    using base_quad_type = QuadratureBase<std::vector<double>, std::vector<double>>;
    std::unique_ptr<base_quad_type> base_quad_default = nullptr;

    SECTION("Becke") {
      //using quad_type = UniformTrapezoid<double,double>;
      error = intxx_quad_init(&quad, INTXX_RADQ_BECKE);
      name = "BECKE";
      //base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("MHL") {
      //using quad_type = UniformTrapezoid<double,double>;
      error = intxx_quad_init(&quad, INTXX_RADQ_MHL);
      name = "MURRAY_HANDY_LAMING";
      //base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("TA") {
      //using quad_type = UniformTrapezoid<double,double>;
      error = intxx_quad_init(&quad, INTXX_RADQ_TA);
      name = "TREUTLER_AHLRICHS";
      //base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("MK") {
      //using quad_type = UniformTrapezoid<double,double>;
      error = intxx_quad_init(&quad, INTXX_RADQ_MK);
      name = "MURA_KNOWLES";
      //base_quad = std::make_unique<quad_type>(base_npts);
    }

    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad._state == NULL);
    REQUIRE(!strcmp(quad.info->name, name));

    // Get before set
    int npts;
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_INVALID_OUT);
    REQUIRE(npts == -1);

    // Set NPTS
    error = intxx_quad_set_npts(&quad, base_npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(quad.npoints == base_npts);

    // Get NPTS
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(npts == base_npts);
  }

  intxx_quad_end(&quad);

}
