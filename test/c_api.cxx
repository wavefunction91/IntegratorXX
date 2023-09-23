#include "catch2/catch_all.hpp"
#include <integratorxx/quadratures/all.hpp>
#include <integratorxx/c_api.h>

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

    SECTION("Uniform") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_UNIFORM);
      name = "UNIFORM";
    }

    SECTION("Gauss-Legendre") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSLEG);
      name = "GAUSS_LEGENDRE";
    }

    SECTION("Gauss-Lobatto") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSLOB);
      name = "GAUSS_LOBATTO";
    }

    SECTION("Gauss-Chebyshev 1") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_1);
      name = "GAUSS_CHEBYSHEV_1";
    }

    SECTION("Gauss-Chebyshev 2") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_2);
      name = "GAUSS_CHEBYSHEV_2";
    }

    SECTION("Gauss-Chebyshev 2MOD") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_2MOD);
      name = "GAUSS_CHEBYSHEV_2MOD";
    }

    SECTION("Gauss-Chebyshev 3") {
      error = intxx_quad_init(&quad, INTXX_PRMQ_GAUSSCHEB_3);
      name = "GAUSS_CHEBYSHEV_3";
    }

    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad.info->ext_params.n == 0);
    REQUIRE(!strcmp(quad.info->name, name));

    // Set NPTS
    error = intxx_quad_set_npts(&quad, base_npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(quad.npoints == base_npts);

    // Get NPTS
    int npts;
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(npts == base_npts);
  }

  intxx_quad_end(&quad);

}
