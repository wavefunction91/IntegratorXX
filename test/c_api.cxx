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

    int quad_num;
    SECTION("Uniform") {
      using quad_type = UniformTrapezoid<double,double>;
      quad_num = INTXX_PRMQ_UNIFORM;
      name = "UNIFORM";
      base_quad = std::make_unique<quad_type>(base_npts, 0.0, 1.0);
    }

    SECTION("Gauss-Legendre") {
      using quad_type = GaussLegendre<double,double>;
      quad_num = INTXX_PRMQ_GAUSSLEG;
      name = "GAUSS_LEGENDRE";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Lobatto") {
      using quad_type = GaussLobatto<double,double>;
      quad_num = INTXX_PRMQ_GAUSSLOB;
      name = "GAUSS_LOBATTO";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 1") {
      using quad_type = GaussChebyshev1<double,double>;
      quad_num = INTXX_PRMQ_GAUSSCHEB_1;
      name = "GAUSS_CHEBYSHEV_1";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 2") {
      using quad_type = GaussChebyshev2<double,double>;
      quad_num = INTXX_PRMQ_GAUSSCHEB_2;
      name = "GAUSS_CHEBYSHEV_2";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 2MOD") {
      using quad_type = GaussChebyshev2Modified<double,double>;
      quad_num = INTXX_PRMQ_GAUSSCHEB_2MOD;
      name = "GAUSS_CHEBYSHEV_2MOD";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    SECTION("Gauss-Chebyshev 3") {
      using quad_type = GaussChebyshev3<double,double>;
      quad_num = INTXX_PRMQ_GAUSSCHEB_3;
      name = "GAUSS_CHEBYSHEV_3";
      base_quad = std::make_unique<quad_type>(base_npts);
    }

    // Initialize
    error = intxx_quad_init(&quad, quad_num);
    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad._state_quad == NULL);
    REQUIRE(quad.info->number == quad_num);
    REQUIRE(quad.info->kind == INTXX_PRM_QUAD);
    REQUIRE(quad.info->dim == 1);
    REQUIRE(quad.info->generate != NULL);
    REQUIRE(quad.info->destroy != NULL);
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
      REQUIRE(quad._state_quad != NULL);

      // Check validity of the state
      auto state_as_quad = reinterpret_cast<base_quad_type*>(quad._state_quad);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE_THAT(state_as_quad->points()[i], Matchers::WithinAbs(name, base_quad->points()[i], 1e-15));
        REQUIRE_THAT(state_as_quad->weights()[i], Matchers::WithinAbs(name, base_quad->weights()[i], 1e-15));
      }

      intxx_destroy_quad(&quad);
      REQUIRE(quad._state_quad == NULL);
    }
  }


  SECTION("Radial") {
    const char* name;
    const int base_npts = 100;
    const double RSCAL = 2.0;
    using base_quad_type = QuadratureBase<std::vector<double>, std::vector<double>>;
    std::unique_ptr<base_quad_type> base_quad_default = nullptr;
    std::unique_ptr<base_quad_type> base_quad_scaled  = nullptr;

    int quad_num;
    SECTION("Becke") {
      using quad_type = Becke<double,double>;
      using traits_type = typename quad_type::traits_type;
      quad_num = INTXX_RADQ_BECKE;
      name = "BECKE";
      base_quad_default = std::make_unique<quad_type>(base_npts);
      base_quad_scaled  = std::make_unique<quad_type>(base_npts,traits_type(RSCAL));
    }

    SECTION("MHL") {
      using quad_type = MurrayHandyLaming<double,double>;
      using traits_type = typename quad_type::traits_type;
      quad_num = INTXX_RADQ_MHL;
      name = "MURRAY_HANDY_LAMING";
      base_quad_default = std::make_unique<quad_type>(base_npts);
      base_quad_scaled  = std::make_unique<quad_type>(base_npts,traits_type(RSCAL));
    }

    SECTION("TA") {
      using quad_type = TreutlerAhlrichs<double,double>;
      using traits_type = typename quad_type::traits_type;
      quad_num = INTXX_RADQ_TA;
      name = "TREUTLER_AHLRICHS";
      base_quad_default = std::make_unique<quad_type>(base_npts);
      base_quad_scaled  = std::make_unique<quad_type>(base_npts,traits_type(RSCAL));
    }

    SECTION("MK") {
      using quad_type = MuraKnowles<double,double>;
      using traits_type = typename quad_type::traits_type;
      quad_num = INTXX_RADQ_MK;
      name = "MURA_KNOWLES";
      base_quad_default = std::make_unique<quad_type>(base_npts);
      base_quad_scaled  = std::make_unique<quad_type>(base_npts,traits_type(RSCAL));
    }

    // Initialize
    error = intxx_quad_init(&quad, quad_num);
    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad._state_quad == NULL);
    REQUIRE(quad.info->number == quad_num);
    REQUIRE(quad.info->kind == INTXX_RAD_QUAD);
    REQUIRE(quad.info->dim == 1);
    REQUIRE(quad.info->generate != NULL);
    REQUIRE(quad.info->destroy != NULL);
    REQUIRE(!strcmp(quad.info->name, name));


    // Check default parameters
    REQUIRE(quad.info->ext_params.n == 1);
    REQUIRE(!strcmp(quad.info->ext_params.names[0], "RAD_SCAL"));
    REQUIRE(!strcmp(quad.info->ext_params.descriptions[0], "Radial Scaling Factor"));

    double R;
    intxx_get_ext_param_name(&quad, "RAD_SCAL", &R);
    REQUIRE(R == 1.0);

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
    if(base_quad_default) {
      intxx_generate_quad(&quad);
      REQUIRE(quad._state_quad != NULL);

      // Check validity of the state
      auto state_as_quad = reinterpret_cast<base_quad_type*>(quad._state_quad);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE_THAT(state_as_quad->points()[i], Matchers::WithinAbs(name, base_quad_default->points()[i], 1e-15));
        REQUIRE_THAT(state_as_quad->weights()[i], Matchers::WithinAbs(name, base_quad_default->weights()[i], 1e-15));
      }

      intxx_destroy_quad(&quad);
      REQUIRE(quad._state_quad == NULL);
    }

    // Check Scaling Works
    if(base_quad_scaled) {
      intxx_set_ext_param_name(&quad, "RAD_SCAL", RSCAL);
    
      // Check change stuck
      intxx_get_ext_param_name(&quad, "RAD_SCAL", &R);
      REQUIRE(R == RSCAL);

      // Regenerate and check
      intxx_generate_quad(&quad);

      auto state_as_quad = reinterpret_cast<base_quad_type*>(quad._state_quad);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE_THAT(state_as_quad->points()[i], Matchers::WithinAbs(name, base_quad_scaled->points()[i], 1e-15));
        REQUIRE_THAT(state_as_quad->weights()[i], Matchers::WithinAbs(name, base_quad_scaled->weights()[i], 1e-15));
      }


      // Set without destroying to make sure quadratures are regerated
      intxx_set_ext_param_name(&quad, "RAD_SCAL", 1.0);
      state_as_quad = reinterpret_cast<base_quad_type*>(quad._state_quad);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE_THAT(state_as_quad->points()[i], Matchers::WithinAbs(name, base_quad_default->points()[i], 1e-15));
        REQUIRE_THAT(state_as_quad->weights()[i], Matchers::WithinAbs(name, base_quad_default->weights()[i], 1e-15));
      }
    
     
    }
  }

  SECTION("Angular") {
    const char* name;
    using base_quad_type = QuadratureBase<std::vector<std::array<double,3>>, std::vector<double>>;
    std::unique_ptr<base_quad_type> base_quad_default = nullptr;
    //std::unique_ptr<base_quad_type> base_quad_other   = nullptr;

    int default_npts, other_npts;
    int quad_num;
    SECTION("LebedevLaikov") {
      using quad_type = LebedevLaikov<double>;
      quad_num = INTXX_ANGQ_LEB;
      name = "LEBEDEV_LAIKOV";
      default_npts = 302;
      other_npts   = 974;
    }

    SECTION("Delley") {
      using quad_type = Delley<double>;
      quad_num = INTXX_ANGQ_DEL;
      name = "DELLEY";
      default_npts = 302;
      other_npts   = 974;
      base_quad_default = std::make_unique<quad_type>(default_npts);
      //base_quad_other = std::make_unique<quad_type>(other_npts);
    }

    SECTION("AB") {
      using quad_type = AhrensBeylkin<double>;
      quad_num = INTXX_ANGQ_AB;
      name = "AHRENS_BEYLKIN";
      default_npts = 372;
      other_npts   = 972;
      base_quad_default = std::make_unique<quad_type>(default_npts);
      //base_quad_other = std::make_unique<quad_type>(other_npts);
    }

    SECTION("WOM") {
      using quad_type = Womersley<double>;
      quad_num = INTXX_ANGQ_WOM;
      name = "WOMERSLEY";
      default_npts = 393;
      other_npts   = 969;
      base_quad_default = std::make_unique<quad_type>(default_npts);
      //base_quad_other = std::make_unique<quad_type>(other_npts);
    }

    // Initialize
    error = intxx_quad_init(&quad, quad_num);
    REQUIRE(error == INTXX_SUCCESS);

    // Meta data
    REQUIRE(quad.info != NULL);
    REQUIRE(quad.npoints == -1);
    REQUIRE(quad._state_quad == NULL);
    REQUIRE(quad.info->number == quad_num);
    REQUIRE(quad.info->kind == INTXX_ANG_QUAD);
    REQUIRE(quad.info->dim == 3);
    REQUIRE(quad.info->generate != NULL);
    REQUIRE(quad.info->destroy != NULL);
    REQUIRE(!strcmp(quad.info->name, name));


    // Check default parameters
    REQUIRE(quad.info->ext_params.n == 0);

    // Get before set
    int npts;
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_INVALID_OUT);
    REQUIRE(npts == -1);

    // Set NPTS (incorrect)
    error = intxx_quad_set_npts(&quad, default_npts+1);
    REQUIRE(error == INTXX_INVALID_ARG);
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_INVALID_OUT);
    REQUIRE(npts == -1);

    // Set NPTS (correct)
    error = intxx_quad_set_npts(&quad, default_npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(quad.npoints == default_npts);

    // Get NPTS
    error = intxx_quad_get_npts(&quad, &npts);
    REQUIRE(error == INTXX_SUCCESS);
    REQUIRE(npts == default_npts);

    // Check Quadrature Generation and Destruction
    if(base_quad_default) {
      intxx_generate_quad(&quad);
      REQUIRE(quad._state_quad != NULL);

      // Check validity of the state
      auto state_as_quad = reinterpret_cast<base_quad_type*>(quad._state_quad);
      REQUIRE(state_as_quad->npts() == npts);
      for(auto i = 0; i < npts; ++ i) {
        REQUIRE(state_as_quad->points()[i] == base_quad_default->points()[i]);
        REQUIRE(state_as_quad->weights()[i] == base_quad_default->weights()[i]);
      }

      intxx_destroy_quad(&quad);
      REQUIRE(quad._state_quad == NULL);
    }

  }

  intxx_quad_end(&quad);

}
