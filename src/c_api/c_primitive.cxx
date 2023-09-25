#include <integratorxx/quadratures/primitive.hpp>

#include "c_internal.h"
#include "c_api_util.hpp"

extern "C" {

// C++ types for internal quadrature states
using uniform_quad_type     = IntegratorXX::UniformTrapezoid<double,double>;
using gaussleg_quad_type    = IntegratorXX::GaussLegendre<double,double>;
using gausslob_quad_type    = IntegratorXX::GaussLobatto<double,double>;
using gausscheb1_quad_type  = IntegratorXX::GaussChebyshev1<double,double>;
using gausscheb2_quad_type  = IntegratorXX::GaussChebyshev2<double,double>;
using gausscheb3_quad_type  = IntegratorXX::GaussChebyshev3<double,double>;
using gausscheb2m_quad_type = IntegratorXX::GaussChebyshev2Modified<double,double>;

/*************************************************************/
/****** Forward declare C-API for Primitive Quadratures ******/
/*************************************************************/

#define FWD_PRM(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p); \
int intxx_generate_##cname(intxx_quad_type* p); \
int intxx_destroy_##cname(intxx_quad_type* p); 

FWD_PRM(uniform)
FWD_PRM(gaussleg)
FWD_PRM(gausslob)
FWD_PRM(gausscheb1)
FWD_PRM(gausscheb2)
FWD_PRM(gausscheb2m)
FWD_PRM(gausscheb3)

#undef FWD_PRM

/*************************************************************/
/****** Runtime Generator for Primitive Quadrature Info ******/
/*************************************************************/

int intxx_get_prmq_info(intxx_quad_info_type* p, int quad) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  quad = quad & INTXX_PRMQ_MASK;
  p->number = quad;
  p->kind   = INTXX_PRM_QUAD;
  switch(quad) {
    case INTXX_PRMQ_UNIFORM:
      p->name = "UNIFORM";
      return intxx_get_uniform_info(p);
    case INTXX_PRMQ_GAUSSLEG:
      p->name = "GAUSS_LEGENDRE";
      return intxx_get_gaussleg_info(p);
    case INTXX_PRMQ_GAUSSCHEB_1:
      p->name = "GAUSS_CHEBYSHEV_1";
      return intxx_get_gausscheb1_info(p);
    case INTXX_PRMQ_GAUSSCHEB_2:
      p->name = "GAUSS_CHEBYSHEV_2";
      return intxx_get_gausscheb2_info(p);
    case INTXX_PRMQ_GAUSSCHEB_2MOD:
      p->name = "GAUSS_CHEBYSHEV_2MOD";
      return intxx_get_gausscheb2m_info(p);
    case INTXX_PRMQ_GAUSSCHEB_3:
      p->name = "GAUSS_CHEBYSHEV_3";
      return intxx_get_gausscheb3_info(p);
    case INTXX_PRMQ_GAUSSLOB:
      p->name = "GAUSS_LOBATTO";
      return intxx_get_gausslob_info(p);
    default:
      return INTXX_INVALID_QUAD;
  }
}

/*******************************************************/
/****** Info generation for Primitive Quadratures ******/
/*******************************************************/

#define INTXX_NOPARAM_GET_INFO_IMPL(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p) { \
  return intxx_noparam_info(p, &intxx_generate_##cname, \
    &intxx_destroy_##cname);                            \
}

INTXX_NOPARAM_GET_INFO_IMPL(uniform);
INTXX_NOPARAM_GET_INFO_IMPL(gaussleg);
INTXX_NOPARAM_GET_INFO_IMPL(gausslob);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb1);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb2);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb2m);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb3);

#undef INTXX_NOPARAM_GET_INFO_IMPL

/*************************************************************/
/****** Allocate/Deallocate Primitive Quadrature States ******/
/*************************************************************/

#define INTXX_GD_BASIC_IMPL(cname, qname, ...)   \
int intxx_generate_##cname(intxx_quad_type* p) { \
  return intxx_generate_noparam_impl<qname>(p, ##__VA_ARGS__);  \
}                                                \
int intxx_destroy_##cname(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);   \
}

INTXX_GD_BASIC_IMPL(uniform,      uniform_quad_type, 0.0, 1.0);
INTXX_GD_BASIC_IMPL(gaussleg,      gaussleg_quad_type     );
INTXX_GD_BASIC_IMPL(gausslob,      gausslob_quad_type     );
INTXX_GD_BASIC_IMPL(gausscheb1,    gausscheb1_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb2,    gausscheb2_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb3,    gausscheb3_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb2m, gausscheb2m_quad_type);

#undef INTXX_GD_BASIC_IMPL

} // extern C
