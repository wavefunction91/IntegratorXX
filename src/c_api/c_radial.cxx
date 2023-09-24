#include "c_internal.h"
#include <cstddef>

extern "C" {

int intxx_get_uniform_info(intxx_quad_info_type* p);
int intxx_get_gaussleg_info(intxx_quad_info_type* p);
int intxx_get_gausslob_info(intxx_quad_info_type* p);
int intxx_get_gausscheb1_info(intxx_quad_info_type* p);
int intxx_get_gausscheb2_info(intxx_quad_info_type* p);
int intxx_get_gausscheb2m_info(intxx_quad_info_type* p);
int intxx_get_gausscheb3_info(intxx_quad_info_type* p);

void intxx_default_quad_info(intxx_quad_info_type* p) {
  p->number = 0; // Invalid
  p->kind   = 0; // Invalid
  p->name   = "Default"; // Default state
  //p->init   = NULL;
}

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

int intxx_noparam_info(intxx_quad_info_type* p) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  p->ext_params.n = 0; /// No External Parameters

  return INTXX_SUCCESS;
}

int intxx_get_uniform_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gausslob_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gaussleg_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gausscheb1_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gausscheb2_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gausscheb2m_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}
int intxx_get_gausscheb3_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p);
}

}
